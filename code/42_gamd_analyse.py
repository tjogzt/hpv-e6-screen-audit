#!/usr/bin/env python3
"""42_gamd_analyse.py — GaMD 轨迹分析（远端 worker + 本地拉取）。

产出 §3.9 / Figure 5 需要的三个量：
  1. **配体 RMSD**（对蛋白 CA 做最优叠合后，配体重原子相对起始 pose 的 RMSD）
     —— 与 §3.8「pose retention」用的是同一个坐标，可直接对照。
  2. **最小重原子距离** d_min（配体 ↔ 受体）—— 与 §3.6 engagement 指标同源。
  3. **逐残基接触占据率**：每个受体残基在多少比例的帧里与配体重原子 ≤4.5 Å。

PBC 处理（关键，本机体系的盒子约 7–8 nm）：
  先用 `gmx trjconv -pbc mol` 把每个分子补全（蛋白与配体各自完整），
  再在 Python 里把配体整体按**最小镜像**平移到离蛋白最近的那个像，
  否则配体一旦跨周期边界，RMSD 会被算成几十 Å 的假值。

重加权不在这里做：boost 势 ΔV 来自 `gamd-reweight.dat`，需要与轨迹按 step
对齐后在本地做（见 `code/43_fig5_gamd.R` 前的 Python 重加权步骤）。

用法：
    python3 code/42_gamd_analyse.py deploy    # 上传 worker + trjconv 预处理
    python3 code/42_gamd_analyse.py run       # 对已完成作业跑分析（后台）
    python3 code/42_gamd_analyse.py status
    python3 code/42_gamd_analyse.py pull      # 拉回 CSV 到 results/gamd/
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import paramiko

# 云端实例停机/重启后 SSH 端口会轮换：用环境变量覆盖，别再改源码
HOST = os.environ.get("HPV_SSH_HOST", "hz-4.matpool.com")
PORT = int(os.environ.get("HPV_SSH_PORT", "26255"))
USER = os.environ.get("HPV_SSH_USER", "root")
PW = os.environ.get("HPV_SSH_PW", "wmKG=3YRsZaHmij0")
BASE = "/root/autodl-tmp/HPV_structure/candidates_run"
# GAMD_TAG 用于区分两轮 GaMD：gamd（σ0=6 强 boost 版）/ gamd_v2（标定后可重加权版）
TAG = os.environ.get("GAMD_TAG", "gamd").strip() or "gamd"
GBASE = f"{BASE}/{TAG}"
GMXRC = "/usr/local/gromacs-gamd/bin/GMXRC"

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
OUT = REBUILD / "results" / TAG

SYSTEMS = ["CHEMBL848_pose1", "CHEMBL1129_pose1", "LxxLL_peptide"]
REPLICAS = [1, 2, 3]
# 生产起点：按远端 gamd.in 实测的四个自适应段累计
#   ntcmdprep 100000 + ntcmd 1000000 + ntebprep 200000 + nteb 1000000
#   = 2,300,000 步 × 0.002 ps = 4600 ps
# （旧值 4000 是早期按 ntcmd+nteb 两段估的，把 0.6 ns 的 nteb 平衡算进了生产）
PROD_START_PS = 4600.0

# 界面残基（4GIZ 结构编号）→ gro resid 需 −1（见 44 文件头的编号映射说明）。
# 部署时由 do_deploy 读 results/29_interface_positions.csv 后注入 worker，
# 避免在 worker 里二次硬编码生物学信息。
F_IF = REBUILD / "results" / "29_interface_positions.csv"

WORKER = r'''#!/usr/bin/env python3
"""GaMD 单作业分析 worker（在远端 gamd 环境里跑）。

为什么不用 mda.Universe(tpr)：MDAnalysis 2.9 的 TPRParser 不认 GROMACS 2025
的 tpx（version 138）。改用：
    - 拓扑：解析 start.gro（含 resname/resid/atomname）
    ⚠️ **参考构象不用 start.gro**：-pbc mol 处理后的轨迹与原始 gro 的
    分子可能落在**不同的周期像**上（实测 LxxLL 相差整整一个盒子向量，
    RMSD 直接被算成 63 Å）。改为把轨迹**第 0 帧**作为参考，
    单位与像约定都天然一致，frame0 的 RMSD 恒为 0（可作自检）。
  - 分子边界：解析 topol.top 的 [ moleculetype ] 原子数 + [ molecules ] 顺序
    （LxxLL 阳性对照的"配体"是 system2，残基名全是普通氨基酸，**不能**按
     resname LIG 找）。解析完必须用「分子总原子数 == .gro 原子数」断言，
    否则差一个原子就会把蛋白的最后一个原子算进配体，
    最小距离会假到 0.6 Å 量级（本机踩过）。
  - 轨迹：MDAnalysis.coordinates.XTC.XTCReader 独立使用（不需要 Universe）
"""
import sys, os, time
import numpy as np
from MDAnalysis.coordinates.XTC import XTCReader
from MDAnalysis.lib.distances import distance_array

job = sys.argv[1]
d = f"/root/autodl-tmp/HPV_structure/candidates_run/__TAG__/{job}"
CONTACT_CUT = 4.5          # Å，重原子接触判据


def read_gro(path):
    lines = open(path).read().splitlines()
    n = int(lines[1])
    resids, resnames, names, pos = [], [], [], []
    for i in range(2, 2 + n):
        L = lines[i]
        resids.append(int(L[0:5]))
        resnames.append(L[5:10].strip())
        names.append(L[10:15].strip())
        pos.append((float(L[20:28]), float(L[28:36]), float(L[36:44])))
    box = [float(v) for v in lines[2 + n].split()]
    return (np.array(resids), np.array(resnames), np.array(names),
            np.array(pos, dtype=np.float64), box, n)


def molecule_ranges(topol):
    """返回 [(molname, start, stop)]（按 [ molecules ] 的排列顺序）。"""
    lines = open(topol).read().splitlines()
    counts = {}                       # moleculetype -> 原子数
    cur_mt, n_at, state, order = None, 0, None, []
    for L in lines:
        s = L.split(";")[0].strip()
        if s.startswith("["):
            sec = s.strip("[] ").lower()
            if sec == "moleculetype":
                if cur_mt is not None:
                    counts.setdefault(cur_mt, n_at)
                cur_mt, n_at = None, 0
            elif sec == "molecules":
                if cur_mt is not None:
                    counts.setdefault(cur_mt, n_at)
                cur_mt, n_at = None, 0
            state = sec
            continue
        if not s:
            continue
        if state == "moleculetype" and cur_mt is None:
            cur_mt = s.split()[0]
        elif state == "atoms":
            # 只数以原子序号开头的行，避免把分节内的杂行算进去
            f = s.split()
            if f and f[0].lstrip("-").isdigit():
                n_at += 1
        elif state == "molecules":
            f = s.split()
            if len(f) >= 2 and f[1].isdigit():
                order.append((f[0], int(f[1])))
    if cur_mt is not None:
        counts.setdefault(cur_mt, n_at)
    ranges, k = [], 0
    for name, cnt in order:
        na = counts.get(name)
        if na is None:
            raise SystemExit(f"moleculetype {name} 的原子数未解析到")
        for _ in range(cnt):
            ranges.append((name, k, k + na))
            k += na
    return ranges


def kabsch(mob, ref):
    mc, rc = mob.mean(0), ref.mean(0)
    P, Q = mob - mc, ref - rc
    C = P.T @ Q
    V, S, Wt = np.linalg.svd(C)
    dd = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, dd]) @ Wt
    return R, mc, rc


gro = f"{d}/start.gro"
xtc = f"{d}/gamd_nj.xtc"                     # -pbc mol 后再 -pbc nojump
if not os.path.exists(xtc):
    xtc = f"{d}/gamd_mol.xtc"
if not os.path.exists(xtc):
    xtc = f"{d}/gamd.xtc"

resids, resnames, names, pos0, box0, n_atom = read_gro(gro)
pos0 = pos0 * 10.0                    # nm -> Å（与 XTCReader 的单位一致）
ranges = molecule_ranges(f"{d}/topol.top")
tot = ranges[-1][2] if ranges else -1
if tot != n_atom:
    raise SystemExit(f"分子边界解析不一致：topol 给出 {tot} 原子，gro 有 {n_atom}")
mol0, mol1 = ranges[0], ranges[1]
print(f"[{job}] molecules: " + ", ".join(f"{n}[{a}:{b}]" for n, a, b in ranges[:4])
      + f"  ... total={tot}/{n_atom} OK")

heavy = np.array([not nm.startswith("H") for nm in names])
prot_all = np.zeros(n_atom, dtype=bool)
prot_all[mol0[1]:mol0[2]] = True
lig_all = np.zeros(n_atom, dtype=bool)
lig_all[mol1[1]:mol1[2]] = True

prot_idx = np.where(prot_all & heavy)[0]
lig_idx = np.where(lig_all & heavy)[0]
ca_idx = np.where(prot_all & (names == "CA"))[0]

# 界面残基 CA（gro resid，部署时注入）：用于**界面局部坐标**下的配体 RMSD。
# 为什么必须加这套坐标：E6 是双结构域蛋白，LxxLL 界面所在结构域相对全蛋白
# 坐标系会整体漂移（实测 848 rep1 达 18.6 ± ? Å），此时「全蛋白 CA 叠合后的
# 配体 RMSD」把蛋白构象变化一并算进配体位移，系统性高估 4–5 Å。
IF_RESIDS = [__IF_RESIDS__]
if_idx = np.where(prot_all & (names == "CA") & np.isin(resids, IF_RESIDS))[0]
if if_idx.size < 12:
    raise SystemExit(f"[{job}] 界面 CA 只匹配到 {if_idx.size} 个，检查 resid 编号映射")

# 参考构象 = 轨迹第 0 帧（= batch 2 二十纳秒末态、GaMD 的起点），
# 与轨迹同一单位/同一周期像约定；见文件头说明。
ref_lig = None
ref_ca = None
ref_if = None

# 逐残基：只在蛋白范围内按 (resid) 分槽
pr = resids[prot_idx]
uniq, slot = np.unique(pr, return_inverse=True)
rn = {}
for i in prot_idx:
    rn.setdefault(int(resids[i]), resnames[i])
resnames_out = [rn.get(int(x), "UNK") for x in uniq]
n_res = len(uniq)

print(f"[{job}] ligand heavy={len(lig_idx)} protein heavy={len(prot_idx)} "
      f"CA={len(ca_idx)} residues={n_res} atoms={n_atom}")

reader = XTCReader(xtc)
rows = []
mat = []          # 逐帧接触位图（与 rows 同序），供 boost 重加权占据率用
occ = np.zeros(n_res, dtype=np.int64)
nprod = 0
t0 = time.time()
for ts in reader:
    tps = float(ts.time)
    # 参考帧必须最先处理：它在 t=0，早于生产阶段起点，
    # 因此生产过滤要放在参考捕获之后（否则 ref 永远是 None）。
    box = ts.dimensions
    X = ts.positions
    P = X[prot_idx]
    L = X[lig_idx].copy()

    # --- PBC 处理全部交给 trjconv 预处理（-pbc mol → -pbc nojump）---
    # 曾试过在 Python 里搬动配体，三种搬法都错：
    #   (a) 按质心差做最小镜像 → 把配体搬进蛋白内部（min_dist 假到 0.66 Å）；
    #   (b) 按「最近原子对」的镜像位移平移 → 该原子对帧间翻转，配体被放进
    #       不同的周期像，RMSD 凭空多 60–80 Å（出现 dmin=2.6 Å 却 RMSD=63 Å 的矛盾）；
    #   (c) 叠合后对质心差做最小镜像 → 配体一旦真的解离扩散就被折回盒内，
    #       位移被系统性低估（41% 的帧触发修正即为证据）。
    # 正解：坐标连续化（-pbc nojump），此后 RMSD 直接算、无需任何盒修正；
    # 距离/接触仍用 distance_array 的 box 参数（最小镜像，对近距接触恒正确）。
    D = distance_array(L, P, box=box)        # (n_lig, n_prot)，含 PBC

    if tps == 0.0:
        # 第 0 帧即为起始 pose，同时作为 RMSD 参考（单位与周期像都与轨迹一致）
        ref_lig = L.copy()
        ref_ca = X[ca_idx].copy()
        ref_if = X[if_idx].copy()
        print(f"[{job}] reference = trajectory frame 0 "
              f"(ligand {ref_lig.shape[0]} heavy atoms, interface CA {if_idx.size})",
              flush=True)
        continue
    if tps < PROD_START_PS:                  # 只统计生产阶段
        continue
    if ref_lig is None:                      # 保险：轨迹缺第 0 帧
        raise SystemExit(f"[{job}] 轨迹缺少 t=0 帧，无法确定参考构象")

    # --- (a) 全蛋白 CA 叠合后的配体 RMSD（保留，供跨体系/文献对照）---
    R, mc, rc = kabsch(X[ca_idx], ref_ca)
    Lal = (R @ (L - mc).T).T + rc
    # 坐标已连续化（-pbc nojump），直接算，不做任何整盒折回
    rmsd = float(np.sqrt(((Lal - ref_lig) ** 2).sum(1).mean()))

    # --- (b) **界面残基 CA 叠合**后的配体 RMSD（姿态保持的主指标）---
    R2, m2, r2 = kabsch(X[if_idx], ref_if)
    Lal2 = (R2 @ (L - m2).T).T + r2
    rmsd_if = float(np.sqrt(((Lal2 - ref_lig) ** 2).sum(1).mean()))

    # --- (c) 界面区相对全蛋白坐标系的漂移（蛋白构象变化，解释 a 与 b 的差）---
    IFnow = (R @ (X[if_idx] - mc).T).T + rc
    if_drift = float(np.sqrt(((IFnow - ref_if) ** 2).sum(1).mean()))

    # --- 最小重原子距离 + 逐残基接触 ---
    dmin = float(D.min())
    per_atom_min = D.min(axis=0)
    resmin = np.full(n_res, np.inf)
    np.minimum.at(resmin, slot, per_atom_min)
    hit = resmin <= CONTACT_CUT
    occ += hit

    rows.append((int(round(tps)), rmsd, rmsd_if, if_drift, dmin, int(hit.sum())))
    # 逐帧接触位图：resocc.csv 只给计数，**无法做 boost 重加权**；
    # 有了逐帧 0/1 才能在 43 里按 exp(beta*dV) 重算占据率（v2 的 beta*sigma 已到 2.3，
    # 未重加权的占据率作为平衡系综估计是有偏的）。
    mat.append("".join("1" if h else "0" for h in hit))
    nprod += 1
    if nprod == 1:
        print(f"  [{job}] frame0: rmsd={rmsd:.3f} dmin={dmin:.3f} "
              f"ncontact={int(hit.sum())}", flush=True)
    if nprod % 2000 == 0:
        print(f"  [{job}] {nprod} frames, {time.time()-t0:.0f}s", flush=True)

print(f"[{job}] production frames = {nprod}")
with open(f"{d}/frames.csv", "w") as fh:
    fh.write("time_ps,lig_rmsd,lig_rmsd_if,if_drift,min_dist,n_res_contact\n")
    for t, r, ri, dr, m, k in rows:
        fh.write(f"{t},{r:.6f},{ri:.6f},{dr:.6f},{m:.6f},{k}\n")
with open(f"{d}/resocc.csv", "w") as fh:
    fh.write("resid,resname,n_frames_contact,occupancy\n")
    for i, rid in enumerate(uniq):
        fh.write(f"{int(rid)},{resnames_out[i]},{int(occ[i])},"
                 f"{occ[i]/max(nprod,1):.6f}\n")
# rescontact.csv：逐帧接触位图（行序与 frames.csv 一致），供重加权用
with open(f"{d}/rescontact.csv", "w") as fh:
    fh.write("resid," + ",".join(str(int(x)) for x in uniq) + "\n")
    fh.write("resname," + ",".join(resnames_out) + "\n")
    for (t, r, ri, dr, m, k), bits in zip(rows, mat):
        fh.write(f"{t},{','.join(bits)}\n")
print(f"[{job}] DONE")
'''


def cli() -> paramiko.SSHClient:
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, port=PORT, username=USER, password=PW, timeout=30)
    return c


def sh(c, cmd, timeout=180):
    _, o, e = c.exec_command(cmd, timeout=timeout)
    return o.read().decode(errors="replace").strip(), e.read().decode(errors="replace").strip()


def run_bg(c, cmd):
    # bash -c 包一层：保住 cmd 自带的 > log 重定向（否则被外层 /dev/null 覆盖）
    c.exec_command(f"setsid nohup bash -c '{cmd}' < /dev/null > /dev/null 2>&1 &", timeout=10)


def p(*a):
    print(*a, flush=True)


def jobs():
    return [f"{s}_rep{r}" for s in SYSTEMS for r in REPLICAS]


def interface_resids():
    """从 results/29_interface_positions.csv 取 24 个 LxxLL 界面残基，转成 gro resid。

    4GIZ 结构编号 → gro resid 需 −1（构建体取 4GIZ 2–143 位）。
    生物学信息只在这一处（CSV）定义，worker 不重复硬编码。
    """
    import csv
    if not F_IF.exists():
        raise SystemExit(f"缺少界面残基表：{F_IF}")
    out = []
    with open(F_IF, newline="") as fh:
        for row in csv.DictReader(fh):
            v = str(row.get("e6ap_lxxll_interface", "")).strip().lower()
            if v in ("true", "1", "yes"):
                out.append(int(float(row["struct_num"])) - 1)
    if len(out) < 12:
        raise SystemExit(f"界面残基只解析到 {len(out)} 个，检查 {F_IF}")
    return sorted(out)


def do_deploy(c):
    ifr = interface_resids()
    sftp = c.open_sftp()
    with sftp.open(f"{GBASE}/gamd_analyse_worker.py", "w") as fh:
        fh.write(WORKER.replace("PROD_START_PS", str(PROD_START_PS))
                       .replace("__TAG__", TAG)
                       .replace("__IF_RESIDS__", ", ".join(str(i) for i in ifr)))
    sftp.close()
    out, _ = sh(c, f"source /root/miniconda3/etc/profile.d/conda.sh; conda activate gamd; "
                   f"python -c \"import MDAnalysis; print(MDAnalysis.__version__)\" 2>&1 | tail -1")
    p("MDAnalysis:", out)
    p(f"worker 已上传（界面残基 {len(ifr)} 个，生产起点 {PROD_START_PS/1000:.1f} ns）")
    return 0


def do_run(c, force=False):
    """对已有完整 xtc 的作业跑分析（后台，串行）。

    force=True：删掉旧分析产物与 trjconv 缓存，9 条全部重算
    （换过指标定义/生产起点、或怀疑旧产物预处理不一致时用）。
    """
    joblist = " ".join(jobs())
    script = ["#!/bin/bash",
              "source /root/miniconda3/etc/profile.d/conda.sh",
              "conda activate gamd",
              f"source {GMXRC} > /dev/null 2>&1",
              f"for j in {joblist}; do",
              f'  echo "===== $j ====="',
              f"  cd {GBASE}/$j || {{ echo '  目录不存在，跳过'; continue; }}",
              ]
    if force:
        script += ["  echo '  [force] 清除旧产物与 trjconv 缓存'",
                   "  rm -f frames.csv resocc.csv rescontact.csv analyse.log "
                   "gamd_mol.xtc gamd_nj.xtc gamd_mol.tmp.xtc gamd_nj.tmp.xtc trjconv.log"]
    script += [
              # 跳过条件必须验**帧数**：只判断文件存在曾被 40 帧的陈旧半成品骗过
              # （mdrun 重跑尚未产出 xtc 的窗口期留下的产物 + 后来出现的
              #  'Finished mdrun' 一起满足了存在性检查）。生产 100 ns / 10 ps = 10001 帧，
              # 取 9000 作下限（容忍 trjconv 边界丢帧）。
              "  nf=$(wc -l < frames.csv 2>/dev/null || echo 0)",
              "  if [ \"${nf:-0}\" -ge 9000 ] && [ -s rescontact.csv ] && "
              "  grep -q 'Finished mdrun' gamd.log 2>/dev/null; then "
              "  echo '  已分析且帧数完整，跳过'; continue; fi",
              # trjconv 两遍：先 -pbc mol 把分子补全，再 -pbc nojump 让坐标连续。
              # 防护1：xtc 比缓存新就重跑（防拿到半成品轨迹时留下的旧缓存）。
              # 防护2：先写 .tmp 再 mv（防进程被杀留下残缺缓存被当成完整）。
              "  if [ ! -f gamd_mol.xtc ] || [ gamd.xtc -nt gamd_mol.xtc ]; then "
              "  echo System | gmx trjconv -f gamd.xtc -s gamd.tpr -pbc mol "
              "  -o gamd_mol.tmp.xtc >> trjconv.log 2>&1 && "
              "  mv gamd_mol.tmp.xtc gamd_mol.xtc; fi",
              "  if [ ! -f gamd_nj.xtc ] || [ gamd_mol.xtc -nt gamd_nj.xtc ]; then "
              "  echo System | gmx trjconv -f gamd_mol.xtc -s gamd.tpr -pbc nojump "
              "  -o gamd_nj.tmp.xtc >> trjconv.log 2>&1 && "
              "  mv gamd_nj.tmp.xtc gamd_nj.xtc; fi",
              f"  python {GBASE}/gamd_analyse_worker.py $j >> analyse.log 2>&1",
              "done"]
    script.append("echo ANALYSE_ALLDONE")
    rp = f"{GBASE}/runner_analyse.sh"
    sftp = c.open_sftp()
    with sftp.open(rp, "w") as fh:
        fh.write("\n".join(script) + "\n")
    sftp.close()
    run_bg(c, f"bash {rp} > {GBASE}/runner_analyse.log 2>&1")
    p("分析已启动（含 trjconv -pbc mol 预处理）")
    return 0


def do_status(c):
    for j in jobs():
        out, _ = sh(c, f"cd {GBASE}/{j} 2>/dev/null && {{ "
                       f"tail -2 analyse.log 2>/dev/null | cut -c1-110; "
                       f"wc -l frames.csv 2>/dev/null | awk '{{print $1\" frames\"}}'; "
                       f"}} || echo 'no dir'", 60)
        p(f"{j}: {out.strip()}")
    out, _ = sh(c, f"tail -3 {GBASE}/runner_analyse.log 2>/dev/null", 60)
    p("runner:", out)
    return 0


def do_pull(c):
    OUT.mkdir(parents=True, exist_ok=True)
    sftp = c.open_sftp()
    n = 0
    for j in jobs():
        for f in ("frames.csv", "resocc.csv", "rescontact.csv",
                  "gamd-reweight.dat", "gamd-para.dat"):
            try:
                sftp.get(f"{GBASE}/{j}/{f}", str(OUT / f"{j}__{f}"))
                n += 1
            except FileNotFoundError:
                p(f"  缺失 {j}/{f}")
    sftp.close()
    p(f"已拉取 {n} 个文件 → {OUT}")
    return 0


def main() -> int:
    cmd = sys.argv[1] if len(sys.argv) > 1 else "status"
    c = cli()
    try:
        if cmd == "deploy":
            return do_deploy(c)
        if cmd == "run":
            return do_run(c)
        if cmd == "rerun":                    # 强制 9 条全部重算（清旧产物+缓存）
            return do_run(c, force=True)
        if cmd == "status":
            return do_status(c)
        if cmd == "pull":
            return do_pull(c)
        print(__doc__)
        return 2
    finally:
        c.close()


if __name__ == "__main__":
    sys.exit(main())
