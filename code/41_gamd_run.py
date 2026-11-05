#!/usr/bin/env python3
"""41_gamd_run.py — 在「修好的受体」上部署 GaMD（§3.9，Figure 5）。

对象：CHEMBL848（主，§3.8 唯一三副本一致的候选）、CHEMBL1129（对照）、
LxxLL_peptide（阳性对照参考）。起点取 batch 2 各副本 20 ns cMD 的**末态**
`prod.gro` —— 那套体系就是锌位点已修复的受体（Zn–S 0.239–0.241 nm），
因此无需重新平衡体系本身，只需按 GaMD 规范重跑自适应阶段。

二进制：`/usr/local/gromacs-gamd`（GROMACS 2025.4 + GaMD v2.1，由
`code/40_gamd_build.py` 构建）。原版 /opt/gmx2025.4 **不含** GaMD。

关键约束（fork README，违反会在 dynamics 之前直接被拒）：
  - mdp 必须 `integrator = md`、`gamd = yes`、`nstcalcenergy = 1`、`nstfout = 0`
  - NPT 只能用 C-rescale（isotropic / semiisotropic）
  - GPU 路径必须 `-ntmpi 1 -nb gpu -pme gpu -bonded gpu -update gpu`
  - `gamd.in` 由 mdrun 从**当前工作目录**读取，grompp 不把它写进 tpr
  - 每个 replica 必须用独立目录（GaMD 输出文件名是固定的）

用法：
    python3 code/41_gamd_run.py prep              # 建目录 + grompp
    python3 code/41_gamd_run.py smoke             # 2 ns 冒烟，量 ns/day
    python3 code/41_gamd_run.py launch [NSTEPS]   # 正式跑（默认 52e6 = 104 ns）
    python3 code/41_gamd_run.py status
    python3 code/41_gamd_run.py pull              # 拉回小型分析结果
"""
from __future__ import annotations

import os
import sys
import time

import paramiko

# 云端实例停机/重启后 SSH 端口会轮换：用环境变量覆盖，别再改源码
HOST = os.environ.get("HPV_SSH_HOST", "hz-4.matpool.com")
PORT = int(os.environ.get("HPV_SSH_PORT", "26255"))
USER = os.environ.get("HPV_SSH_USER", "root")
PW = os.environ.get("HPV_SSH_PW", "wmKG=3YRsZaHmij0")
BASE = "/root/autodl-tmp/HPV_structure/candidates_run"
GBASE = f"{BASE}/gamd"
GMXRC = "/usr/local/gromacs-gamd/bin/GMXRC"

# (体系, 副本) —— 与 batch 2 一致的三副本设计
SYSTEMS = ["CHEMBL848_pose1", "CHEMBL1129_pose1", "LxxLL_peptide"]
REPLICAS = [1, 2, 3]

# 每个副本的 ld_seed：显式写死，保证可复现（与 batch 2 的 seed 无关，
# 因为这里是从末态结构重新赋速，目的是三个独立副本）。
SEEDS = {1: 49411, 2: 49422, 3: 49433}

# 自适应阶段（步）——ntave 的整数倍，且 ntcmdprep < ntcmd、ntebprep < nteb
NTAVEGAMD = 50000
NTCMDPREP = 100000      # 0.2 ns  常规 MD 准备
NTCMD = 1000000         # 2.0 ns  常规 MD 统计
NTEBPREP = 200000       # 0.4 ns  boost 准备
NTEB = 1000000          # 2.0 ns  boost 统计 + 自适应参数
# 生产阶段从 step > NTCMD + NTEB = 2e6（4 ns）开始

NSTX = 5000             # 轨迹每 10 ps 一帧；reweight_nst 与之对齐
REWEIGHT_NST = 5000

MDP = """integrator               = md
dt                       = 0.002
nsteps                   = {nsteps}
nstxout-compressed       = {nstx}
nstenergy                = {nstx}
nstlog                   = {nstx}
nstcalcenergy            = 1
nstfout                  = 0
gamd                     = yes
nstlist                  = 10
cutoff-scheme            = Verlet
rlist                    = 1.0
rcoulomb                 = 1.0
rvdw                     = 1.0
coulombtype              = PME
tcoupl                   = v-rescale
tc_grps                  = System
tau_t                    = 0.1
ref_t                    = 300
pcoupl                   = C-rescale
pcoupltype               = isotropic
nstpcouple               = 100
tau_p                    = 5.0
compressibility          = 4.5e-5
ref_p                    = 1.0
constraints              = h-bonds
; 从 batch 2 末态结构重新赋速 -> 必须 continuation = no，
; 否则 grompp 报 "Generating velocities is inconsistent with attempting to
; continue a previous run"。代价是 grompp 会对起始坐标再约束一次，这正是我们要的。
continuation             = no
gen_vel                  = yes
gen_temp                 = 300
ld_seed                  = {seed}
"""

GAMD_IN = f"""igamd 3
irest_gamd 0
iE 1
iEP 1
iED 1
ntcmdprep {NTCMDPREP}
ntcmd {NTCMD}
ntebprep {NTEBPREP}
nteb {NTEB}
ntave {NTAVEGAMD}
reweight_nst {REWEIGHT_NST}
para_nst {REWEIGHT_NST}
sigma0P 6.0
sigma0D 6.0
"""

MDRUN = ("gmx mdrun -s gamd.tpr -deffnm gamd -ntmpi 1 -ntomp 8 "
         "-pin on -pinoffset {pinoff} -nb gpu -pme gpu -bonded gpu "
         "-update gpu -gpu_id 0")


def cli() -> paramiko.SSHClient:
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, port=PORT, username=USER, password=PW, timeout=30)
    return c


def sh(c, cmd, timeout=180):
    _, o, e = c.exec_command(cmd, timeout=timeout)
    return o.read().decode(errors="replace").strip(), e.read().decode(errors="replace").strip()


def run_bg(c, cmd):
    """后台执行；绝不读该 channel（paramiko 会阻塞）。
    用 bash -c 包一层：这样 cmd 自带的 `> log 2>&1` 在子 shell 内生效，
    不会被外层 `> /dev/null` 覆盖（之前直接拼接导致 runner 日志全空）。"""
    c.exec_command(f"setsid nohup bash -c '{cmd}' < /dev/null > /dev/null 2>&1 &", timeout=10)


def p(*a):
    print(*a, flush=True)


def jobs():
    return [(s, r) for s in SYSTEMS for r in REPLICAS]


def do_prep(c, nsteps=52_000_000):
    mdp = MDP.format(nsteps=nsteps, nstx=NSTX, seed="{seed}")
    sftp = c.open_sftp()

    def put(path, text):
        with sftp.open(path, "w") as fh:
            fh.write(text)

    for sysname, rep in jobs():
        d = f"{GBASE}/{sysname}_rep{rep}"
        sh(c, f"mkdir -p {d}")
        sh(c, f"cp {BASE}/{sysname}/rep{rep}/prod.gro {d}/start.gro && "
               f"cp {BASE}/{sysname}/rep{rep}/topol.top {d}/topol.top && "
               f"ls -la {d}/start.gro {d}/topol.top", 300)
        put(f"{d}/gamd.mdp", mdp.format(seed=SEEDS[rep]))
        put(f"{d}/gamd.in", GAMD_IN)
        out, err = sh(c, f"source {GMXRC} > /dev/null 2>&1; cd {d} && "
                         f"gmx grompp -f gamd.mdp -c start.gro -p topol.top "
                         f"-o gamd.tpr -maxwarn 2 > grompp.log 2>&1; "
                         f"echo rc=$?; tail -6 grompp.log", 900)
        p(f"{sysname}_rep{rep}: {out.splitlines()[-1] if out else ''}")
        if "rc=0" not in out:
            p("  !! grompp 失败：", out[-800:])
            return 1
    sftp.close()
    p("prep 完成：", len(jobs()), "个体系")
    return 0


def runner_script(gpu, todo, nsteps):
    """每张 GPU 一个串行 runner。"""
    lines = ["#!/bin/bash", f"source {GMXRC} > /dev/null 2>&1", "set -x"]
    for sysname, rep in todo:
        d = f"{GBASE}/{sysname}_rep{rep}"
        lines += [
            f"echo '===== {sysname}_rep{rep} on GPU{gpu} ====='",
            f"cd {d}",
            f"export CUDA_VISIBLE_DEVICES={gpu}",
            MDRUN.format(pinoff=8 * gpu) + f" >> {d}/mdrun.log 2>&1",
            f"echo 'rc='$? >> {d}/mdrun.log",
        ]
    lines.append(f"echo ALLDONE_GPU{gpu}")
    return "\n".join(lines) + "\n"


def clean(c):
    """清掉冒烟/上一轮残留。GaMD 输出名固定，mdrun 会拒绝覆盖而生成 #backup#。"""
    for sysname, rep in jobs():
        d = f"{GBASE}/{sysname}_rep{rep}"
        sh(c, f"cd {d} 2>/dev/null && rm -f gamd.xtc gamd.log gamd.edr gamd.cpt "
               f"gamd.gro gamd-reweight.dat gamd-para.dat gamd-restart.dat "
               f"gamd-restart-prev.dat gamd-restart.dat.tmp mdrun.log "
               f"mdrun_smoke.log '#gamd.xtc.1#' 2>/dev/null; true", 120)
    p("已清理旧输出")


def launch(c, nsteps, tag="gamd"):
    allj = jobs()
    # 9 个作业 / 4 张 GPU，round-robin 使最长队列 = ceil(9/4) = 3
    buckets = {g: [] for g in range(4)}
    for i, j in enumerate(allj):
        buckets[i % 4].append(j)
    clean(c)
    do_prep(c, nsteps=nsteps)
    sftp = c.open_sftp()
    for gpu, todo in buckets.items():
        if not todo:
            continue
        script = runner_script(gpu, todo, nsteps)
        rp = f"{GBASE}/runner_gpu{gpu}_{tag}.sh"
        with sftp.open(rp, "w") as fh:
            fh.write(script)
        sh(c, f"chmod +x {rp}")
        run_bg(c, f"bash {rp} > {GBASE}/runner_gpu{gpu}_{tag}.log 2>&1")
        p(f"GPU{gpu}: {todo}")
    sftp.close()
    p(f"已启动，nsteps={nsteps}（{nsteps * 0.002 / 1000:.0f} ns）")


def gamd_in(igamd, s0p, s0d, irest=0):
    """生成 gamd.in（fork 解析器很严格：多余键/未知键/尾随字符都会被拒）。"""
    lines = [f"igamd {igamd}", f"irest_gamd {irest}", "iE 1", "iEP 1", "iED 1",
             f"ntcmdprep {NTCMDPREP}", f"ntcmd {NTCMD}",
             f"ntebprep {NTEBPREP}", f"nteb {NTEB}", f"ntave {NTAVEGAMD}",
             f"reweight_nst {REWEIGHT_NST}", f"para_nst {REWEIGHT_NST}"]
    if igamd in (1, 3):
        lines.append(f"sigma0P {s0p}")
    if igamd in (2, 3):
        lines.append(f"sigma0D {s0d}")
    return "\n".join(lines) + "\n"


def do_pilot2(c, args):
    """单个 boost 标定 pilot：python3 41 pilot2 <tag> <igamd> <s0P> <s0D> <gpu>"""
    tag, ig, s0p, s0d, gpu = args[0], int(args[1]), float(args[2]), float(args[3]), int(args[4])
    nsteps = 6_000_000
    d = f"{GBASE}/pilot_{tag}"
    sh(c, f"rm -rf {d}; mkdir -p {d} && cp {BASE}/CHEMBL848_pose1/rep1/prod.gro {d}/start.gro "
           f"&& cp {BASE}/CHEMBL848_pose1/rep1/topol.top {d}/topol.top", 300)
    sftp = c.open_sftp()
    with sftp.open(f"{d}/gamd.mdp", "w") as fh:
        fh.write(MDP.format(nsteps=nsteps, nstx=NSTX, seed=SEEDS[1]))
    with sftp.open(f"{d}/gamd.in", "w") as fh:
        fh.write(gamd_in(ig, s0p, s0d))
    ok, _ = sh(c, f"source {GMXRC} > /dev/null 2>&1; cd {d} && gmx grompp -f gamd.mdp "
                  f"-c start.gro -p topol.top -o gamd.tpr -maxwarn 2 > grompp.log 2>&1; echo rc=$?", 900)
    p(f"{tag}: grompp {ok.splitlines()[-1]}")
    rp = f"{GBASE}/runner_pilot_{tag}.sh"
    with sftp.open(rp, "w") as fh:
        fh.write("#!/bin/bash\n"
                 f"source {GMXRC} > /dev/null 2>&1\n"
                 f"cd {d}\nexport CUDA_VISIBLE_DEVICES={gpu}\n"
                 + MDRUN.format(pinoff=8 * gpu) + f" > {d}/mdrun.log 2>&1\n"
                 f"echo rc=$? >> {d}/mdrun.log\n")
    sftp.close()
    sh(c, f"chmod +x {rp}")
    run_bg(c, f"bash {rp} > {GBASE}/runner_pilot_{tag}.log 2>&1")
    p(f"pilot {tag}: igamd={ig} sigma0P={s0p} sigma0D={s0d} → GPU{gpu}")
    return 0


def do_pilot(c):
    """boost 标定 pilot：实测不同 (igamd, sigma0) 下的 ΔV 分布。

    为什么需要：正式 9 条 100 ns 跑完后发现 σ_ΔV ≈ 6.5 kcal/mol，
    β·σ ≈ 11 —— 远超累积量重加权可用区间（βσ ≲ 1），
    二阶累积量与指数平均给出互相矛盾的 PMF（20–34 vs 3–16 kcal/mol），
    说明 boost 太大、PMF 不可用。必须先标定到可重加权的 boost 再重跑。
    """
    cfgs = [("dual_s1.0", 3, 1.0, 1.0), ("pot_s2.0", 1, 2.0, 0.0),
            ("dual_s2.0", 3, 2.0, 2.0)]
    nsteps = 6_000_000          # 12 ns：4 ns 自适应 + 8 ns 生产
    sftp = c.open_sftp()
    for k, (tag, ig, s0p, s0d) in enumerate(cfgs):
        gpu = 1 + k             # GPU0 留给还在跑的 LxxLL_rep3
        d = f"{GBASE}/pilot_{tag}"
        sh(c, f"mkdir -p {d} && cp {BASE}/CHEMBL848_pose1/rep1/prod.gro {d}/start.gro "
               f"&& cp {BASE}/CHEMBL848_pose1/rep1/topol.top {d}/topol.top", 300)
        with sftp.open(f"{d}/gamd.mdp", "w") as fh:
            fh.write(MDP.format(nsteps=nsteps, nstx=NSTX, seed=SEEDS[1]))
        with sftp.open(f"{d}/gamd.in", "w") as fh:
            fh.write(gamd_in(ig, s0p, s0d))
        sh(c, f"source {GMXRC} > /dev/null 2>&1; cd {d} && gmx grompp -f gamd.mdp "
               f"-c start.gro -p topol.top -o gamd.tpr -maxwarn 2 > grompp.log 2>&1; "
               f"echo rc=$?", 900)
        rp = f"{GBASE}/runner_pilot_{tag}.sh"
        with sftp.open(rp, "w") as fh:
            fh.write("#!/bin/bash\n"
                     f"source {GMXRC} > /dev/null 2>&1\n"
                     f"cd {d}\nexport CUDA_VISIBLE_DEVICES={gpu}\n"
                     + MDRUN.format(pinoff=8 * gpu) + f" > {d}/mdrun.log 2>&1\n"
                     f"echo rc=$? >> {d}/mdrun.log\n")
        sh(c, f"chmod +x {rp}")
        run_bg(c, f"bash {rp} > {GBASE}/runner_pilot_{tag}.log 2>&1")
        p(f"pilot {tag}: igamd={ig} sigma0P={s0p} sigma0D={s0d} → GPU{gpu}")
    sftp.close()
    p("3 个 pilot 已启动（12 ns 各），约 40–60 min")
    return 0


def do_smoke(c):
    """2 ns 冒烟：验证 GPU-resident 路径 + 量 ns/day。"""
    nsteps = 1_100_000  # 2.2 ns
    ok = do_prep(c, nsteps=nsteps)
    if ok:
        return 1
    d = f"{GBASE}/{SYSTEMS[0]}_rep1"
    sftp = c.open_sftp()
    rp = f"{GBASE}/runner_smoke.sh"
    with sftp.open(rp, "w") as fh:
        fh.write("#!/bin/bash\n"
                 f"source {GMXRC} > /dev/null 2>&1\n"
                 f"cd {d}\nexport CUDA_VISIBLE_DEVICES=0\n"
                 + MDRUN.format(pinoff=0) + f" > {d}/mdrun_smoke.log 2>&1\n"
                 f"echo rc=$? >> {d}/mdrun_smoke.log\n")
    sftp.close()
    sh(c, f"chmod +x {rp}")
    run_bg(c, f"bash {rp} > {GBASE}/runner_smoke.log 2>&1")
    p("冒烟已启动（GPU0）")
    for k in range(80):
        time.sleep(30)
        out, _ = sh(c, f"tail -3 {d}/mdrun_smoke.log 2>/dev/null | cut -c1-120; "
                       f"grep -c 'rc=' {d}/mdrun_smoke.log 2>/dev/null", 60)
        p(f"  [{(k + 1) * 30}s] {out.strip()}")
        if out.strip().endswith("1"):
            break
    out, _ = sh(c, f"grep -E 'Performance|ns/day|GaMD execution mode|CUDA Graph' "
                   f"{d}/mdrun_smoke.log | head -8", 60)
    p(out)
    return 0


def do_status(c):
    out, _ = sh(c, f"ls {GBASE} 2>&1 | head -30; echo '--- running ---'; "
                   f"pgrep -af 'mdrun' | head -8", 60)
    p(out)
    p("-" * 50)
    for sysname, rep in jobs():
        d = f"{GBASE}/{sysname}_rep{rep}"
        out, _ = sh(c, f"cd {d} 2>/dev/null && {{ "
                       f"grep -oE 'Step *[0-9]+' gamd.log 2>/dev/null | tail -1; "
                       f"grep -E 'Performance' gamd.log 2>/dev/null | tail -1; "
                       f"ls -la gamd.xtc 2>/dev/null | awk '{{print $5}}'; "
                       f"tail -1 gamd.log 2>/dev/null | cut -c1-100; }} "
                       f"|| echo 'no dir'", 60)
        p(f"{sysname}_rep{rep}: {out.strip()}")
    out, _ = sh(c, "nvidia-smi --query-gpu=index,utilization.gpu,memory.used "
                   "--format=csv,noheader", 60)
    p(out)
    return 0


def do_launch2(c, args):
    """第二版 GaMD：用 pilot 标定出的可重加权 boost 参数重跑。

    用法：python3 code/41_gamd_run.py launch2 <igamd> <sigma0P> <sigma0D> [nsteps] [gpus]
    例：  python3 code/41_gamd_run.py launch2 3 1.0 1.0 52000000 1,2,3
    """
    ig = int(args[0]); s0p = float(args[1]); s0d = float(args[2])
    nsteps = int(args[3]) if len(args) > 3 else 52_000_000
    gpus = [int(x) for x in args[4].split(",")] if len(args) > 4 else [1, 2, 3]
    GBASE2 = f"{BASE}/gamd_v2"
    allj = jobs()
    buckets = {g: [] for g in gpus}
    for i, j in enumerate(allj):
        buckets[gpus[i % len(gpus)]].append(j)
    sftp = c.open_sftp()
    for gpu, todo in buckets.items():
        lines = ["#!/bin/bash", f"source {GMXRC} > /dev/null 2>&1", "set -x"]
        for sysname, rep in todo:
            d = f"{GBASE2}/{sysname}_rep{rep}"
            sh(c, f"mkdir -p {d} && cp {BASE}/{sysname}/rep{rep}/prod.gro {d}/start.gro "
                   f"&& cp {BASE}/{sysname}/rep{rep}/topol.top {d}/topol.top", 300)
            with sftp.open(f"{d}/gamd.mdp", "w") as fh:
                fh.write(MDP.format(nsteps=nsteps, nstx=NSTX, seed=SEEDS[rep]))
            with sftp.open(f"{d}/gamd.in", "w") as fh:
                fh.write(gamd_in(ig, s0p, s0d))
            ok, _ = sh(c, f"source {GMXRC} > /dev/null 2>&1; cd {d} && gmx grompp "
                          f"-f gamd.mdp -c start.gro -p topol.top -o gamd.tpr "
                          f"-maxwarn 2 > grompp.log 2>&1; echo rc=$?", 900)
            p(f"{sysname}_rep{rep}: {ok.splitlines()[-1]}")
            lines += [f"echo '===== {sysname}_rep{rep} GPU{gpu} ====='",
                      f"cd {d}", f"export CUDA_VISIBLE_DEVICES={gpu}",
                      MDRUN.format(pinoff=8 * gpu) + f" >> {d}/mdrun.log 2>&1",
                      f"echo 'rc='$? >> {d}/mdrun.log"]
        lines.append(f"echo ALLDONE_GPU{gpu}")
        rp = f"{GBASE2}/runner_gpu{gpu}_v2.sh"
        sh(c, f"mkdir -p {GBASE2}")
        with sftp.open(rp, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        sh(c, f"chmod +x {rp}")
        run_bg(c, f"bash {rp} > {GBASE2}/runner_gpu{gpu}_v2.log 2>&1")
        p(f"GPU{gpu}: {[f'{a}_rep{b}' for a, b in todo]}")
    sftp.close()
    p(f"第二版已启动：igamd={ig} sigma0P={s0p} sigma0D={s0d} nsteps={nsteps} "
      f"（{nsteps * 0.002 / 1000:.0f} ns）→ {GBASE2}")
    return 0


def main() -> int:
    cmd = sys.argv[1] if len(sys.argv) > 1 else "status"
    c = cli()
    try:
        if cmd == "launch2":
            return do_launch2(c, sys.argv[2:])
        if cmd == "prep":
            return do_prep(c)
        if cmd == "smoke":
            return do_smoke(c)
        if cmd == "pilot":
            return do_pilot(c)
        if cmd == "pilot2":
            return do_pilot2(c, sys.argv[2:])
        if cmd == "launch":
            nsteps = int(sys.argv[2]) if len(sys.argv) > 2 else 52_000_000
            return launch(c, nsteps) or 0
        if cmd == "status":
            return do_status(c)
        print(__doc__)
        return 2
    finally:
        c.close()


if __name__ == "__main__":
    sys.exit(main())
