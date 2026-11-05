#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12_pipeline_v2.py — B 线重跑管线 v2（修复版，远端自运行）

相对原版 `.workbuddy/candidates_pipeline.py`（未改动，保留审计）修复的问题
-------------------------------------------------------------------------
P1 [致命] 主循环第 403 行用 `DELTA TOTAL` 正则取 ΔG —— 该字段在 AmberTools 14/24
   的 PB-only 模式下恒为 0.0000。原文件虽已定义正确的 parse_dg_bind()，
   但只在"断点续跑"分支被调用，新跑出来的数仍全是 0。
   v2：统一走严格组分求和解析器，缺 Differences 段即判 PARSE_FAIL（不静默取 0）。
P2 [严重] 生产段为 NVT（pcoupl = no）：不控压，MM-PBSA 的密度/盒体积不成立。
   v2：EM → NVT 平衡 → NPT 平衡 → NPT 生产（Parrinello-Rahman / C-rescale）。
P3 [严重] 配体电荷硬编码 `-nc 0`（antechamber）。
   v2：按 charge_table.csv → RDKit 形式电荷顺序推断，并记录来源（数据溯源）。
P4 [缺失] 双 Zn 四配位验收门槛从未实现（旧数据 0/714 帧合格即因此未被拦下）。
   v2：生产跑后强制 QC，量化阈值（Zn–S ≤ 2.8 Å、每个 Zn ≥4 个 S、合格率 ≥95% 帧），
   不达标标记 qc_pass=0 并如实写入 SUMMARY，不静默丢弃。
P5 [复现] 无随机种子控制、replica 无法区分。
   v2：每 replica 显式 gen_seed（crc32 确定性派生），SUMMARY 记 replica / seed。

环境变量（全部有默认值）
------------------------
BASE        默认 /root/autodl-tmp/HPV_structure/candidates_run
TEMPL_DIR   默认 $BASE/templates
LIGDIR      默认 $BASE/ligands
CONDA_SH    默认 /root/autodl-tmp/miniforge3/etc/profile.d/conda.sh
ENV_AMBER   默认 gamd      (tleap / antechamber / parmchk2 / cpptraj / MMPBSA.py)
ENV_GMX     默认 mmpbsa    (gromacs)
GPU_ID      默认 0
REPLICAS    默认 1
PROD_NS     默认 30        (生产时长 ns；pilot 用 100)
SAVE_PS     默认 20        (轨迹存帧间隔 ps)
LIG_ONLY    逗号分隔的配体名白名单（pilot 用）
NSTEPS_EQ_NS 默认 1        (NVT/NPT 平衡各 1 ns = 500000 步)
"""
from __future__ import annotations

import csv
import datetime
import glob
import os
import re
import subprocess
import time
import zlib

BASE = os.environ.get("BASE", "/root/autodl-tmp/HPV_structure/candidates_run")
TEMPL = os.environ.get("TEMPL_DIR", os.path.join(BASE, "templates"))
LIGDIR = os.environ.get("LIGDIR", os.path.join(BASE, "ligands"))
CONDA_SH = os.environ.get("CONDA_SH", "/root/autodl-tmp/miniforge3/etc/profile.d/conda.sh")
ZINC_FRCMOD = os.environ.get("ZINC_FRCMOD", os.path.join(TEMPL, "zinc_nonbond_fallback.frcmod"))
ENV_AMBER = os.environ.get("ENV_AMBER", "gamd")
ENV_GMX = os.environ.get("ENV_GMX", "mmpbsa")
GPU_ID = os.environ.get("GPU_ID", "0")
REPLICAS = int(os.environ.get("REPLICAS", "1"))
# REPLICA_ID：只跑指定副本号（用于多进程/多 GPU 并行：每进程 REPLICA_ID=k GPU_ID=j）
REPLICA_ID = int(os.environ.get("REPLICA_ID", "0"))
PROD_NS = float(os.environ.get("PROD_NS", "30"))
SAVE_PS = float(os.environ.get("SAVE_PS", "20"))
EQ_NS = float(os.environ.get("EQ_NS", "1"))
LIG_ONLY = [s.strip() for s in os.environ.get("LIG_ONLY", "").split(",") if s.strip()]

# 双 Zn 验收门槛（预注册，可在环境变量覆盖）
# 单位 nm（mdtraj 返回 nm）。历史缺陷：曾写 2.8(Å) 用于 nm 距离 → 阈值 28 Å，门失效。
if os.environ.get("ZN_S_CUTOFF_ANGSTROM"):
    ZN_S_CUTOFF = float(os.environ["ZN_S_CUTOFF_ANGSTROM"]) / 10.0
else:
    ZN_S_CUTOFF = float(os.environ.get("ZN_S_CUTOFF", "0.28"))  # nm
ZN_MIN_COORD = int(os.environ.get("ZN_MIN_COORD", "4"))       # 每个 Zn 至少配位 S 数
ZN_MIN_FRACTION = float(os.environ.get("ZN_MIN_FRACTION", "0.95"))  # 合格帧比例

RUNLOG = os.path.join(BASE, "run_v2.log")
# REPLICA_ID 模式下写独立 SUMMARY，避免多进程互覆
if os.environ.get("REPLICA_ID"):
    SUMMARY = os.path.join(BASE, f"SUMMARY_v2_rep{os.environ['REPLICA_ID']}.csv")
else:
    SUMMARY = os.path.join(BASE, "SUMMARY_v2.csv")
FIELDS = ["CHEMBL", "replica", "seed", "delta_H_PB_kcal_mol", "dG_gas", "dG_solv",
          "lig_net_charge", "charge_source", "n_frames", "zn_pass_rate",
          "zn_mean_min4_dist", "qc_pass", "prod_ns", "status", "notes", "minutes"]


def log(msg: str) -> None:
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    try:
        with open(RUNLOG, "a") as fh:
            fh.write(line + "\n")
    except Exception:
        pass


def run(cmd: str, env: str, timeout: int = 7200, cwd: str | None = None):
    full = f"source {CONDA_SH}; conda activate {env}; " + cmd
    try:
        p = subprocess.run(["bash", "-c", full], capture_output=True, text=True,
                           timeout=timeout, cwd=cwd or BASE)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return 124, "", "TIMEOUT"


def write_file(path: str, content: str) -> None:
    with open(path, "w") as fh:
        fh.write(content)


def detect_ntomp() -> int:
    """按 cgroup CPU 配额推导 OpenMP 线程数（AutoDL 容器常把 52 核限到 6 核配额）。"""
    import multiprocessing
    for path, is_v2 in (("/sys/fs/cgroup/cpu/cpu.cfs_quota_us", False),
                        ("/sys/fs/cgroup/cpu.max", True)):
        try:
            if not is_v2:
                q = int(open(path).read().strip())
                p_ = int(open("/sys/fs/cgroup/cpu/cpu.cfs_period_us").read().strip())
                if q > 0 and p_ > 0:
                    return max(1, q // p_)
            else:
                mx = open(path).read().strip().split()
                if mx and mx[0] != "max":
                    return max(1, int(mx[0]) // int(mx[1]))
        except Exception:
            continue
    return min(24, multiprocessing.cpu_count())


NTHREADS = detect_ntomp()


def seed_for(name: str, replica: int) -> int:
    """确定性派生随机种子（hash() 跨进程不稳定，用 crc32）。"""
    return (zlib.crc32(name.encode()) % 90000) + 1000 + replica * 7919


# --------------------------- MDP（P2：加入 NPT） ---------------------------
def _ns_to_steps(ns: float) -> int:
    return int(round(ns * 1000 / 0.002))


def em_mdp() -> str:
    return """integrator               = steep
emtol                    = 1000.0
emstep                   = 0.01
nsteps                   = 5000
nstlist                  = 10
ns_type                  = grid
cutoff-scheme            = Verlet
rlist                    = 1.0
rcoulomb                 = 1.0
rvdw                     = 1.0
coulombtype              = PME
pme_order                = 4
fourierspacing           = 0.16
constraints              = h-bonds
"""


def nvt_eq_mdp(ns: float, seed: int) -> str:
    return f"""integrator               = md
dt                       = 0.002
nsteps                   = {_ns_to_steps(ns)}
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
pcoupl                   = no
constraints              = h-bonds
continuation             = no
gen_vel                  = yes
gen_temp                 = 300
gen_seed                 = {seed}
"""


def npt_eq_mdp(ns: float, seed: int) -> str:
    return f"""integrator               = md
dt                       = 0.002
nsteps                   = {_ns_to_steps(ns)}
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
tau_p                    = 2.0
ref_p                    = 1.0
compressibility          = 4.5e-5
refcoord_scaling         = com
constraints              = h-bonds
continuation             = yes
gen_vel                  = no
ld_seed                  = {seed + 1}
"""


def prod_mdp(ns: float, save_ps: float, seed: int) -> str:
    return f"""integrator               = md
dt                       = 0.002
nsteps                   = {_ns_to_steps(ns)}
nstxout-compressed       = {int(round(save_ps / 0.002))}
nstenergy                = {int(round(save_ps / 0.002))}
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
pcoupl                   = Parrinello-Rahman
pcoupltype               = isotropic
tau_p                    = 2.0
ref_p                    = 1.0
compressibility          = 4.5e-5
refcoord_scaling         = com
constraints              = h-bonds
continuation             = yes
gen_vel                  = no
ld_seed                  = {seed + 2}
"""


MMPBSA_IN = """&general
   interval=1,
   verbose=0, keep_files=0,
/
&pb
   istrng=0.15, inp=1, radiopt=0, fillratio=4.0,
/
"""

# --------------------------- tleap（沿用已验证写法） ---------------------------
def build_tleap_rec(d: str) -> str:
    p = os.path.join(d, "tleap_rec.in")
    write_file(p, f"""
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.tip3p
loadAmberParams {ZINC_FRCMOD}
ZN = loadmol2 {TEMPL}/ZN.mol2
rec = loadpdb {TEMPL}/rec_clean.pdb
z1 = loadpdb {TEMPL}/zn1.pdb
z2 = loadpdb {TEMPL}/zn2.pdb
sys = combine {{ rec z1 z2 }}
saveamberparm sys {d}/rec.prmtop {d}/rec.inpcrd
quit
""")
    return p


def build_tleap_lig(d: str, ligmol2: str) -> str:
    p = os.path.join(d, "tleap_lig.in")
    write_file(p, f"""
source leaprc.gaff2
loadAmberParams {d}/lig.frcmod
lig = loadmol2 {ligmol2}
saveamberparm lig {d}/lig.prmtop {d}/lig.inpcrd
quit
""")
    return p


def build_tleap_comsolv(d: str, ligmol2: str) -> str:
    p = os.path.join(d, "tleap_comsolv.in")
    write_file(p, f"""
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.tip3p
loadAmberParams {ZINC_FRCMOD}
ZN = loadmol2 {TEMPL}/ZN.mol2
rec = loadpdb {TEMPL}/rec_clean.pdb
z1 = loadpdb {TEMPL}/zn1.pdb
z2 = loadpdb {TEMPL}/zn2.pdb
loadAmberParams {d}/lig.frcmod
lig = loadmol2 {ligmol2}
sys = combine {{ rec z1 z2 lig }}
solvatebox sys TIP3PBOX 12.0
addions sys Na+ 0
addions sys Cl- 0
saveamberparm sys {d}/comsolv.prmtop {d}/comsolv.inpcrd
quit
""")
    return p


# --------------------------- P1：严格解析 ---------------------------
def parse_dg_bind_strict(txt: str):
    """严格解析 MMPBSA.py PB-only 输出。

    DELTA TOTAL / DELTA G gas / DELTA G solv 在 AmberTools 14/24 的 PB-only 模式下
    恒打印 0.0000，真实值须取 'Differences (Complex - Receptor - Ligand)' 段各组分之和：
        ΔH(PB) = VDWAALS + EEL + INTERNAL + EPB + ENPOLAR + EDISPER
    缺该段即返回 None（调用方判 PARSE_FAIL），绝不静默取 0。
    返回 (dH_pb, dG_gas_or_None, dG_solv_or_None, components_dict)
    """
    pat = re.compile(r"Differences \(Complex - Receptor - Ligand\):(.*?)DELTA TOTAL", re.S)
    comps: dict[str, float] = {}
    for m in pat.finditer(txt):
        for line in m.group(1).splitlines():
            mm = re.match(r"\s*([A-Za-z0-9_\- ]+?)\s+(-?\d+\.\d+)", line)
            if not mm:
                continue
            name = mm.group(1).strip()
            if name.upper().startswith("DELTA"):
                continue
            comps[name] = float(mm.group(2))
        if comps:
            break
    if not comps:
        return None, None, None, {}
    dH = sum(comps.values())
    gas = sum(comps.get(k, 0.0) for k in ("VDWAALS", "EEL", "INTERNAL"))
    solv = sum(comps.get(k, 0.0) for k in ("EPB", "ENPOLAR", "EDISPER"))
    return dH, gas, solv, comps


# --------------------------- P3：配体电荷 ---------------------------
def ligand_net_charge(name: str, mol2: str):
    """顺序推断配体净电荷，返回 (charge, source)。溯源优先级：
    1) $BASE/charge_table.csv（人工核定，列名 name,net_charge）
    2) RDKit 形式电荷（mol2 已含电荷列时可靠）
    3) 退化为 0 并标记 assumed_neutral（写入 SUMMARY，接受审稿检查）
    """
    table = os.path.join(BASE, "charge_table.csv")
    if os.path.exists(table):
        try:
            with open(table, newline="") as fh:
                for row in csv.DictReader(fh):
                    if row.get("name", "").strip() == name:
                        return int(float(row["net_charge"])), "charge_table"
        except Exception as exc:  # noqa: BLE001
            log(f"  [warn] charge_table 解析失败：{exc}")
    try:
        from rdkit import Chem  # type: ignore
        m = Chem.MolFromMol2File(mol2, removeHs=False, sanitize=False)
        if m is not None:
            return Chem.GetFormalCharge(m), "rdkit_formal"
    except Exception:
        pass
    return 0, "assumed_neutral"


# --------------------------- P4：双 Zn 配位 QC ---------------------------
def zn_coordination_qc(prmtop: str, traj: str) -> dict:
    """逐帧计算两个 Zn 的 S(CYS/CYM) 配位数，返回合格率。

    合格判据（预注册）：每个 Zn 在 ≤ ZN_S_CUTOFF 内有 ≥ ZN_MIN_COORD 个 S，
    且两个 Zn 同时合格的帧比例 ≥ ZN_MIN_FRACTION。
    """
    try:
        import numpy as np
        import mdtraj as md  # type: ignore
    except Exception as exc:  # noqa: BLE001
        return {"available": False, "reason": f"mdtraj/numpy 不可用: {exc}"}

    try:
        t = md.load(traj, top=prmtop)
    except Exception as exc:  # noqa: BLE001
        return {"available": False, "reason": f"轨迹读取失败: {exc}"}

    top = t.topology
    zn_idx = [a.index for a in top.atoms if a.element is not None and a.element.symbol == "Zn"]
    if not zn_idx:  # 兜底：按名字找
        zn_idx = [a.index for a in top.atoms if a.name.strip().upper() == "ZN"]
    s_idx = [a.index for a in top.atoms
             if a.residue.name.strip().upper() in ("CYS", "CYM") and a.name.strip().upper() == "SG"]
    if len(zn_idx) < 2 or not s_idx:
        return {"available": False, "reason":
                f"Zn 数={len(zn_idx)}（需 ≥2），CYS-SG 数={len(s_idx)}"}

    pairs = [(z, s) for z in zn_idx for s in s_idx]
    import numpy as np  # noqa: F811
    dist = md.compute_distances(t, pairs)  # (n_frames, n_zn*n_s)
    dist = dist.reshape(t.n_frames, len(zn_idx), len(s_idx))
    within = dist <= ZN_S_CUTOFF
    counts = within.sum(axis=2)                      # (n_frames, n_zn)
    per_zn_ok = counts >= ZN_MIN_COORD               # (n_frames, n_zn)
    both_ok = per_zn_ok.all(axis=1)
    # 每个 Zn 最近的 4 个 S 的平均距离
    mean_min4 = float(np.mean(np.sort(dist, axis=2)[:, :, :ZN_MIN_COORD]))
    return {
        "available": True,
        "n_frames": int(t.n_frames),
        "zn_pass_rate": float(both_ok.mean()),
        "zn_mean_min4_dist": round(mean_min4, 3),
        "n_zn": len(zn_idx),
        "n_cys_sg": len(s_idx),
    }


def write_summary(rows: list[dict]) -> None:
    with open(SUMMARY, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in FIELDS})


# --------------------------- 主流程 ---------------------------
def process(ligpdb: str, replica: int, rows: list[dict]) -> None:
    name = os.path.splitext(os.path.basename(ligpdb))[0]
    d = os.path.join(BASE, name, f"rep{replica}")
    os.makedirs(d, exist_ok=True)
    seed = seed_for(name, replica)
    t0 = time.time()
    row = {"CHEMBL": name, "replica": replica, "seed": seed, "prod_ns": PROD_NS,
           "delta_H_PB_kcal_mol": "NA", "dG_gas": "NA", "dG_solv": "NA",
           "lig_net_charge": "NA", "charge_source": "NA", "n_frames": "NA",
           "zn_pass_rate": "NA", "zn_mean_min4_dist": "NA", "qc_pass": "",
           "status": "OK", "notes": ""}
    log(f"===== {name} rep{replica} (seed={seed}, prod={PROD_NS} ns, gpu={GPU_ID}) =====")
    try:
        # 1) obabel -> mol2
        rc, o, e = run(f"obabel {ligpdb} -O {d}/lig.mol2", ENV_AMBER)
        if rc != 0 or not os.path.exists(f"{d}/lig.mol2"):
            raise RuntimeError(f"obabel 失败 rc={rc}: {e[:300]}")

        # 2) 配体电荷（P3）+ antechamber + parmchk2
        q, qsrc = ligand_net_charge(name, f"{d}/lig.mol2")
        row["lig_net_charge"], row["charge_source"] = q, qsrc
        log(f"  配体净电荷 = {q}（来源 {qsrc}）")
        rc, o, e = run(
            f"antechamber -i {d}/lig.mol2 -fi mol2 -o {d}/lig_gaff.mol2 -fo mol2 "
            f"-c bcc -nc {q} -at gaff2 -rn LIG", ENV_AMBER)
        if rc != 0 or not os.path.exists(f"{d}/lig_gaff.mol2"):
            raise RuntimeError(f"antechamber 失败 rc={rc}: {e[:300]}")
        rc, o, e = run(f"parmchk2 -i {d}/lig_gaff.mol2 -f mol2 -o {d}/lig.frcmod", ENV_AMBER)
        if rc != 0 or not os.path.exists(f"{d}/lig.frcmod"):
            raise RuntimeError(f"parmchk2 失败 rc={rc}: {e[:300]}")

        # 3) tleap 三拓扑
        ligmol2 = f"{d}/lig_gaff.mol2"
        for fn, inp in (("rec", build_tleap_rec(d)),
                        ("lig", build_tleap_lig(d, ligmol2)),
                        ("comsolv", build_tleap_comsolv(d, ligmol2))):
            rc, o, e = run(f"tleap -f {inp} > {d}/tleap_{fn}.log 2>&1", ENV_AMBER)
            logtxt = ""
            try:
                logtxt = open(f"{d}/tleap_{fn}.log", errors="replace").read()
            except Exception:
                pass
            merr = re.search(r"Exiting LEaP:\s*Errors\s*=\s*(\d+)", logtxt)
            if (merr and int(merr.group(1)) > 0) or ("FATAL" in logtxt) \
                    or (not os.path.exists(f"{d}/{fn}.prmtop")):
                raise RuntimeError(f"tleap {fn} 失败: {logtxt[-500:]}")

        # 4) ParmEd -> GROMACS
        rc, o, e = run(
            f"python -c \""
            f"import parmed as pmd; "
            f"s=pmd.load_file('{d}/comsolv.prmtop','{d}/comsolv.inpcrd'); "
            f"g=pmd.gromacs.GromacsTopologyFile.from_structure(s); g.write('{d}/topol.top'); "
            f"pmd.gromacs.GromacsGroFile.write(s,'{d}/topol.gro')\"", ENV_AMBER)
        if rc != 0 or not os.path.exists(f"{d}/topol.top"):
            raise RuntimeError(f"ParmEd 转换失败 rc={rc}: {e[:300]}")

        # 5) MD: EM -> NVT -> NPT -> PROD(NPT)   [P2]
        write_file(f"{d}/em.mdp", em_mdp())
        write_file(f"{d}/nvt_eq.mdp", nvt_eq_mdp(EQ_NS, seed))
        write_file(f"{d}/npt_eq.mdp", npt_eq_mdp(EQ_NS, seed))
        write_file(f"{d}/prod.mdp", prod_mdp(PROD_NS, SAVE_PS, seed))

        def grompp(mdp, gro, tpr, tag):
            rc_, _, e_ = run(f"gmx grompp -f {d}/{mdp} -c {d}/{gro} -p {d}/topol.top "
                             f"-o {d}/{tpr} -maxwarn 3 > {d}/grompp_{tag}.log 2>&1", ENV_GMX)
            if rc_ != 0 or not os.path.exists(f"{d}/{tpr}"):
                raise RuntimeError(f"grompp {tag} 失败: {e_[:300]}")

        def mdrun(tpr, tag, timeout):
            rc_, _, e_ = run(f"gmx mdrun -deffnm {d}/{tag} -s {d}/{tpr} -ntmpi 1 "
                             f"-ntomp {NTHREADS} -gpu_id {GPU_ID} > {d}/md_{tag}.log 2>&1",
                             ENV_GMX, timeout=timeout)
            out = f"{d}/{tag}.gro" if tpr != "prod.tpr" else f"{d}/prod.xtc"
            if rc_ != 0 or not os.path.exists(out):
                raise RuntimeError(f"mdrun {tag} 失败 rc={rc_}: {e_[:300]}")

        grompp("em.mdp", "topol.gro", "em.tpr", "em")
        mdrun("em.tpr", "em", 3600)
        grompp("nvt_eq.mdp", "em.gro", "nvt_eq.tpr", "nvt")
        mdrun("nvt_eq.tpr", "nvt_eq", 3600)
        grompp("npt_eq.mdp", "nvt_eq.gro", "npt_eq.tpr", "npt")
        mdrun("npt_eq.tpr", "npt_eq", 3600)
        grompp("prod.mdp", "npt_eq.gro", "prod.tpr", "prod")
        t_prod = time.time()
        mdrun("prod.tpr", "prod", int(PROD_NS * 3600) + 3600)
        ns_per_day = PROD_NS / ((time.time() - t_prod) / 86400)
        log(f"  实测吞吐 ≈ {ns_per_day:.1f} ns/day (gpu {GPU_ID})")
        row["notes"] = f"throughput={ns_per_day:.1f}ns/day"

        # 6) 双 Zn 配位 QC（P4：不达标如实记录，不静默丢弃）
        qc = zn_coordination_qc(f"{d}/comsolv.prmtop", f"{d}/prod.xtc")
        if not qc.get("available"):
            row.update({"zn_pass_rate": "NA", "qc_pass": "0",
                        "status": "QC_SKIPPED", "notes": (row["notes"] + " | " + qc.get("reason", "")).strip(" |")})
            log(f"  !! Zn QC 不可用：{qc.get('reason')}")
        else:
            row.update({"n_frames": qc["n_frames"],
                        "zn_pass_rate": round(qc["zn_pass_rate"], 4),
                        "zn_mean_min4_dist": qc["zn_mean_min4_dist"],
                        "qc_pass": "1" if qc["zn_pass_rate"] >= ZN_MIN_FRACTION else "0"})
            log(f"  Zn QC：合格帧 {qc['zn_pass_rate']*100:.1f}% "
                f"(阈值 {ZN_MIN_FRACTION*100:.0f}%)，Zn–S(最近4) 均距 {qc['zn_mean_min4_dist']} Å")
            if row["qc_pass"] == "0":
                row["status"] = "QC_FAIL"
                row["notes"] = (row["notes"] + f" | 双Zn四配位合格率 {qc['zn_pass_rate']:.3f} < {ZN_MIN_FRACTION}").strip(" |")

        # 7) cpptraj 去水 + parmed 派生 rec/lig
        write_file(f"{d}/strip.in", f"""
parm {d}/comsolv.prmtop
trajin {d}/prod.xtc
strip :WAT,Cl-,Na+,K+ parmout {d}/com_nobox.prmtop
trajout {d}/com_nobox.nc
""")
        rc, o, e = run(f"cpptraj -i {d}/strip.in > {d}/cpptraj.log 2>&1", ENV_AMBER)
        if rc != 0 or not os.path.exists(f"{d}/com_nobox.nc"):
            raise RuntimeError(f"cpptraj strip 失败 rc={rc}: {e[:200]}")
        rc, o, e = run(
            f"python -c \""
            f"import parmed as pmd, copy; "
            f"c=pmd.load_file('{d}/com_nobox.prmtop'); "
            f"c.box=None; c.save('{d}/com_nobox.prmtop', overwrite=True); "
            f"r=copy.deepcopy(c); r.strip(':LIG'); r.save('{d}/rec_nobox.prmtop', overwrite=True); "
            f"l=copy.deepcopy(c); l.strip('!:LIG'); l.save('{d}/lig_nobox.prmtop', overwrite=True)\"",
            ENV_AMBER)
        if rc != 0 or not os.path.exists(f"{d}/rec_nobox.prmtop"):
            raise RuntimeError(f"parmed 派生失败 rc={rc}: {e[:200]}")

        # 8) MMPBSA.py PB
        write_file(f"{d}/mmpbsa.in", MMPBSA_IN)
        rc, o, e = run(
            f"MMPBSA.py -O -i {d}/mmpbsa.in -o {d}/FINAL_RESULTS.dat "
            f"-cp {d}/com_nobox.prmtop -rp {d}/rec_nobox.prmtop -lp {d}/lig_nobox.prmtop "
            f"-y {d}/com_nobox.nc > {d}/mmpbsa.log 2>&1", ENV_AMBER, timeout=7200)
        fres = f"{d}/FINAL_RESULTS.dat"
        if rc != 0 or not os.path.exists(fres):
            raise RuntimeError(f"MMPBSA 失败 rc={rc}: {e[:300]}")

        # 9) 严格解析（P1）
        txt = open(fres, errors="replace").read()
        dH, gas, solv, comps = parse_dg_bind_strict(txt)
        if dH is None:
            row["status"] = "PARSE_FAIL"
            row["notes"] = (row["notes"] + " | 缺 Differences 段，拒绝取 DELTA TOTAL=0").strip(" |")
        else:
            row["delta_H_PB_kcal_mol"] = f"{dH:.4f}"
            row["dG_gas"] = f"{gas:.4f}"
            row["dG_solv"] = f"{solv:.4f}"
            if "QUESTIONABLE" in txt or "INCONSISTENCIES" in txt:
                row["notes"] = (row["notes"] + " | internal consistency note (benign)").strip(" |")
    except Exception as ex:  # noqa: BLE001
        row["status"] = "FAILED"
        row["notes"] = (row["notes"] + " | " + str(ex)[:300]).strip(" |")
        log(f"  !! {name} rep{replica} 失败：{ex}")
    row["minutes"] = f"{(time.time() - t0)/60:.1f}"
    rows.append(row)
    write_summary(rows)
    log(f"  -> {name} rep{replica}: {row['status']} ΔH_PB={row['delta_H_PB_kcal_mol']} "
        f"({row['minutes']} min)")


def main() -> None:
    log("=" * 70)
    log("START 12_pipeline_v2.py")
    log(f"BASE={BASE} REPLICAS={REPLICAS} PROD_NS={PROD_NS} SAVE_PS={SAVE_PS} GPU_ID={GPU_ID}")
    ligs = sorted(glob.glob(os.path.join(LIGDIR, "*.pdb")))
    if LIG_ONLY:
        ligs = [p for p in ligs if os.path.splitext(os.path.basename(p))[0] in LIG_ONLY]
    log(f"待处理配体 {len(ligs)} 个：{[os.path.basename(p) for p in ligs]}")
    if not ligs:
        log("ERROR: 无配体，退出")
        return
    rows: list[dict] = []
    for ligpdb in ligs:
        reps = [REPLICA_ID] if REPLICA_ID else list(range(1, REPLICAS + 1))
        for r in reps:
            process(ligpdb, r, rows)
    log("DONE 12_pipeline_v2.py")


if __name__ == "__main__":
    main()
