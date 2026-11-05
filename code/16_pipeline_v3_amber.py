#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""16_pipeline_v3_amber.py — B 线重跑管线 v3：AmberTools 单栈（远端自运行）。

为什么不再用 GROMACS
--------------------
MatPool 实例实测 GitHub / ftp.gromacs.org 不可达，无法拿到 GROMACS 源码编译 GPU 版；
而 conda-forge 的 AmberTools 自带 **pmemd.cuda**（GPU MD 引擎）。
因此改为单栈：

    tleap(拓扑) → pmemd.cuda(EM/NVT/NPT/生产) → cpptraj(去水) → MMPBSA.py(PB)

比 v2 少一步 ParmEd→GROMACS 格式转换，也少一处出错点。

沿用 v2 的全部修复
------------------
P1 严格 ΔH 解析（拒绝 DELTA TOTAL=0.0000）
P2 生产在 NPT（v3：ntp=1 / barostat=2），不是 NVT
P3 配体电荷三级溯源（charge_table.csv → RDKit → 退化 0 并标记）
P4 双 Zn 四配位 QC（mdtraj 逐帧，阈值 Zn–S ≤2.8 Å、≥4 个 S、合格帧 ≥95%）
P5 每 replica 显式随机种子（crc32 确定性派生）

环境变量（同 v2）
-----------------
BASE, TEMPL_DIR, LIGDIR, CONDA_SH, ENV_AMBER, GPU_ID, REPLICAS, PROD_NS, SAVE_PS, EQ_NS, LIG_ONLY
ZN_S_CUTOFF, ZN_MIN_COORD, ZN_MIN_FRACTION
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
LIGDIR = os.environ.get("LIGDIR", os.path.join(BASE, "ligands_vina"))
CONDA_SH = os.environ.get("CONDA_SH", "/root/miniconda3/etc/profile.d/conda.sh")
ZINC_FRCMOD = os.environ.get("ZINC_FRCMOD", os.path.join(TEMPL, "zinc_nonbond_fallback.frcmod"))
ENV_AMBER = os.environ.get("ENV_AMBER", "gamd")
GPU_ID = os.environ.get("GPU_ID", "0")
REPLICAS = int(os.environ.get("REPLICAS", "1"))
PROD_NS = float(os.environ.get("PROD_NS", "30"))
SAVE_STEPS = int(round(float(os.environ.get("SAVE_PS", "20")) / 0.002))
HEAT_NS = float(os.environ.get("HEAT_NS", "0.2"))
EQ_NS = float(os.environ.get("EQ_NS", "1"))

# 单位 nm（mdtraj 返回 nm）。历史缺陷：曾写 2.8(Å) 用于 nm 距离 → 阈值 28 Å，门失效。
if os.environ.get("ZN_S_CUTOFF_ANGSTROM"):
    ZN_S_CUTOFF = float(os.environ["ZN_S_CUTOFF_ANGSTROM"]) / 10.0
else:
    ZN_S_CUTOFF = float(os.environ.get("ZN_S_CUTOFF", "0.28"))  # nm
ZN_MIN_COORD = int(os.environ.get("ZN_MIN_COORD", "4"))
ZN_MIN_FRACTION = float(os.environ.get("ZN_MIN_FRACTION", "0.95"))

RUNLOG = os.path.join(BASE, "run_v3.log")
SUMMARY = os.path.join(BASE, "SUMMARY_v3.csv")
FIELDS = ["CHEMBL", "replica", "seed", "delta_H_PB_kcal_mol", "dG_gas", "dG_solv",
          "lig_net_charge", "charge_source", "n_frames", "zn_pass_rate",
          "zn_mean_min4_dist", "qc_pass", "prod_ns", "ns_per_day", "status", "notes", "minutes"]


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
    full = f"source {CONDA_SH}; conda activate {env}; export CUDA_VISIBLE_DEVICES={GPU_ID}; " + cmd
    try:
        p = subprocess.run(["bash", "-c", full], capture_output=True, text=True,
                           timeout=timeout, cwd=cwd or BASE)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return 124, "", "TIMEOUT"


def write_file(path: str, content: str) -> None:
    with open(path, "w") as fh:
        fh.write(content)


def seed_for(name: str, replica: int) -> int:
    return (zlib.crc32(name.encode()) % 90000) + 1000 + replica * 7919


def ns_to_steps(ns: float) -> int:
    return int(round(ns * 1000 / 0.002))


# ---------------------------- pmemd 输入（P2：NPT 生产） ----------------------------
def min_in(maxcyc: int, restraint: bool, outroot: str) -> str:
    ntr = 1 if restraint else 0
    return f"""Minimization ({'restrained' if restraint else 'free'})
 &cntrl
  imin=1, maxcyc={maxcyc}, ncyc={max(1, maxcyc // 2)}, ntmin=2,
  ntb=1, cut=10.0, ntpr=100,
  ntr={ntr}, restraint_wt=5.0, restraintmask='!:WAT,Cl-,Na+,K+ & !@H=',
 &end
"""


def md_in(nsteps: int, ntt_is_md: bool, npt: bool, istep_interval: int, seed: int, restart: bool) -> str:
    ntb = 2 if npt else 1
    return f"""MD {'NPT' if npt else 'NVT'}
 &cntrl
  imin=0, irest={1 if restart else 0}, ntx={1 if restart else 1},
  nstlim={nsteps}, dt=0.002,
  ntb={ntb}, ntp={1 if npt else 0}, barostat=2, mcbarint=100,
  taup=2.0, pres0=1.0, comp=44.6e-6,
  ntt=3, temp0=300.0, tempi=300.0, gamma_ln=2.0, ig={seed},
  ntc=2, ntf=2, cut=10.0,
  ntpr={istep_interval}, ntwx={istep_interval}, ntwr={max(istep_interval * 5, 5000)},
  ioutfm=1, ntxo=1, ntwprt=0,
 &end
"""


# ---------------------------- tleap ----------------------------
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


# ---------------------------- P1 严格解析 ----------------------------
def parse_dg_bind_strict(txt: str):
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


# ---------------------------- P3 配体电荷 ----------------------------
def ligand_net_charge(name: str, mol2: str):
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


# ---------------------------- P4 双 Zn QC ----------------------------
def zn_coordination_qc(prmtop: str, traj: str) -> dict:
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
    zn_idx = [a.index for a in top.atoms if a.element is not None and a.element.symbol == "Zn"] \
        or [a.index for a in top.atoms if a.name.strip().upper() == "ZN"]
    s_idx = [a.index for a in top.atoms
             if a.residue.name.strip().upper() in ("CYS", "CYM") and a.name.strip().upper() == "SG"]
    if len(zn_idx) < 2 or not s_idx:
        return {"available": False, "reason": f"Zn 数={len(zn_idx)}（需 ≥2），CYS-SG={len(s_idx)}"}
    pairs = [(z, s) for z in zn_idx for s in s_idx]
    dist = md.compute_distances(t, pairs).reshape(t.n_frames, len(zn_idx), len(s_idx))
    counts = (dist <= ZN_S_CUTOFF).sum(axis=2)
    both_ok = (counts >= ZN_MIN_COORD).all(axis=1)
    mean_min4 = float(np.mean(np.sort(dist, axis=2)[:, :, :ZN_MIN_COORD]))
    return {"available": True, "n_frames": int(t.n_frames),
            "zn_pass_rate": float(both_ok.mean()),
            "zn_mean_min4_dist": round(mean_min4, 3),
            "n_zn": len(zn_idx), "n_cys_sg": len(s_idx)}


def write_summary(rows: list[dict]) -> None:
    with open(SUMMARY, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in FIELDS})


# ---------------------------- 主流程 ----------------------------
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
           "ns_per_day": "NA", "status": "OK", "notes": ""}
    log(f"===== {name} rep{replica} (seed={seed}, {PROD_NS} ns, GPU {GPU_ID}) =====")
    try:
        rc, o, e = run(f"obabel {ligpdb} -O {d}/lig.mol2", ENV_AMBER)
        if rc != 0 or not os.path.exists(f"{d}/lig.mol2"):
            raise RuntimeError(f"obabel 失败 rc={rc}: {e[:300]}")

        q, qsrc = ligand_net_charge(name, f"{d}/lig.mol2")
        row["lig_net_charge"], row["charge_source"] = q, qsrc
        log(f"  配体净电荷={q}（来源 {qsrc}）")

        rc, o, e = run(f"antechamber -i {d}/lig.mol2 -fi mol2 -o {d}/lig_gaff.mol2 -fo mol2 "
                       f"-c bcc -nc {q} -at gaff2 -rn LIG", ENV_AMBER)
        if rc != 0 or not os.path.exists(f"{d}/lig_gaff.mol2"):
            raise RuntimeError(f"antechamber 失败 rc={rc}: {e[:300]}")
        rc, o, e = run(f"parmchk2 -i {d}/lig_gaff.mol2 -f mol2 -o {d}/lig.frcmod", ENV_AMBER)
        if rc != 0 or not os.path.exists(f"{d}/lig.frcmod"):
            raise RuntimeError(f"parmchk2 失败 rc={rc}: {e[:300]}")

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

        # ---- pmemd.cuda 四步：min1(约束) → min2(自由) → heat(NVT) → eq(NPT) → prod(NPT) ----
        write_file(f"{d}/min1.in", min_in(2500, True, "min1"))
        write_file(f"{d}/min2.in", min_in(2500, False, "min2"))
        write_file(f"{d}/heat.in", md_in(ns_to_steps(HEAT_NS), True, False, 0, seed, False))
        write_file(f"{d}/eq.in", md_in(ns_to_steps(EQ_NS), True, True, 0, seed + 1, True))
        write_file(f"{d}/prod.in", md_in(ns_to_steps(PROD_NS), True, True, SAVE_STEPS, seed + 2, True))

        def pmemd(inp: str, crd_in: str, out: str, ref: str | None, tag: str, timeout: int):
            refc = f"-ref {ref}" if ref else ""
            cmd = (f"pmemd.cuda -O -i {d}/{inp} -o {d}/{tag}.out -p {d}/comsolv.prmtop "
                   f"-c {crd_in} {refc} -r {d}/{tag}.rst7 -x {d}/{tag}.nc "
                   f"-inf {d}/{tag}.mdinfo > {d}/pmemd_{tag}.log 2>&1")
            rc_, _, e_ = run(cmd, ENV_AMBER, timeout=timeout)
            if rc_ != 0 or not os.path.exists(f"{d}/{tag}.rst7"):
                tail = ""
                try:
                    tail = open(f"{d}/{tag}.out", errors="replace").read()[-800:]
                except Exception:
                    pass
                raise RuntimeError(f"pmemd {tag} 失败 rc={rc_}: {tail or e_[:300]}")

        pmemd("min1.in", f"{d}/comsolv.inpcrd", "min1", f"{d}/comsolv.inpcrd", "min1", 3600)
        pmemd("min2.in", f"{d}/min1.rst7", "min2", None, "min2", 3600)
        pmemd("heat.in", f"{d}/min2.rst7", "heat", None, "heat", 3600)
        pmemd("eq.in", f"{d}/heat.rst7", "eq", None, "eq", 7200)
        t_prod = time.time()
        pmemd("prod.in", f"{d}/eq.rst7", "prod", None, "prod", int(PROD_NS * 3600 * 3) + 3600)
        npd = PROD_NS / ((time.time() - t_prod) / 86400)
        row["ns_per_day"] = round(npd, 1)
        log(f"  实测吞吐 ≈ {npd:.1f} ns/day")
        row["notes"] = f"throughput={npd:.1f}ns/day"

        qc = zn_coordination_qc(f"{d}/comsolv.prmtop", f"{d}/prod.nc")
        if not qc.get("available"):
            row.update({"qc_pass": "0", "status": "QC_SKIPPED",
                        "notes": (row["notes"] + " | " + qc.get("reason", "")).strip(" |")})
            log(f"  !! Zn QC 不可用：{qc.get('reason')}")
        else:
            row.update({"n_frames": qc["n_frames"], "zn_pass_rate": round(qc["zn_pass_rate"], 4),
                        "zn_mean_min4_dist": qc["zn_mean_min4_dist"],
                        "qc_pass": "1" if qc["zn_pass_rate"] >= ZN_MIN_FRACTION else "0"})
            log(f"  Zn QC：合格帧 {qc['zn_pass_rate']*100:.1f}%（阈值 {ZN_MIN_FRACTION*100:.0f}%），"
                f"Zn–S(最近4) 均距 {qc['zn_mean_min4_dist']} Å")
            if row["qc_pass"] == "0":
                row["status"] = "QC_FAIL"
                row["notes"] = (row["notes"] +
                                f" | 双Zn四配位合格率 {qc['zn_pass_rate']:.3f} < {ZN_MIN_FRACTION}").strip(" |")

        write_file(f"{d}/strip.in", f"""
parm {d}/comsolv.prmtop
trajin {d}/prod.nc
strip :WAT,Cl-,Na+,K+ parmout {d}/com_nobox.prmtop
trajout {d}/com_nobox.nc netcdf
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

        write_file(f"{d}/mmpbsa.in", """&general
   interval=1, verbose=0, keep_files=0,
/
&pb
   istrng=0.15, inp=1, radiopt=0, fillratio=4.0,
/
""")
        rc, o, e = run(f"MMPBSA.py -O -i {d}/mmpbsa.in -o {d}/FINAL_RESULTS.dat "
                       f"-cp {d}/com_nobox.prmtop -rp {d}/rec_nobox.prmtop "
                       f"-lp {d}/lig_nobox.prmtop -y {d}/com_nobox.nc > {d}/mmpbsa.log 2>&1",
                       ENV_AMBER, timeout=7200)
        fres = f"{d}/FINAL_RESULTS.dat"
        if rc != 0 or not os.path.exists(fres):
            raise RuntimeError(f"MMPBSA 失败 rc={rc}: {e[:300]}")

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
    log(f"  -> {name} rep{replica}: {row['status']} ΔH_PB={row['delta_H_PB_kcal_mol']} ({row['minutes']} min)")


def main() -> None:
    log("=" * 70)
    log("START 16_pipeline_v3_amber.py")
    log(f"BASE={BASE} LIGDIR={LIGDIR} REPLICAS={REPLICAS} PROD_NS={PROD_NS} GPU_ID={GPU_ID}")
    ligs = sorted(glob.glob(os.path.join(LIGDIR, "*.pdb")))
    lig_only = [s.strip() for s in os.environ.get("LIG_ONLY", "").split(",") if s.strip()]
    if lig_only:
        ligs = [p for p in ligs if os.path.splitext(os.path.basename(p))[0] in lig_only]
    log(f"待处理配体 {len(ligs)} 个")
    if not ligs:
        log("ERROR: 无配体，退出")
        return
    rows: list[dict] = []
    for ligpdb in ligs:
        for r in range(1, REPLICAS + 1):
            process(ligpdb, r, rows)
    log("DONE 16_pipeline_v3_amber.py")


if __name__ == "__main__":
    main()
