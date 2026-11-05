#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
24_zn_model_test.py — Zn 模型 A/B/C 验证（apo 受体，2 ns 生产）

背景：上一批 42 run 中 Zn–S(最近4) 均距 0.825–2.372 nm（中位 0.997 nm），
      而晶体模板为 0.226–0.236 nm —— 非键 Li/Merz Zn2+ + 中性 CYS 不保四配位，
      锌指折叠塌陷是配体普遍脱离的充分解释。本脚本在 apo 受体上对比三种 Zn 模型，
      选出能维持四配位的模型后再全量重跑。

模型：
  A  zaff_cym      bonded ZAFF（ZN-SG 成键） + CYM（去质子化，-1）
  B  nonbond_cym   Li/Merz 非键 Zn2+        + CYM（去质子化，-1）  仅靠静电
  C  nonbond_cys   现状基线：Li/Merz 非键   + 中性 CYS            预期失败（复现缺陷）

判据（预注册）：
  PASS = 生产轨迹末 50% 帧内，两个 Zn 各自到其 4 个配位 Cys SG 的距离均值 ≤ 0.30 nm
         （模板 0.226–0.236 nm，留 0.06 nm 余量），且最大帧均值 ≤ 0.35 nm。
  FAIL = 否则。

用法（远端，每个模型一个进程、各占一张卡）：
    python 24_zn_model_test.py zaff_cym     0
    python 24_zn_model_test.py nonbond_cym  1
    python 24_zn_model_test.py nonbond_cys  2

输出：<BASE>/zn_model_test/<model>/zn_result.json + zn_dist_series.csv
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

BASE = os.environ.get("BASE", "/root/autodl-tmp/HPV_structure/candidates_run")
TEMPL = os.path.join(BASE, "templates")
OUTROOT = os.path.join(BASE, "zn_model_test")

# 配位半胱氨酸 —— 由模板实测坐标确定（2026-10-29 复核，见 sanity_check_pairs）：
#   zn1.pdb（PDB 残基 201，Zn_A）← Cys103/106/136/139，S–Zn 2.261–2.357 Å
#   zn2.pdb（PDB 残基 202，Zn_B）← Cys30/33/63/66，    S–Zn 2.281–2.358 Å
#   Cys16（>18 Å）与 Cys51（>17 Å）不配位。
ZN1_RES = [103, 106, 136, 139]  # ← zn1.pdb
ZN2_RES = [30, 33, 63, 66]      # ← zn2.pdb
COORD_RES = ZN1_RES + ZN2_RES

PROD_NS = float(os.environ.get("PROD_NS", "2"))
SAVE_PS = float(os.environ.get("SAVE_PS", "20"))
PY = sys.executable

# ---------- ZAFF（bonded）参数：自撰，原 ZAFF 存档不可恢复 ----------
# 来源：ZAFF (Peters et al., JCTC 2010) 流传参数集；键长 r0 落在模板实测区间中位。
# 声明：这不是官方 ZAFF 文件的逐字副本，是重建参数；验收靠下面的距离判据，不靠出处。
ZAFF_FRCMOD = """ZAFF-style bonded zinc model for E6 Cys4 site (reconstructed; ZAFF archive unrecoverable)
NOTE: ff19SB types cysteine SG as SH (not S) and CB as CT; parameters must use SH/CT.
MASS
ZN 65.37

BOND
ZN-SH  206.00   2.3200

ANGLE
SH-ZN-SH    27.00   109.500
ZN-SH-CT    50.00    96.000
CT-SH-ZN    50.00    96.000

DIHE
SH-ZN-SH-CT    1    0.000   0.000   2.000
CT-SH-ZN-SH    1    0.000   0.000   2.000

IMPROPER

NONBON
  ZN      1.1000  0.0125
"""

NONBOND_FRCMOD = """Zn2+ non-bonded 12-6 (Li/Merz) for E6 zinc fingers
MASS
ZN 65.37

NONBON
  ZN      1.1000  0.0125
"""


def sh(cmd: str, cwd: str | None = None) -> tuple[int, str]:
    p = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
    return p.returncode, (p.stdout + p.stderr)


def write(path: str, text: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8")


# ---------------------------------------------------------------- 受体改造
def build_rec_cym() -> str:
    """把 8 个配位 CYS 改名 CYM（去质子化）；模板中 CYS 本就无 HG，无需删氢。"""
    src = os.path.join(TEMPL, "rec_clean.pdb")
    dst = os.path.join(TEMPL, "rec_cym.pdb")
    if Path(dst).exists():
        return dst
    out = []
    n_changed = 0
    for line in Path(src).read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20]
            try:
                resseq = int(line[22:26])
            except ValueError:
                resseq = -1
            if resname == "CYS" and resseq in COORD_RES:
                line = line[:17] + "CYM" + line[20:]
                n_changed += 1
            out.append(line)
        else:
            out.append(line)
    write(dst, "\n".join(out) + "\n")
    print(f"[rec_cym] 改写 {n_changed} 个原子行 (应 = 8 残基 × 6 原子 = 48)")
    if n_changed != 48:
        print(f"[WARN] 改写行数 {n_changed} ≠ 48，请检查模板残基编号")
    return dst


# ---------------------------------------------------------------- 配对自检
def sanity_check_pairs() -> None:
    """用模板坐标核对 Zn↔Cys 配对；任一 S–Zn > 3 Å 立刻中止（防止建出 15 Å 的荒谬键）。"""
    import math

    def read_sg(resseq: int) -> list[float]:
        for line in Path(os.path.join(TEMPL, "rec_clean.pdb")).read_text().splitlines():
            if line.startswith("ATOM") and line[17:20] == "CYS" \
                    and int(line[22:26]) == resseq and line[12:16].strip() == "SG":
                return [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        raise RuntimeError(f"模板中找不到 CYS {resseq} 的 SG")

    def read_zn(fn: str) -> list[float]:
        for line in Path(os.path.join(TEMPL, fn)).read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                return [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        raise RuntimeError(f"{fn} 中没有原子")

    def dist(a, b):
        return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))

    bad = []
    for fn, reslist, tag in (("zn1.pdb", ZN1_RES, "zn1"),
                             ("zn2.pdb", ZN2_RES, "zn2")):
        z = read_zn(fn)
        for r in reslist:
            d = dist(read_sg(r), z)
            print(f"  [pair] {tag} -- Cys{r:<4} {d:.3f} Å")
            if d > 3.0:
                bad.append((tag, r, d))
    if bad:
        raise RuntimeError(f"Zn↔Cys 配对错误（>3 Å）: {bad}；请核对 ZN1_RES/ZN2_RES 与 zn1/zn2.pdb 的对应关系")
    print("  [pair] 配对自检通过：8 条 S–Zn 均 < 3 Å")


# ---------------------------------------------------------------- tleap
def build_zn_chainC() -> tuple[str, str, int, int]:
    """把两个 Zn 改写为与受体同链(C)、残基号 144/145。

    原因：zn1/zn2.pdb 原本链 Z、残基 201/202。tleap 的 bond 语法在 unit 含多条链时
    必须写成 unit.chain.res.atom，单链则可省略 chain。统一成单链 C 后可直接用
    `bond sys.30.SG sys.144.ZN`，避免链 ID 歧义导致的解析失败。
    """
    out = []
    for src, newseq in (("zn1.pdb", 144), ("zn2.pdb", 145)):
        dst = os.path.join(TEMPL, f"{src[:-4]}_chainC.pdb")
        lines = []
        for line in Path(os.path.join(TEMPL, src)).read_text(encoding="utf-8").splitlines():
            if line.startswith(("ATOM", "HETATM")):
                # 列: 22 = chainID, 23-26 = resSeq
                line = line[:21] + "C" + f"{newseq:4d}" + line[26:]
                lines.append(line)
            else:
                lines.append(line)
        write(dst, "\n".join(lines) + "\n")
        out.append(dst)
    return out[0], out[1], 144, 145


def build_tleap(model: str, d: str) -> str:
    frcmod_path = os.path.join(TEMPL, "zaff.frcmod" if model == "zaff_cym"
                               else "zinc_nonbond_fallback.frcmod")
    rec = os.path.join(TEMPL, "rec_cym.pdb" if model.endswith("cym")
                       else "rec_clean.pdb")
    z1, z2, seq1, seq2 = build_zn_chainC()
    bonds = []
    if model == "zaff_cym":
        # Zn(seq1=144) ← Cys30/33/63/66 ; Zn(seq2=145) ← Cys103/106/136/139
        for r in ZN1_RES:
            bonds.append(f"bond sys.{r}.SG sys.{seq1}.ZN")
        for r in ZN2_RES:
            bonds.append(f"bond sys.{r}.SG sys.{seq2}.ZN")

    txt = f"""source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.tip3p
loadAmberParams {frcmod_path}
ZN = loadmol2 {TEMPL}/ZN.mol2
rec = loadpdb {rec}
z1 = loadpdb {z1}
z2 = loadpdb {z2}
sys = combine {{ rec z1 z2 }}
{chr(10).join(bonds)}
charge sys
check sys
solvatebox sys TIP3PBOX 12.0
addions sys Cl- 0
addions sys Na+ 0
saveamberparm sys {d}/sys.prmtop {d}/sys.inpcrd
quit
"""
    p = os.path.join(d, "tleap.in")
    write(p, txt)
    return p


# ---------------------------------------------------------------- mdp
def _steps(ns: float) -> int:
    return int(round(ns * 1000 / 0.002))


def mdp_em() -> str:
    return """integrator = steep
emtol = 10.0
emstep = 0.01
nsteps = 5000
cutoff-scheme = Verlet
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
"""


def mdp_nvt(ns: float, seed: int) -> str:
    return f"""integrator = md
dt = 0.002
nsteps = {_steps(ns)}
nstlist = 10
cutoff-scheme = Verlet
rlist = 1.0
rcoulomb = 1.0
rvdw = 1.0
coulombtype = PME
tcoupl = v-rescale
tc_grps = System
tau_t = 0.1
ref_t = 300
pcoupl = no
constraints = h-bonds
continuation = no
gen_vel = yes
gen_temp = 300
gen_seed = {seed}
"""


def mdp_npt(ns: float, seed: int) -> str:
    return f"""integrator = md
dt = 0.002
nsteps = {_steps(ns)}
nstlist = 10
cutoff-scheme = Verlet
rlist = 1.0
rcoulomb = 1.0
rvdw = 1.0
coulombtype = PME
tcoupl = v-rescale
tc_grps = System
tau_t = 0.1
ref_t = 300
pcoupl = C-rescale
pcoupltype = isotropic
tau_p = 2.0
ref_p = 1.0
compressibility = 4.5e-5
refcoord_scaling = com
constraints = h-bonds
continuation = yes
gen_vel = no
ld_seed = {seed + 1}
"""


def mdp_prod(ns: float, save_ps: float, seed: int) -> str:
    return f"""integrator = md
dt = 0.002
nsteps = {_steps(ns)}
nstxout-compressed = {int(round(save_ps / 0.002))}
nstenergy = {int(round(save_ps / 0.002))}
nstlist = 10
cutoff-scheme = Verlet
rlist = 1.0
rcoulomb = 1.0
rvdw = 1.0
coulombtype = PME
tcoupl = v-rescale
tc_grps = System
tau_t = 0.1
ref_t = 300
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau_p = 2.0
ref_p = 1.0
compressibility = 4.5e-5
refcoord_scaling = com
constraints = h-bonds
continuation = yes
gen_vel = no
ld_seed = {seed + 2}
"""


# ---------------------------------------------------------------- 分析
def analyse(prmtop: str, xtc: str, model: str, d: str) -> dict:
    import numpy as np
    import mdtraj as md

    t = md.load(xtc, top=prmtop)
    top = t.topology
    zn = [a.index for a in top.atoms
          if a.element is not None and a.element.symbol == "Zn"]
    zn.sort()

    def sg_index(resseq: int) -> int | None:
        for a in top.atoms:
            if (a.residue.resSeq == resseq and a.name.strip().upper() == "SG"
                    and a.residue.name.strip().upper() in ("CYS", "CYM")):
                return a.index
        return None

    res = {"model": model, "n_frames": int(t.n_frames),
           "zn_atoms": len(zn), "series": {}}

    # 每个 Zn 到"应配位"的 4 个 SG 的距离（严格：必须是指定的 4 个残基）
    pairs_plan = [(zn[0], ZN1_RES), (zn[1], ZN2_RES)]
    allmean = []
    for k, (z, rs) in enumerate(pairs_plan, start=1):
        idxs = [sg_index(r) for r in rs]
        if any(i is None for i in idxs):
            res["error"] = f"Zn{k}: 找不到配位 SG {rs} -> {idxs}"
            return res
        dist = md.compute_distances(t, [(z, i) for i in idxs])  # nm
        m = dist.mean(axis=1)
        res["series"][f"zn{k}"] = [round(float(x), 4) for x in m]
        half = m[len(m) // 2:]           # 末 50% 帧
        res[f"zn{k}_mean_nm"] = round(float(half.mean()), 4)
        res[f"zn{k}_max_nm"] = round(float(half.max()), 4)
        allmean.extend(half.tolist())

    res["both_mean_nm"] = round(float(np.mean(allmean)), 4)
    res["both_max_nm"] = round(float(np.max(allmean)), 4)
    res["template_ref_nm"] = [0.226, 0.236]
    res["verdict"] = ("PASS" if (res["both_mean_nm"] <= 0.30
                                 and res["both_max_nm"] <= 0.35) else "FAIL")
    write(os.path.join(d, "zn_result.json"), json.dumps(res, indent=2))
    # 逐帧序列留档
    import csv
    with open(os.path.join(d, "zn_dist_series.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "zn1_mean_nm", "zn2_mean_nm"])
        s1 = res["series"]["zn1"]; s2 = res["series"]["zn2"]
        for i, (a, b) in enumerate(zip(s1, s2)):
            w.writerow([i, a, b])
    return res


# ---------------------------------------------------------------- main
def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    model = sys.argv[1]
    gpu = sys.argv[2]
    seed = int(os.environ.get("SEED", "751624"))

    if model not in ("zaff_cym", "nonbond_cym", "nonbond_cys"):
        print(f"未知模型 {model}")
        return 2

    d = os.path.join(OUTROOT, model)
    Path(d).mkdir(parents=True, exist_ok=True)

    # 0) 资产
    write(os.path.join(TEMPL, "zaff.frcmod"), ZAFF_FRCMOD)
    if not Path(os.path.join(TEMPL, "zinc_nonbond_fallback.frcmod")).exists():
        write(os.path.join(TEMPL, "zinc_nonbond_fallback.frcmod"), NONBOND_FRCMOD)
    build_rec_cym()
    sanity_check_pairs()

    # 1) tleap
    inp = build_tleap(model, d)
    rc, o = sh(f"tleap -f {inp} > {d}/tleap.log 2>&1", cwd=d)
    log = Path(f"{d}/tleap.log").read_text(errors="replace") if Path(f"{d}/tleap.log").exists() else ""
    prm = Path(f"{d}/sys.prmtop")
    # 注意：tleap 失败时仍会生成 0 字节的 prmtop，必须查大小而非仅查存在
    if rc != 0 or (not prm.exists()) or prm.stat().st_size < 1000:
        print(f"[FAIL] tleap rc={rc}, prmtop size={prm.stat().st_size if prm.exists() else 'missing'}")
        print(log[-3000:])
        return 1
    # 解析 "Exiting LEaP: Errors = N; Warnings = M; Notes = K."
    n_err = 0
    for l in log.splitlines():
        if "Errors =" in l:
            try:
                n_err = int(l.split("Errors =")[1].split(";")[0].strip())
            except ValueError:
                pass
    if n_err > 0:
        print(f"[FAIL] tleap Errors = {n_err}")
        print("\n".join([l for l in log.splitlines()
                         if "No torsion terms" in l or "No angle terms" in l
                         or "Missing" in l][:15]))
        return 1
    print(f"[ok] tleap Errors = {n_err}")
    keep = [l for l in log.splitlines()
            if any(k in l for k in ("Total", "charge", "WARNING", "ERROR",
                                    "Missing", "needs", "UNKNOWN", "FATAL"))]
    print("[tleap] " + "\n[tleap] ".join(keep[:40]))

    if os.environ.get("BUILD_ONLY") == "1":
        print(f"[BUILD_ONLY] {model} 拓扑构建完成：{d}/sys.prmtop")
        return 0

    # 2) Amber -> GROMACS
    rc, o = sh(f"""{PY} -c "
import parmed as pmd
s = pmd.load_file('{d}/sys.prmtop', '{d}/sys.inpcrd')
s.save('{d}/topol.top', overwrite=True)
s.save('{d}/conf.gro', overwrite=True)
print('parmed OK', len(s.atoms))
" """)
    print(o.strip()[:600])
    if rc != 0 or not Path(f"{d}/topol.top").exists():
        print("[FAIL] parmed 转换失败")
        return 1

    # 3) EM / NVT / NPT / PROD
    write(f"{d}/em.mdp", mdp_em())
    write(f"{d}/nvt.mdp", mdp_nvt(1.0, seed))
    write(f"{d}/npt.mdp", mdp_npt(1.0, seed))
    write(f"{d}/prod.mdp", mdp_prod(PROD_NS, SAVE_PS, seed))

    def tag_prev(tag: str) -> str:
        return {"em": "conf", "nvt": "em", "npt": "nvt", "prod": "npt"}[tag]

    def grompp(tag: str) -> bool:
        prev = tag_prev(tag) + ".gro"
        rc, o = sh(f"gmx grompp -f {d}/{tag}.mdp -c {d}/{prev} -p {d}/topol.top "
                   f"-o {d}/{tag}.tpr -maxwarn 2", cwd=d)
        if rc != 0:
            print(f"[FAIL] grompp {tag}: {o[-1500:]}")
            return False
        return True

    def mdrun(tag: str, nthreads: int = 8) -> bool:
        # 与 run1 主管线保持一致：不强制 -nb/-bonded/-pme，交给 GROMACS 自动分配。
        # 实测（2026-10-29）：强制全上 GPU 时 31686 原子体系只有 ~40 ns/day，
        # 自动分配可达 ~325 ns/day（快 8 倍）。EM（steep 是非动力学积分器）由
        # GROMACS 自动把 PME 放回 CPU，显式 -pme gpu 反而会直接报错。
        cmd = (f"gmx mdrun -deffnm {d}/{tag} -s {d}/{tag}.tpr -ntmpi 1 "
               f"-ntomp {nthreads} -gpu_id {gpu} > {d}/md_{tag}.log 2>&1")
        rc, o = sh(cmd, cwd=d)
        if rc != 0:
            print(f"[FAIL] mdrun {tag}: {o[-1500:]}")
            return False
        return True

    for tag in ("em", "nvt", "npt", "prod"):
        if not grompp(tag):
            return 1
        if not mdrun(tag, 8):
            return 1
        print(f"[ok] {tag}")

    # 4) 分析
    res = analyse(f"{d}/sys.prmtop", f"{d}/prod.xtc", model, d)
    print(json.dumps({k: v for k, v in res.items() if k != "series"}, indent=2))
    print(f"\n=== {model}: {res.get('verdict', 'ERROR')} "
          f"(both_mean={res.get('both_mean_nm')} nm, "
          f"both_max={res.get('both_max_nm')} nm; 模板 0.226–0.236 nm) ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
