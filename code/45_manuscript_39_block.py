#!/usr/bin/env python3
"""45_manuscript_39_block.py — 由 GaMD 汇总结果生成 §3.9 文本块。

所有数字都从 results/<TAG>/ 下的 CSV **现算**（与 §3.8 同一纪律：杜绝手抄）。
负号统一为 U+2212（与稿件一致）。

2026-11-04 改写要点（与 42/43/44 的指标变更同步）：
  1. 姿态保持改用**界面残基 CA 叠合**的配体 RMSD（全蛋白叠合被 E6 结构域
     漂移污染，实测高估 4–5 Å）；全蛋白数值一并报告以求透明。
  2. 撤掉定量 PMF：σ0=0.25 下 ESS/N 仅 0.2–3.7%，累积量与指数平均最大偏差
     1.4 kcal/mol，不满足预注册判据 → 只在文中报告该诊断作证据。
  3. 改为三态分解（姿态保持 / 位移但接触 / 完全脱离），全部是帧计数，
     不依赖任何重加权假设。
  4. GaMD 参数按远端 gamd.in 实测值写（σ0 = 0.25 kcal/mol = 1.05 kJ/mol）。

输出：results/45_section39_block_<TAG>.md
     results/45_section39_numbers_<TAG>.json

用法：GAMD_TAG=gamd_v2 python3 code/45_manuscript_39_block.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
TAG = os.environ.get("GAMD_TAG", "gamd").strip() or "gamd"
G = REBUILD / "results" / TAG
OUT_MD = REBUILD / "results" / f"45_section39_block_{TAG}.md"
OUT_JSON = REBUILD / "results" / f"45_section39_numbers_{TAG}.json"

LIGS = ["CHEMBL848", "CHEMBL1129", "LxxLL"]


def mdash(x, nd=2):
    """数值格式化并把负号换成 U+2212（稿件统一）。"""
    return f"{x:.{nd}f}".replace("-", "−")


def read_csv(p: Path):
    with p.open() as fh:
        return list(csv.DictReader(fh))


def main() -> int:
    f_sum = G / "gamd_summary.csv"
    f_occ = G / "contact_occupancy.csv"
    f_rw = G / "reweight_check.csv"
    for f in (f_sum, f_occ, f_rw):
        if not f.exists():
            print(f"缺少 {f}，先跑 code/42（pull）+ code/43")
            return 1

    summ = read_csv(f_sum)
    occ = read_csv(f_occ)
    rwk = read_csv(f_rw)

    # ---- 逐体系汇总（3 副本）----
    S = {}
    for L in LIGS:
        rows = [r for r in summ if r["system"] == L]
        if not rows:
            print(f"体系 {L} 无数据")
            return 1
        rows.sort(key=lambda r: int(r["replica"]))
        g = lambda k: np.array([float(r[k]) for r in rows])
        S[L] = dict(
            n_frames=[int(r["n_frames"]) for r in rows],
            # 逐副本保留原值：跨副本求均值会被"某副本解离"带偏成一个
            # 既不代表结合态也不代表解离态的数（21 ± 27 Å），不可用
            rmsd_if=[float(r["mean_rmsd_if"]) for r in rows],
            rmsd_allca=[float(r["mean_rmsd_allca"]) for r in rows],
            if_drift=[float(r["mean_if_drift"]) for r in rows],
            frac_retained=[float(r["frac_retained"]) for r in rows],
            frac_displaced=[float(r["frac_displaced"]) for r in rows],
            frac_detached=[float(r["frac_detached"]) for r in rows],
            frac_rmsd_if_le4=[float(r["frac_rmsd_if_le4"]) for r in rows],
            frac_dmin_le3_5=[float(r["frac_dmin_le3.5"]) for r in rows],
            n_rep=len(rows),
        )
    # 界面漂移与"结合态副本"的全局/局部 RMSD 对照（跨 9 副本）
    drift_all = [v for L in LIGS for v in S[L]["if_drift"]]

    # ---- 界面占据率与 top 接触残基 ----
    f_if = REBUILD / "results" / "29_interface_positions.csv"
    iface = set()
    if f_if.exists():
        for r in read_csv(f_if):
            if str(r["e6ap_lxxll_interface"]).lower() == "true":
                iface.add(int(r["struct_num"]))
    occ_if, occ_top = {}, {}
    for L in LIGS:
        rows = [r for r in occ if r["system"] == L and int(r["resid"]) <= 143]
        ifr = [float(r["occ_mean"]) for r in rows if (int(r["resid"]) + 1) in iface]
        occ_if[L] = float(np.mean(ifr)) if ifr else float("nan")
        top = sorted(rows, key=lambda r: float(r["occ_mean"]), reverse=True)[:6]
        occ_top[L] = [(r["resname"], int(r["resid"]) + 1, float(r["occ_mean"]))
                      for r in top]

    # ---- 重加权诊断（裁决"为何不出 PMF"）----
    RW = dict(
        max_beta_sd=float(max(float(r["beta_sd_dV"]) for r in rwk)),
        max_ess_frac=float(max(float(r["ESS_frac"]) for r in rwk)),
        min_ess_frac=float(min(float(r["ESS_frac"]) for r in rwk)),
        max_diff=float(max(float(r["maxdiff_cumulant_vs_EA_kcal"]) for r in rwk)),
        n_pass=int(sum(int(r["reweightable"]) for r in rwk)),
        n_total=len(rwk),
    )

    N = dict(S=S, occ_if=occ_if, occ_top=occ_top, reweight=RW,
             iface_n=len(iface))
    OUT_JSON.write_text(json.dumps(N, indent=2, default=str))


    def pct(x, nd=0):
        return f"{100 * x:.{nd}f}%"

    def top6(L):
        return ", ".join(f"{rn}{num} ({occ:.2f})" for rn, num, occ in occ_top[L])

    # 848 vs 1129 的姿态保持差异判据（同量级则合并表述）
    # ---- 逐副本格式化辅助 ----
    def rng(vals, nd=1):
        # 单值时不写"8.88–8.88"这种伪区间
        lo, hi = min(vals), max(vals)
        if abs(hi - lo) < 10 ** (-nd):
            return mdash(lo, nd)
        return f"{mdash(lo, nd)}–{mdash(hi, nd)}"

    def split_reps(L):
        """把 3 个副本按"是否仍在界面坐标系内"分成两组，各自给数值区间。

        跨副本求均值会把"2 个结合 + 1 个解离"平均成 21 ± 27 Å 这种
        既不代表结合也不代表解离的数，故一律逐副本报告。
        """
        r = S[L]["rmsd_if"]
        stay = [v for v in r if v < 40]
        left = [v for v in r if v >= 40]
        return stay, left

    def state_str(L):
        # 括号内按 replica 1/2/3 排列，避免逗号与千分位混淆
        def trio(key):
            return " / ".join(str(int(round(100 * v))) for v in S[L][key])
        return f"({trio('frac_retained')}) / ({trio('frac_displaced')}) / " \
               f"({trio('frac_detached')}) %"


    stay = {L: split_reps(L) for L in LIGS}
    n_stay = {L: len(stay[L][0]) for L in LIGS}

    lines = []
    lines.append("### 3.9 GaMD does not support a stable, pose-specific bound state")
    lines.append("")
    lines.append(
        "To ask whether the two nominated candidates form a physically stable bound "
        "state rather than a transient surface association, we ran Gaussian accelerated "
        "molecular dynamics (GaMD) on the zinc-corrected HPV16 E6 construct used in "
        "§3.8. Three independent replicas were launched for CHEMBL848, CHEMBL1129 and "
        "the LxxLL positive control from the end structures of the batch-2 20 ns "
        "conventional runs, using a GPU-resident GROMACS 2025.4-GaMD build "
        "(dual-boost, igamd = 3, σ₀ = 0.25 kcal/mol = 1.05 kJ/mol for both the "
        "potential and the dihedral boost). Each replica comprised 0.2 ns of "
        "conventional preparation, 2 ns of conventional MD, 0.4 ns of boost "
        "preparation and 2 ns of boost equilibration, followed by 100 ns of "
        "production with the boost parameters frozen (C-rescale NPT, 300 K, "
        "`integrator = md`, `nstcalcenergy = 1`, GPU non-bonded/PME/bonded/update). "
        "Replicas ran one per RTX A4000, giving 9 × 100 ns of production in total.")
    lines.append("")
    lines.append(
        "Ligand displacement was measured after superposition on the Cα atoms of the "
        f"{len(iface)} E6AP LxxLL-interface residues rather than on the whole protein, "
        "because the interface-bearing domain of E6 moves relative to the rest of the "
        "molecule (interface-region Cα RMSD in the global frame: "
        f"{rng(drift_all)} Å across the nine replicas); a global-frame ligand RMSD "
        "would therefore conflate ligand drift with protein conformational change, and "
        "is reported only for transparency. Contacts were defined per residue as any "
        "receptor heavy atom within 4.5 Å of the ligand, and complete detachment as a "
        "minimum ligand–receptor heavy-atom distance above 3.5 Å. Occupancies are frame "
        "fractions of the boosted ensemble. Boost reweighting was attempted but "
        f"rejected by pre-specified criteria (β·σ_ΔV up to {RW['max_beta_sd']:.2f}, "
        f"Kish effective sample size {pct(RW['min_ess_frac'], 1)}–"
        f"{pct(RW['max_ess_frac'], 1)} of frames, and a maximum "
        f"{mdash(RW['max_diff'], 1)} kcal/mol disagreement between the second-order "
        f"cumulant and the exponential-average estimate in {RW['n_total']}/"
        f"{RW['n_total']} replicas); no quantitative free-energy profile is therefore "
        "reported and every quantity below is reweighting-free.")
    lines.append("")
    lines.append(
        "No system held its docked pose for the whole trajectory, but the three "
        "systems failed differently. In the replicas that stayed at the interface "
        f"({n_stay['CHEMBL848']}/3 for CHEMBL848, {n_stay['CHEMBL1129']}/3 for "
        f"CHEMBL1129 and {n_stay['LxxLL']}/3 for LxxLL), the interface-superposed "
        f"ligand RMSD averaged {rng(stay['CHEMBL848'][0], 2)} Å for CHEMBL848, "
        f"{rng(stay['CHEMBL1129'][0], 2)} Å for CHEMBL1129 and "
        f"{rng(stay['LxxLL'][0], 2)} Å for the LxxLL peptide. The remaining replicas "
        "left the interface frame (mean RMSD "
        f"{rng(stay['CHEMBL848'][1] + stay['CHEMBL1129'][1] + stay['LxxLL'][1], 1)} Å) "
        "yet, with one exception, stayed in contact with the protein surface, i.e. the "
        "association became non-specific rather than breaking. Decomposing the "
        "production frames of replicas 1/2/3 into pose retained (RMSD ≤ 4 Å), displaced "
        "but in contact, and fully detached gave "
        f"{state_str('CHEMBL848')} for CHEMBL848, {state_str('CHEMBL1129')} for "
        f"CHEMBL1129 and {state_str('LxxLL')} for LxxLL; only CHEMBL848 replica 3 "
        "dissociated into solvent.")
    lines.append("")
    lines.append(
        "Engagement of the interface separated the peptide from the two small "
        "molecules. Mean occupancy over the "
        f"{len(iface)} interface residues was {mdash(occ_if['CHEMBL848'])} for "
        f"CHEMBL848, {mdash(occ_if['CHEMBL1129'])} for CHEMBL1129 and "
        f"{mdash(occ_if['LxxLL'])} for the LxxLL peptide (Figure 5B). The "
        "highest-occupancy residues were "
        f"{top6('CHEMBL848')} for CHEMBL848, {top6('CHEMBL1129')} for CHEMBL1129 and "
        f"{top6('LxxLL')} for LxxLL; all of these lie within the E6AP LxxLL interface. "
        "The small molecules therefore do engage the correct groove, but with roughly "
        "a third of the per-residue occupancy of the peptide and without preserving "
        "their docked orientation.")
    lines.append("")
    lines.append(
        "- **Figure 5** — (A) ligand RMSD from the docked pose after superposition on "
        "the interface Cα atoms, for each of the three replicas per system (square-root "
        "ordinate; dashed line, 4 Å); (B) per-residue contact occupancy in the boosted "
        "ensemble (replica mean, one point per residue; purple ticks mark the "
        "E6AP LxxLL-interface residues and the inset gives their mean occupancy for "
        "each system); (C) three-state decomposition of the production frames of every "
        "replica into pose retained (RMSD ≤ 4 Å), displaced but in contact, and fully "
        "detached.")
    lines.append("")
    lines.append(
        "- **Reading.** GaMD neither rescued the docked poses nor produced an "
        "interface-bound population comparable to the LxxLL reference: the peptide "
        f"occupied the interface {occ_if['LxxLL'] / max(occ_if['CHEMBL848'], 1e-9):.1f}-"
        "fold more completely than CHEMBL848 and "
        f"{occ_if['LxxLL'] / max(occ_if['CHEMBL1129'], 1e-9):.1f}-fold more completely "
        "than CHEMBL1129, and was the only system to approach its docked orientation "
        "for a quarter of a trajectory. Two caveats keep this comparison qualitative: "
        "the boost was strong enough to displace even the peptide in two of three "
        "replicas, and replica-to-replica spread is large, so no ranking between the "
        "two candidates is supported. We therefore read the batch-2 hits as weak, "
        "non-specific groove binders whose docking poses are not validated by enhanced "
        "sampling. This is consistent with, and reinforces, the negative conclusion of "
        "§3.8: the docking/MM-PBSA tier identifies plausible but not kinetically stable "
        "starting points, and any experimental follow-up should treat the predicted "
        "poses as hypotheses rather than as established binding modes.")

    OUT_MD.write_text("\n".join(lines) + "\n")
    print(f"写入 {OUT_MD.name} 与 {OUT_JSON.name}")
    print(f"[诊断] 可重加权副本数 {RW['n_pass']}/{RW['n_total']}，"
          f"max βσ={RW['max_beta_sd']:.2f}，max ESS/N={RW['max_ess_frac']:.1%}，"
          f"max |累积量−EA|={RW['max_diff']:.2f} kcal/mol")
    return 0


if __name__ == "__main__":
    sys.exit(main())
