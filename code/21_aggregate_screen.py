#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""21_aggregate_screen.py — 全量筛选结果汇总（MM-PBSA + 姿态稳定性 + Zn QC）

用法
----
    python 21_aggregate_screen.py [结果目录] [-o 输出目录]

默认结果目录： ../results_fullscreen_20ns
输出：
    results/21_screen_summary.csv      逐 (配体, replica) 明细
    results/21_screen_byligand.csv     逐配体汇总（含跨副本均值±SD、保留帧率、Zn 合格率）
    results/21_screen_report.md         可读报告（含 PASS/FAIL 判据）

判据（与管线一致，可预注册）
--------------------------
  - qc_pass=1 需同时满足：保留帧 ≥ MIN_KEEP_FRAMES(500) 且 双Zn四配位合格率 ≥ 0.95
  - 配体级结论：「3/3 副本 QC 通过」= STABLE；「1–2 副本通过」= MARGINAL；「0 副本」= UNSTABLE
  - ΔH 仅在 qc_pass=1 的副本上取均值；UNSTABLE 配体不给 ΔH（记为 NA）
"""
from __future__ import annotations

import csv
import glob
import os
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
DEFAULT_SRC = REBUILD / "results_fullscreen_20ns"

FIELDS_NUM = ("delta_H_PB_kcal_mol", "dG_gas", "dG_solv", "zn_pass_rate",
              "zn_mean_min4_dist", "n_frames")


def fnum(v):
    try:
        if v in ("", "NA", None):
            return None
        return float(v)
    except Exception:
        return None


def load_rows(src: Path):
    rows = []
    for f in sorted(glob.glob(str(src / "SUMMARY_v3_*_rep*.csv"))):
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh):
                if not r.get("CHEMBL"):
                    continue
                rr = dict(r)
                for k in FIELDS_NUM:
                    rr[k + "__n"] = fnum(r.get(k))
                # 从 notes 里回收 kept/N 帧数
                kept = None
                ntot = None
                for tok in (r.get("notes") or "").split("|"):
                    tok = tok.strip()
                    if tok.startswith("kept "):
                        try:
                            a, b = tok.split()[1].split("/")
                            kept, ntot = int(a), int(b)
                        except Exception:
                            pass
                rr["kept_frames"] = kept
                rr["total_frames"] = ntot
                # 单位显式化：管线写出的 zn_mean_min4_dist 是 nm（mdtraj 单位）。
                # 旧批次（阈值误为 28 Å）该值 0.825–2.372 nm（中位 0.997）；模板参考 0.226–0.236 nm。
                _zn = rr.get("zn_mean_min4_dist__n")
                rr["zn_mean_min4_dist_A"] = (round(_zn * 10.0, 2)
                                             if _zn is not None else "NA")
                rr["kept_frac"] = (kept / ntot) if (kept and ntot) else None
                rows.append(rr)
    return rows


def main() -> int:
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_SRC
    out = REBUILD / "results"
    out.mkdir(parents=True, exist_ok=True)
    rows = load_rows(src)
    if not rows:
        print(f"未在 {src} 找到 SUMMARY_v3_*_rep*.csv")
        return 1

    # 防御性自检：Zn 门是否真的在起作用（单位 bug 曾使阈值变成 28 Å 而永不触发）
    znm = [r["zn_mean_min4_dist__n"] for r in rows
           if r.get("zn_mean_min4_dist__n") is not None]
    if znm:
        mzn = statistics.fmean(znm)
        if mzn > 0.5:
            print(f"[WARN] Zn–S(最近4) 均距均值 {mzn:.3f} nm ({mzn*10:.1f} Å)，"
                  f"远超晶体模板 0.226–0.236 nm → 锌指未被四配位。"
                  f"请检查 ZN_S_CUTOFF 单位(nm) 与 Zn 参数模型；本批次结合结果不可解释。")
        else:
            print(f"[OK] Zn–S(最近4) 均距均值 {mzn:.3f} nm ({mzn*10:.2f} Å)，"
                  f"落在四配位合理区间。")

    # 逐 (配体, replica) 明细
    det = out / "21_screen_summary.csv"
    cols = ["CHEMBL", "replica", "seed", "status", "qc_pass", "delta_H_PB_kcal_mol",
            "dG_gas", "dG_solv", "kept_frames", "total_frames", "kept_frac",
            "zn_pass_rate", "zn_mean_min4_dist", "zn_mean_min4_dist_A",
            "prod_ns", "minutes", "notes"]
    with open(det, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # 逐配体汇总
    bylig = {}
    for r in rows:
        bylig.setdefault(r["CHEMBL"], []).append(r)
    agg = []
    for lig, rs in sorted(bylig.items()):
        ok = [r for r in rs if r.get("qc_pass") == "1"]
        dh = [r["delta_H_PB_kcal_mol__n"] for r in ok
              if r["delta_H_PB_kcal_mol__n"] is not None]
        kf = [r["kept_frac"] for r in rs if r["kept_frac"] is not None]
        zn = [r["zn_pass_rate__n"] for r in rs if r["zn_pass_rate__n"] is not None]
        n_rep = len(rs)
        n_ok = len(ok)
        verdict = ("STABLE" if n_ok == n_rep and n_rep >= 3
                   else "MARGINAL" if n_ok >= 1 else "UNSTABLE")
        agg.append({
            "CHEMBL": lig,
            "n_replicas": n_rep,
            "n_qc_pass": n_ok,
            "verdict": verdict,
            "delta_H_mean": round(statistics.fmean(dh), 3) if dh else "NA",
            "delta_H_sd": round(statistics.stdev(dh), 3) if len(dh) > 1 else "NA",
            "delta_H_values": ";".join(f"{x:.2f}" for x in dh) if dh else "NA",
            "kept_frac_mean": round(statistics.fmean(kf), 3) if kf else "NA",
            "kept_frac_values": ";".join(f"{x:.3f}" for x in kf) if kf else "NA",
            "zn_pass_rate_min": round(min(zn), 3) if zn else "NA",
            "status_set": ";".join(sorted({r["status"] for r in rs})),
        })
    fa = out / "21_screen_byligand.csv"
    with open(fa, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(agg[0].keys()))
        w.writeheader()
        for a in agg:
            w.writerow(a)

    # 可读报告
    md = out / "21_screen_report.md"
    n_stable = sum(1 for a in agg if a["verdict"] == "STABLE")
    n_marg = sum(1 for a in agg if a["verdict"] == "MARGINAL")
    n_un = sum(1 for a in agg if a["verdict"] == "UNSTABLE")
    with open(md, "w") as fh:
        fh.write("# 全量筛选汇总（14 配体 × 3 replica × 20 ns）\n\n")
        fh.write(f"- 明细行数 {len(rows)}，配体 {len(agg)} 个\n")
        fh.write(f"- STABLE（3/3 副本 QC 通过）**{n_stable}**；"
                 f"MARGINAL（1–2/3）**{n_marg}**；UNSTABLE（0/3）**{n_un}**\n\n")
        fh.write("| 配体 | 副本 | QC通过 | 判定 | ΔH 均值±SD (kcal/mol) | 保留帧率均值 | Zn合格率(最低) | 状态 |\n")
        fh.write("|---|---|---|---|---|---|---|---|\n")
        for a in agg:
            dh = (f"{a['delta_H_mean']}±{a['delta_H_sd']}"
                  if a["delta_H_mean"] != "NA" else "NA")
            fh.write(f"| {a['CHEMBL']} | {a['n_replicas']} | {a['n_qc_pass']} | "
                     f"{a['verdict']} | {dh} | {a['kept_frac_mean']} | "
                     f"{a['zn_pass_rate_min']} | {a['status_set']} |\n")
        fh.write("\n## 判据\n")
        fh.write("- 副本级 qc_pass=1：保留帧 ≥500（配体 RMSD<4 Å）且 双Zn四配位合格率 ≥0.95\n")
        fh.write("- 配体级：3/3 通过=STABLE；1–2/3=MARGINAL；0/3=UNSTABLE（不报 ΔH）\n")

    print(f"明细 {len(rows)} 行 -> {det}")
    print(f"逐配体 {len(agg)} 行 -> {fa}")
    print(f"报告 -> {md}")
    print(f"\nSTABLE {n_stable} / MARGINAL {n_marg} / UNSTABLE {n_un}")
    for a in agg:
        dh = (f"{a['delta_H_mean']}±{a['delta_H_sd']}"
              if a["delta_H_mean"] != "NA" else "NA")
        print(f"  {a['CHEMBL']:20s} {a['verdict']:9s} "
              f"QC {a['n_qc_pass']}/{a['n_replicas']}  ΔH={dh}  保留帧率={a['kept_frac_mean']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
