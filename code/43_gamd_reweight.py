#!/usr/bin/env python3
"""43_gamd_reweight.py — GaMD 汇总 + 重加权可行性诊断（§3.9）。

输入（code/42_gamd_analyse.py pull 拉回 results/<TAG>/）：
  <job>__frames.csv   time_ps, lig_rmsd, lig_rmsd_if, if_drift, min_dist,
                      n_res_contact
                      · lig_rmsd     全蛋白 CA 叠合后的配体 RMSD（对照用）
                      · lig_rmsd_if  **界面残基 CA 叠合**后的配体 RMSD（主指标）
                      · if_drift     界面区相对全蛋白坐标系的漂移（蛋白构象变化）
  <job>__gamd-reweight.dat   step, V_P, V_D, ..., dV_P, dV_D（kcal/mol）
  <job>__resocc.csv          逐残基接触占据率（偏置系综）
  <job>__rescontact.csv      逐帧接触位图（供 exp(βΔV) 重加权试算）

**为什么主指标改成界面局部叠合**
E6 是双结构域蛋白，LxxLL 界面所在结构域相对全蛋白坐标系会整体漂移
（实测 848 rep1 达 18.6 Å）。此时"全蛋白 CA 叠合 + 配体 RMSD"把蛋白构象
变化算进了配体位移，系统性高估 4–5 Å（10.7 → 6.0 Å）。判断"姿态是否保持"
必须在结合位点局部坐标系里做，这是本脚本 r = lig_rmsd_if 的原因。

**重加权：只做诊断，不出 PMF**
GaMD 抬高势能 ΔV：P(x) ∝ P*(x)·exp(+βΔV)。二阶累积量与指数平均两种估计
在 ΔV 近似高斯时应一致。本项目 σ₀=0.25 双 boost 实测：
  · β·σ_ΔV ≈ 1.1–1.3（看似达标）
  · 但 ESS/N 仅 0.2–3.7%（ΔV 右尾重，指数平均被少数帧主导）
  · 累积量 vs EA 最大偏差达 1.4 kcal/mol
故 **PMF 定量结果不作报告**，只在报告中给出这三项诊断作为"不可重加权"的证据。
占据率同样标注为偏置系综帧比例（定性 engagement 指标），不作平衡布居解读。

输出（results/<TAG>/）：
  gamd_summary.csv          每副本汇总量（含界面/全局两套 RMSD）
  contact_occupancy.csv     逐残基占据率（未重加权 + 重加权试算）
  reweight_check.csv        每副本重加权可行性诊断
  gamd_report.md            人读摘要 + 可重加权性裁决
"""
from __future__ import annotations

import csv
import math
import os
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
TAG = os.environ.get("GAMD_TAG", "gamd").strip() or "gamd"
G = REBUILD / "results" / TAG      # 输入（42 拉回）

SYSTEMS = {
    "CHEMBL848_pose1": "CHEMBL848",
    "CHEMBL1129_pose1": "CHEMBL1129",
    "LxxLL_peptide": "LxxLL",
}
REPLICAS = [1, 2, 3]

T = 300.0
KB = 0.001987204259          # kcal/mol/K
KT = KB * T                  # 0.59616 kcal/mol
BETA = 1.0 / KT

DT = 0.002                   # ps
PROD_START_PS = 4600.0       # ntcmdprep+ntcmd+ntebprep+nteb = 2.3e6 步
BIN_RMSD = (0.0, 40.0, 1.0)  # Å（诊断分箱用）

# 可重加权性判据（预先定义，不得事后放宽）
ESS_FRAC_MIN = 0.15          # ESS/N ≥ 15%
BETA_SD_MAX = 1.5            # β·σ_ΔV ≤ 1.5
MAXDIFF_MAX = 1.0            # 累积量 vs EA ≤ 1 kcal/mol


def read_reweight(path: Path) -> dict[int, float]:
    """返回 step -> ΔV_total (kcal/mol)（双 boost 两列相加）。"""
    out = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 9:
                continue
            try:
                step = int(f[1])
                dvp = float(f[7])
                dvd = float(f[8])
            except ValueError:
                continue
            out[step] = dvp + dvd
    return out


def read_rescontact(path: Path):
    """读逐帧接触位图，返回 (resids, resnames, M[n_frames, n_res] uint8)。"""
    with path.open() as fh:
        lines = fh.read().splitlines()
    resids = [int(x) for x in lines[0].split(",")[1:]]
    resnames = lines[1].split(",")[1:]
    M = np.array([[int(c) for c in ln.split(",")[1:]] for ln in lines[2:]],
                 dtype=np.uint8)
    return resids, resnames, M


def read_frames(path: Path) -> dict[str, np.ndarray]:
    """读 frames.csv；缺 lig_rmsd_if 列时回退到 lig_rmsd（旧版产物）。"""
    cols = defaultdict(list)
    with path.open() as fh:
        for row in csv.DictReader(fh):
            for k in ("time_ps", "lig_rmsd", "lig_rmsd_if", "if_drift",
                      "min_dist", "n_res_contact"):
                if k in row and row[k] not in ("", None):
                    cols[k].append(float(row[k]))
    out = {k: np.array(v) for k, v in cols.items()}
    if "lig_rmsd_if" not in out:
        print(f"  [warn] {path.name} 无 lig_rmsd_if 列，回退用全蛋白叠合 RMSD")
        out["lig_rmsd_if"] = out["lig_rmsd"]
    if "if_drift" not in out:
        out["if_drift"] = np.full_like(out["lig_rmsd"], np.nan)
    return out


def pmf_bins(x: np.ndarray, dv: np.ndarray, binspec, min_count: int = 20):
    """二阶累积量 + EA 重加权 PMF（仅用于一致性诊断，不报告曲线）。"""
    lo, hi, w = binspec
    nb = int(round((hi - lo) / w))
    idx = np.clip(((x - lo) / w).astype(int), 0, nb - 1)
    centers = lo + (np.arange(nb) + 0.5) * w
    a_rw = np.full(nb, np.nan)
    a_ea = np.full(nb, np.nan)
    counts = np.zeros(nb, dtype=int)
    for b in range(nb):
        sel = idx == b
        n = int(sel.sum())
        counts[b] = n
        if n < min_count:
            continue
        d = dv[sel]
        a_rw[b] = -KT * math.log(n) - d.mean() - 0.5 * BETA * d.var()
        mx = d.max()
        a_ea[b] = -KT * (mx + math.log(np.exp(BETA * (d - mx)).mean()) + math.log(n))
    return centers, a_rw, a_ea, counts


def load_all():
    data = defaultdict(list)
    occ = defaultdict(lambda: defaultdict(list))
    resname = {}
    summary = []
    mats = {}
    for sysname, label in SYSTEMS.items():
        for rep in REPLICAS:
            job = f"{sysname}_rep{rep}"
            fp = G / f"{job}__frames.csv"
            rp = G / f"{job}__gamd-reweight.dat"
            op = G / f"{job}__resocc.csv"
            if not fp.exists() or not rp.exists():
                print(f"  [skip] {job} 缺少输入")
                continue
            fr = read_frames(fp)
            t = fr["time_ps"]
            n_frame_all = len(t)
            rw = read_reweight(rp)
            steps = np.round(t / DT).astype(np.int64)
            dv = np.array([rw.get(int(s), np.nan) for s in steps])
            ok = ~np.isnan(dv)
            cp = G / f"{job}__rescontact.csv"
            if cp.exists():
                try:
                    rids, rns, Mall = read_rescontact(cp)
                except Exception as ex:                      # pragma: no cover
                    print(f"  [warn] {job} rescontact 解析失败：{ex}")
                    rids, rns, Mall = None, None, None
                if Mall is not None and Mall.shape[0] == n_frame_all:
                    mats[job] = (rids, rns, Mall[ok].astype(np.float64))
                elif Mall is not None:
                    print(f"  [warn] {job} rescontact 行数 {Mall.shape[0]} "
                          f"!= frames {n_frame_all}，跳过重加权占据率")
            t = t[ok]
            r_if = fr["lig_rmsd_if"][ok]
            r_all = fr["lig_rmsd"][ok]
            drift = fr["if_drift"][ok]
            m = fr["min_dist"][ok]
            dv = dv[ok]
            if len(t) == 0:
                print(f"  [skip] {job} 无可用帧")
                continue
            data[label].append((job, t, r_if, m, dv))
            if op.exists():
                with op.open() as fh:
                    for row in csv.DictReader(fh):
                        key = int(row["resid"])
                        occ[label][key].append(float(row["occupancy"]))
                        resname[(label, key)] = row.get("resname", "")
            summary.append({
                "system": label, "replica": rep, "job": job,
                "n_frames": len(t),
                "prod_ns": (t.max() - t.min()) / 1000.0,
                # 主指标：界面局部叠合
                "mean_rmsd_if": float(r_if.mean()),
                "sd_rmsd_if": float(r_if.std(ddof=1)),
                "med_rmsd_if": float(np.median(r_if)),
                "frac_rmsd_if_le4": float((r_if <= 4.0).mean()),
                "frac_rmsd_if_gt40": float((r_if > 40.0).mean()),
                # 对照：全蛋白叠合（含蛋白构象变化）
                "mean_rmsd_allca": float(r_all.mean()),
                "mean_if_drift": float(np.nanmean(drift)),
                # 接触
                "mean_dmin": float(m.mean()), "sd_dmin": float(m.std(ddof=1)),
                "frac_dmin_le3.5": float((m <= 3.5).mean()),
                # 三态分解（Figure 5C，全部为帧计数，不涉及重加权）：
                #   retained  姿态保持：界面叠合 RMSD ≤ 4 Å 且仍在接触
                #   displaced 仍在接触但姿态已偏离（滑移/换 pose）
                #   detached  完全脱离（最小重原子距离 > 3.5 Å）
                "frac_retained": float(((r_if <= 4.0) & (m <= 3.5)).mean()),
                "frac_displaced": float(((r_if > 4.0) & (m <= 3.5)).mean()),
                "frac_detached": float((m > 3.5).mean()),
                # boost
                "mean_dV": float(dv.mean()),
                "sd_dV": float(dv.std(ddof=1)),
                "beta_sd_dV": float(dv.std(ddof=1) * BETA),
                "max_dV": float(dv.max()),
            })
    return data, occ, resname, summary, mats


def weighted_occupancy(data, mats):
    """exp(βΔV) 重加权占据率 + 每副本 Kish ESS（诊断用）。"""
    occw = defaultdict(lambda: defaultdict(list))
    resname = {}
    ess = {}
    for label, reps in data.items():
        for job, t, r, m, dv in reps:
            if job not in mats:
                continue
            rids, rns, M = mats[job]
            w = np.exp(BETA * (dv - dv.max()))
            w /= w.sum()
            ess[job] = float(1.0 / (w ** 2).sum())
            pw = (w[:, None] * M).sum(0)
            for i, rid in enumerate(rids):
                occw[label][rid].append(float(pw[i]))
                resname[(label, rid)] = rns[i]
    return occw, resname, ess


def consistency_pass(data):
    """二阶累积量 vs 指数平均 PMF 的最大偏差（沿界面叠合 RMSD）。"""
    out = {}
    for label, reps in data.items():
        for job, t, r, m, dv in reps:
            c, a_rw, a_ea, cnt = pmf_bins(r, dv, BIN_RMSD)
            ok = (cnt >= 20) & ~np.isnan(a_rw) & ~np.isnan(a_ea)
            if ok.sum() < 2:
                out[job] = float("nan")
                continue
            dr = a_rw[ok] - np.nanmin(a_rw[ok])
            de = a_ea[ok] - np.nanmin(a_ea[ok])
            out[job] = float(np.nanmax(np.abs(dr - de)))
    return out


def main() -> int:
    if not G.exists():
        print(f"缺少目录 {G}；先跑 code/42_gamd_analyse.py pull")
        return 1
    data, occ, resname, summary, mats = load_all()
    if not data:
        print("没有任何可用数据")
        return 1

    occw, resname_w, ess = weighted_occupancy(data, mats)
    if occw:
        resname.update(resname_w)

    cons = consistency_pass(data)
    for s in summary:
        e = ess.get(s["job"], float("nan"))
        s["ESS_frames"] = e
        s["ESS_frac"] = (e / s["n_frames"]) if s["n_frames"] else float("nan")
        s["maxdiff_cumulant_vs_EA_kcal"] = cons.get(s["job"], float("nan"))

    # ---- 每副本汇总 ----
    with (G / "gamd_summary.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)
    print(f"写入 gamd_summary.csv（{len(summary)} 个副本）")

    # ---- 重加权可行性诊断表 ----
    with (G / "reweight_check.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["system", "replica", "n_frames", "mean_dV", "sd_dV",
                    "beta_sd_dV", "ESS_frames", "ESS_frac",
                    "maxdiff_cumulant_vs_EA_kcal",
                    "pass_beta_sd", "pass_ess", "pass_maxdiff", "reweightable"])
        for s in summary:
            pb = s["beta_sd_dV"] <= BETA_SD_MAX
            pe = s["ESS_frac"] >= ESS_FRAC_MIN
            pm = s["maxdiff_cumulant_vs_EA_kcal"] <= MAXDIFF_MAX
            w.writerow([s["system"], s["replica"], s["n_frames"],
                        f"{s['mean_dV']:.4f}", f"{s['sd_dV']:.4f}",
                        f"{s['beta_sd_dV']:.4f}", f"{s['ESS_frames']:.1f}",
                        f"{s['ESS_frac']:.5f}",
                        f"{s['maxdiff_cumulant_vs_EA_kcal']:.4f}",
                        int(pb), int(pe), int(pm), int(pb and pe and pm)])
    print("写入 reweight_check.csv")

    # ---- 逐残基占据率（偏置系综 + 重加权试算）----
    with (G / "contact_occupancy.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["system", "resid", "resname", "occ_mean", "occ_sd",
                    "occ_rw_mean", "occ_rw_sd", "n_rep"])
        for label in sorted(occ):
            for rid in sorted(occ[label]):
                v = np.array(occ[label][rid])
                uw = np.array(occw[label][rid]) if occw and rid in occw[label] else None
                sd = "" if len(v) < 2 else f"{v.std(ddof=1):.6f}"
                if uw is None or len(uw) == 0:
                    mw, sw = "", ""
                else:
                    mw = f"{uw.mean():.6f}"
                    sw = "" if len(uw) < 2 else f"{uw.std(ddof=1):.6f}"
                w.writerow([label, rid, resname.get((label, rid), ""),
                            f"{v.mean():.6f}", sd, mw, sw, len(v)])
    print("写入 contact_occupancy.csv")

    # ---- 人读摘要 ----
    worst_beta = max(s["beta_sd_dV"] for s in summary)
    worst_ess = max(s["ESS_frac"] for s in summary)
    worst_diff = max(s["maxdiff_cumulant_vs_EA_kcal"] for s in summary)
    n_pass = sum(1 for s in summary
                 if s["beta_sd_dV"] <= BETA_SD_MAX
                 and s["ESS_frac"] >= ESS_FRAC_MIN
                 and s["maxdiff_cumulant_vs_EA_kcal"] <= MAXDIFF_MAX)
    rew_ok = n_pass == len(summary)
    lines = [
        "# GaMD §3.9 summary (generated)", "",
        "## 重加权可行性裁决",
        f"- 判据：β·σ_ΔV ≤ {BETA_SD_MAX} **且** ESS/N ≥ {ESS_FRAC_MIN:.0%} "
        f"**且** |累积量 − EA| ≤ {MAXDIFF_MAX} kcal/mol",
        f"- 实测：max β·σ = {worst_beta:.2f}；max ESS/N = {worst_ess:.2%}；"
        f"max |累积量 − EA| = {worst_diff:.2f} kcal/mol",
        f"- **通过副本数 {n_pass}/{len(summary)} → "
        f"{'可重加权' if rew_ok else '不可重加权：不报告定量 PMF'}**", "",
        "## 逐副本汇总量",
        "",
    ]
    cols = [("system", "{}"), ("replica", "{}"), ("n_frames", "{}"),
            ("prod_ns", "{:.1f}"),
            ("mean_rmsd_if", "{:.2f}"), ("frac_rmsd_if_le4", "{:.3f}"),
            ("frac_rmsd_if_gt40", "{:.3f}"),
            ("mean_rmsd_allca", "{:.2f}"), ("mean_if_drift", "{:.2f}"),
            ("mean_dmin", "{:.2f}"), ("frac_dmin_le3.5", "{:.3f}"),
            ("beta_sd_dV", "{:.2f}"), ("ESS_frac", "{:.3f}"),
            ("maxdiff_cumulant_vs_EA_kcal", "{:.2f}")]
    cols = [c for c in cols if c[0] in summary[0]]
    lines.append("| " + " | ".join(c[0] for c in cols) + " |")
    lines.append("|" + "---|" * len(cols))
    for s in summary:
        lines.append("| " + " | ".join(f.format(s[k]) for k, f in cols) + " |")
    (G / "gamd_report.md").write_text("\n".join(lines) + "\n")
    print(f"写入 gamd_report.md（可重加权={'是' if rew_ok else '否'}）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
