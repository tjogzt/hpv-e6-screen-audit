#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""36_pull_full_screen.py — 拉取 batch 2 全量筛选结果并用对照校准判据汇总

用法
----
    python code/36_pull_full_screen.py            # 拉取 + 汇总
    python code/36_pull_full_screen.py --local    # 只汇总已下载到本地的结果

做什么
------
1) SFTP 下载远端 BASE 下全部 SUMMARY_v3_*.csv（14 配体 ×3 replica
   + DMSO ×3 + LxxLL ×3 = 48）到 results/batch2/
2) 解析逐 run 明细：ΔH_PB、保留帧率、Zn 四配位合格率
3) 用阴性对照 DMSO 与阳性对照 LxxLL 肽校准「结合稳定性」分档
4) 输出：
    results/36_screen_full.csv      逐 run 明细（含对照校准档位）
    results/36_screen_byligand.csv  逐配体汇总
    results/36_screen_report.md     可读报告（含对照分离带）

为什么需要对照校准
------------------
管线自带的硬判据「保留帧 ≥500 且 Zn 合格率 ≥0.95」是任意切点：
阳性对照 LxxLL 肽 rep1 只保留 215/1001 帧也会被判 FAIL。
因此本脚本保留原判据作为「严格档」，另以对照实测值给出经验分档：
    阴性上限 NEG_MAX = max(DMSO 保留帧率)
    阳性下限 POS_MIN = min(LxxLL 保留帧率)
    kept_frac <= max(NEG_MAX * 2, 0.05)   -> DMSO-like（与阴性对照不可区分）
    kept_frac >= POS_MIN                  -> LxxLL-like（达阳性对照量级）
    介于两者之间                           -> intermediate
规范：真实数据、无硬编码结论；路径 file.path 等价的 os.path.join；失败优雅降级。
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
OUT_DIR = REBUILD / "results"
PULL_DIR = OUT_DIR / "batch2"

HOST = os.environ.get("MATPOOL_HOST", "hz-4.matpool.com")
PORT = int(os.environ.get("MATPOOL_PORT", "26255"))
USER = os.environ.get("MATPOOL_USER", "root")
PASS = os.environ.get("MATPOOL_PASS", "wmKG=3YRsZaHmij0")
BASE = "/root/autodl-tmp/HPV_structure/candidates_run"

N_EXPECTED = 48          # 14 配体 ×3 + DMSO ×3 + LxxLL ×3
NEG_CTRL = "DMSO_pose1"
POS_CTRL = "LxxLL_peptide"

FIELDS_NUM = ("delta_H_PB_kcal_mol", "dG_gas", "dG_solv",
              "zn_pass_rate", "zn_mean_min4_dist", "n_frames")


def fnum(v):
    """NA/空串 -> None；其余转 float。"""
    try:
        if v in ("", "NA", None):
            return None
        return float(v)
    except Exception:
        return None


def pull() -> int:
    """SFTP 下载远端 SUMMARY_v3_*.csv；返回下载文件数（失败返回 -1）。"""
    try:
        import paramiko
    except ImportError:
        print("[FAIL] 缺少 paramiko，无法拉取；可用 --local 只汇总本地已有结果")
        return -1
    PULL_DIR.mkdir(parents=True, exist_ok=True)
    cli = paramiko.SSHClient()
    cli.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        cli.connect(HOST, port=PORT, username=USER, password=PASS, timeout=25)
    except Exception as e:                                   # 优雅降级
        print(f"[FAIL] SSH 连接失败：{e}")
        print("       检查实例是否开机、端口是否变化；"
              "可用 MATPOOL_HOST/MATPOOL_PORT/MATPOOL_PASS 覆盖。")
        return -1
    sftp = cli.open_sftp()
    try:
        names = [n for n in sftp.listdir(BASE)
                 if n.startswith("SUMMARY_v3_") and n.endswith(".csv")]
    except Exception as e:
        print(f"[FAIL] 列目录失败：{e}")
        sftp.close()
        cli.close()
        return -1
    n = 0
    for nm in sorted(names):
        dst = PULL_DIR / nm
        try:
            sftp.get(os.path.join(BASE, nm), str(dst))
            n += 1
        except Exception as e:
            print(f"[WARN] 下载 {nm} 失败：{e}")
    sftp.close()
    cli.close()
    print(f"[pull] 下载 {n}/{len(names)} 个 SUMMARY -> {PULL_DIR}")
    return n


def load_rows(src: Path):
    """读取 SUMMARY_v3_*_rep*.csv。

    去重是**必须的**：管线收尾会额外写一个合并文件 SUMMARY_v3_all.csv
    （含全部 48 行），若一并读入，每条 run 会被统计两次（曾使跨副本计数失真）。
    这里按 (CHEMBL, replica) 去重，并对被跳过的重复行给出告警。
    """
    rows = []
    seen: dict[tuple, dict] = {}
    n_dup = 0
    for f in sorted(glob.glob(str(src / "SUMMARY_v3_*_rep*.csv"))):
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh):
                if not r.get("CHEMBL"):
                    continue
                rr = dict(r)
                for k in FIELDS_NUM:
                    rr[k + "__n"] = fnum(r.get(k))
                # 从 notes 回收 kept/N 帧数（管线未单列）
                kept = ntot = None
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
                _zn = rr.get("zn_mean_min4_dist__n")
                rr["zn_mean_min4_dist_A"] = (round(_zn * 10.0, 2)
                                             if _zn is not None else "NA")
                rr["kept_frac"] = (kept / ntot) if (kept and ntot) else None
                key = (rr.get("CHEMBL"), str(rr.get("replica")))
                if key in seen:
                    n_dup += 1
                    continue
                seen[key] = rr
                rows.append(rr)
    if n_dup:
        print(f"[dedupe] 跳过 {n_dup} 行重复（SUMMARY_v3_all.csv 等合并文件）")
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--local", action="store_true", help="跳过拉取，只汇总本地")
    args = ap.parse_args()

    if not args.local:
        pull()

    rows = load_rows(PULL_DIR)
    if not rows:
        print(f"未在 {PULL_DIR} 找到 SUMMARY_v3_*_rep*.csv")
        return 1
    print(f"[load] {len(rows)} 行明细"
          f"{'' if len(rows) == N_EXPECTED else f'（未完成，应为 {N_EXPECTED}）'}")

    # ---- 防御性自检 1：Zn 门是否真的在起作用 ----
    znm = [r["zn_mean_min4_dist__n"] for r in rows
           if r.get("zn_mean_min4_dist__n") is not None]
    if znm:
        mzn = statistics.fmean(znm)
        tag = ("[OK]" if mzn <= 0.5 else "[WARN]")
        print(f"{tag} Zn–S(最近4) 均距 {mzn:.3f} nm ({mzn * 10:.2f} Å)；"
              f"晶体模板参考 0.226–0.236 nm")

    # ---- 防御性自检 2：对照是否存在且分离 ----
    neg = [r for r in rows if r["CHEMBL"] == NEG_CTRL and r["kept_frac"] is not None]
    pos = [r for r in rows if r["CHEMBL"] == POS_CTRL and r["kept_frac"] is not None]
    if not neg or not pos:
        print(f"[WARN] 对照缺失（DMSO {len(neg)} 行 / LxxLL {len(pos)} 行）"
              f"→ 跳过对照校准，仅用管线硬判据")
        neg_max = pos_min = None
    else:
        neg_max = max(r["kept_frac"] for r in neg)
        pos_min = min(r["kept_frac"] for r in pos)
        if neg_max >= pos_min:
            print(f"[WARN] 对照未分离：DMSO 上限 {neg_max:.3f} ≥ LxxLL 下限 "
                  f"{pos_min:.3f} → 本批次判据不可用，需排查")
        else:
            print(f"[OK] 对照分离：DMSO ≤{neg_max:.3f} < LxxLL ≥{pos_min:.3f}"
                  f"（空档 {pos_min - neg_max:.3f}）")

    # ---- 对照校准档位 ----
    neg_cut = max(neg_max * 2.0, 0.05) if neg_max is not None else None
    for r in rows:
        kf = r["kept_frac"]
        if kf is None:
            r["ctrl_tier"] = "NA"
        elif neg_cut is None:
            r["ctrl_tier"] = "uncalibrated"
        elif kf <= neg_cut:
            r["ctrl_tier"] = "DMSO-like"
        elif pos_min is not None and kf >= pos_min:
            r["ctrl_tier"] = "LxxLL-like"
        else:
            r["ctrl_tier"] = "intermediate"

    # ---- 逐 run 明细 ----
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    det = OUT_DIR / "36_screen_full.csv"
    cols = ["CHEMBL", "replica", "seed", "status", "qc_pass",
            "delta_H_PB_kcal_mol", "dG_gas", "dG_solv",
            "kept_frames", "total_frames", "kept_frac", "ctrl_tier",
            "zn_pass_rate", "zn_mean_min4_dist", "zn_mean_min4_dist_A",
            "prod_ns", "minutes", "notes"]
    with open(det, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # ---- 逐配体汇总 ----
    bylig: dict[str, list] = {}
    for r in rows:
        bylig.setdefault(r["CHEMBL"], []).append(r)
    agg = []
    for lig, rs in sorted(bylig.items()):
        ok = [r for r in rs if str(r.get("qc_pass")) == "1"]
        dh = [r["delta_H_PB_kcal_mol__n"] for r in ok
              if r["delta_H_PB_kcal_mol__n"] is not None]
        kf = [r["kept_frac"] for r in rs if r["kept_frac"] is not None]
        zn = [r["zn_pass_rate__n"] for r in rs if r["zn_pass_rate__n"] is not None]
        n_rep, n_ok = len(rs), len(ok)
        # 硬判据（管线预注册）
        verdict = ("STABLE" if n_ok == n_rep and n_rep >= 3
                   else "MARGINAL" if n_ok >= 1 else "UNSTABLE")
        # 对照校准档（经验）：以逐副本保留帧率的**中位数**定性，
        # 同时记录有多少副本各自越过阴性/阳性对照线（跨副本一致性）
        kfs = [r["kept_frac"] for r in rs if r["kept_frac"] is not None]
        if not kfs or neg_cut is None:
            ctrl_verdict = "NA"
        else:
            med = statistics.median(kfs)
            if med <= neg_cut:
                ctrl_verdict = "DMSO-like"
            elif pos_min is not None and med >= pos_min:
                ctrl_verdict = "LxxLL-like"
            else:
                ctrl_verdict = "intermediate"
        n_above_neg = sum(1 for x in kfs if neg_cut is not None and x > neg_cut)
        n_above_pos = sum(1 for x in kfs if pos_min is not None and x >= pos_min)
        agg.append({
            "CHEMBL": lig,
            "role": ("negative_control" if lig == NEG_CTRL
                     else "positive_control" if lig == POS_CTRL else "candidate"),
            "n_replicas": n_rep,
            "n_qc_pass": n_ok,
            "verdict_hard": verdict,
            "verdict_ctrl": ctrl_verdict,
            "n_rep_above_neg": n_above_neg,
            "n_rep_above_pos": n_above_pos,
            "kept_frac_median": round(statistics.median(kf), 4) if kf else "NA",
            "kept_frac_values": ";".join(f"{x:.3f}" for x in kf) if kf else "NA",
            "delta_H_mean": round(statistics.fmean(dh), 3) if dh else "NA",
            "delta_H_sd": round(statistics.stdev(dh), 3) if len(dh) > 1 else "NA",
            "delta_H_values": ";".join(f"{x:.2f}" for x in dh) if dh else "NA",
            "zn_pass_rate_min": round(min(zn), 3) if zn else "NA",
            "status_set": ";".join(sorted({r["status"] for r in rs})),
        })
    byl = OUT_DIR / "36_screen_byligand.csv"
    with open(byl, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(agg[0].keys()))
        w.writeheader()
        for a in agg:
            w.writerow(a)

    # ---- 可读报告 ----
    md = OUT_DIR / "36_screen_report.md"
    cands = [a for a in agg if a["role"] == "candidate"]
    n_st = sum(1 for a in cands if a["verdict_hard"] == "STABLE")
    n_mg = sum(1 for a in cands if a["verdict_hard"] == "MARGINAL")
    n_un = sum(1 for a in cands if a["verdict_hard"] == "UNSTABLE")
    n_ll = sum(1 for a in cands if a["verdict_ctrl"] == "LxxLL-like")
    n_im = sum(1 for a in cands if a["verdict_ctrl"] == "intermediate")
    n_dm = sum(1 for a in cands if a["verdict_ctrl"] == "DMSO-like")
    with open(md, "w") as fh:
        fh.write("# batch 2 全量筛选汇总（对照校准）\n\n")
        fh.write(f"- 明细 {len(rows)} 行（应为 {N_EXPECTED}），"
                 f"配体/对照 {len(agg)} 个\n")
        if neg_max is not None and pos_min is not None:
            fh.write(f"- 对照分离带：DMSO ≤{neg_max:.3f}（n={len(neg)}）"
                     f" < LxxLL ≥{pos_min:.3f}（n={len(pos)}）\n")
            fh.write(f"- 校准阈值：DMSO-like ≤ {neg_cut:.3f}；"
                     f"LxxLL-like ≥ {pos_min:.3f}\n")
        fh.write(f"- 硬判据：STABLE **{n_st}** / MARGINAL **{n_mg}** / "
                 f"UNSTABLE **{n_un}**（14 候选中）\n")
        fh.write(f"- 对照校准：LxxLL-like **{n_ll}** / intermediate **{n_im}** / "
                 f"DMSO-like **{n_dm}**\n\n")
        fh.write("| 对象 | 角色 | 副本QC | 硬判据 | 校准档 | 保留帧率中位 | "
                 "ΔH 均值±SD (kcal/mol) | Zn合格率(最低) |\n")
        fh.write("|---|---|---|---|---|---|---|---|\n")
        for a in agg:
            dh = (f"{a['delta_H_mean']}±{a['delta_H_sd']}"
                  if a["delta_H_mean"] != "NA" else "NA")
            fh.write(f"| {a['CHEMBL']} | {a['role']} | "
                     f"{a['n_qc_pass']}/{a['n_replicas']} | {a['verdict_hard']} | "
                     f"{a['verdict_ctrl']} | {a['kept_frac_median']} | {dh} | "
                     f"{a['zn_pass_rate_min']} |\n")
        fh.write("\n## 判据\n")
        fh.write("- 副本级 qc_pass=1：保留帧 ≥500（配体 RMSD<4 Å）"
                 "且双 Zn 四配位合格率 ≥0.95（管线预注册硬判据）\n")
        fh.write("- 对照校准档：由 DMSO（阴性）与 LxxLL 肽（阳性）实测"
                 "保留帧率界定，避免 500 帧切点的任意性\n")
        fh.write("- ΔH_PB 仅在 qc_pass=1 的副本上取均值；UNSTABLE 不给 ΔH\n")

    print(f"明细 -> {det}")
    print(f"逐配体 -> {byl}")
    print(f"报告 -> {md}")
    print(f"\n硬判据（候选 14）：STABLE {n_st} / MARGINAL {n_mg} / UNSTABLE {n_un}")
    print(f"对照校准（候选 14）：LxxLL-like {n_ll} / intermediate {n_im} / "
          f"DMSO-like {n_dm}")
    for a in agg:
        dh = (f"{a['delta_H_mean']}±{a['delta_H_sd']}"
              if a["delta_H_mean"] != "NA" else "NA")
        print(f"  {a['CHEMBL']:20s} {a['verdict_hard']:9s} "
              f"{a['verdict_ctrl']:12s} QC {a['n_qc_pass']}/{a['n_replicas']} "
              f"kf={a['kept_frac_median']}  ΔH={dh}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
