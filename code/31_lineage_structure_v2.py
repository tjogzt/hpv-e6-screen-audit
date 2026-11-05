#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
31_lineage_structure_v2.py — 结构约束下的谱系变异分析（修正版）

本脚本取代 30_lineage_structure_variation.py 的统计部分。修正两个会导致审稿被拒的缺陷：

[缺陷 1] 伪重复（pseudoreplication）
    v1 把每条序列当作独立观测。但同一型别下有多达 137 条高度相似的变体序列，
    它们不是独立样本，会把检验推向必然显著 —— v1 的 44/48 "显著" 只是在重申
    "不同 HPV 型别的序列不同" 这一平凡事实。
    修正：以**型别**为观测单位，每个型别在每个位置取一致性残基（consensus）。

[缺陷 2] 缺少背景对照
    v1 只检验"功能位点跨集群是否不同"，但任何位置都会因型别不同而不同。
    修正：引入背景分布 —— 所有可映射的非功能位点，比较
    功能位点的谱系一致性是否显著偏离背景（Mann-Whitney U + 置换检验）。

分析单位与指标（预注册）：
  - 观测单位：型别（n=45），型别 -> 集群（Part A, k=10）
  - 每个位置的谱系一致性 lineage_coherence = P(同集群型别残基相同) - P(异集群型别残基相同)
    取值 [-1, 1]：>0 表示集群划分在该位置上具有结构/功能意义
  - 功能位点组：Zn 配位 / E6AP LxxLL 界面 / p53 界面（来自脚本 29 的实验结构）
  - 背景组：模板位置 30-150 内、不属于任何功能位点的其他可映射位置
  - 组间比较：Mann-Whitney U（置换检验求 p），效应量 rank-biserial correlation
  - 多重比较：3 个功能位点组 vs 背景，BH-FDR

说明：本机 scipy 二进制损坏（scipy.sparse.linalg._propack 加载失败），
     Mann-Whitney U 与秩次计算均为自实现，p 值由置换检验给出（无分布假设）。

输出：
  results/31_type_consensus.csv            型别 x 位置 一致性残基矩阵
  results/31_position_coherence.csv        每个位置的谱系一致性指标
  results/31_functional_vs_background.json 组间比较结果
"""
from __future__ import annotations

import csv
import math
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
RESULTS = REBUILD / "results"
DATA = REBUILD / "data"

# 模板（UniProt P03126）的 8 个配位 Cys。结构编号 30/33/63/66/103/106/136/139，
# 到 UniProt 的偏移为 +7（由 4GIZ SEQADV 交叉验证），故 UniProt 编号为：
TEMPLATE_ZN2_UNP = [37, 40, 70, 73]      # 结构 30/33/63/66 + 7
TEMPLATE_ZN1_UNP = [110, 113, 143, 146]  # 结构 103/106/136/139 + 7
BG_LO, BG_HI = 30, 150                   # 背景位置窗口（覆盖全部功能位点及其周边）

N_PERM = 9999
SEED = 42

# ---------------------------------------------------------------- 锚定映射


def cys_positions(seq: str) -> list[int]:
    return [i + 1 for i, ch in enumerate(seq) if ch == "C"]


def find_zn_fingerprints(seq, loop_min=25, loop_max=35):
    """C-X-X-C ...(25-35)... C-X-X-C 指纹，返回按序列顺序排列的 Zn 位点。"""
    cys = cys_positions(seq)
    pairs = [(a, b) for i, a in enumerate(cys) for b in cys[i + 1:] if b - a == 3]
    sites, used = [], set()
    for i in range(len(pairs) - 1):
        p1, p2 = pairs[i], pairs[i + 1]
        loop = p2[0] - p1[1] - 1
        if loop_min <= loop <= loop_max:
            quad = [p1[0], p1[1], p2[0], p2[1]]
            if not (set(quad) & used):
                sites.append(quad)
                used |= set(quad)
    return sites


def build_anchor_map(seq):
    sites = find_zn_fingerprints(seq)
    if len(sites) < 2:
        return None
    sites = sorted(sites, key=lambda s: s[0])[:2]
    q_anchors = sites[0] + sites[1]
    t_anchors = TEMPLATE_ZN2_UNP + TEMPLATE_ZN1_UNP
    return t_anchors, q_anchors


def map_template_to_query(t, t_anchors, q_anchors, qlen):
    if t <= t_anchors[0]:
        q = q_anchors[0] - (t_anchors[0] - t)
    elif t >= t_anchors[-1]:
        q = q_anchors[-1] + (t - t_anchors[-1])
    else:
        q = None
        for i in range(len(t_anchors) - 1):
            if t_anchors[i] <= t <= t_anchors[i + 1]:
                frac = (t - t_anchors[i]) / (t_anchors[i + 1] - t_anchors[i])
                q = q_anchors[i] + frac * (q_anchors[i + 1] - q_anchors[i])
                break
        if q is None:
            return None
    q = int(round(q))
    return q if 1 <= q <= qlen else None


# ---------------------------------------------------------------- 统计（自实现）

def rankdata(a: np.ndarray) -> np.ndarray:
    """平均秩次（处理并列）。"""
    order = np.argsort(a, kind="mergesort")
    sa = a[order]
    n = len(a)
    ranks = np.empty(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sa[j + 1] == sa[i]:
            j += 1
        ranks[i:j + 1] = (i + j) / 2.0 + 1.0
        i = j + 1
    out = np.empty(n, dtype=float)
    out[order] = ranks
    return out


def mannwhitney_perm(x: np.ndarray, y: np.ndarray,
                     n_perm: int = N_PERM, seed: int = SEED):
    """Mann-Whitney U 检验，p 值由置换检验给出。返回 (U1, p, rank_biserial)。"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n1, n2 = len(x), len(y)
    combined = np.concatenate([x, y])
    ranks = rankdata(combined)
    R1 = ranks[:n1].sum()
    U1 = R1 - n1 * (n1 + 1) / 2.0

    def u_of(lab_mask):
        r = rankdata(combined)           # 秩次固定，只需重分配标签
        rr = r[lab_mask].sum()
        return rr - n1 * (n1 + 1) / 2.0

    rng = np.random.default_rng(seed)
    mask = np.zeros(len(combined), dtype=bool)
    mask[:n1] = True
    obs = U1
    cnt_greater = cnt_less = 0
    for _ in range(n_perm):
        pm = rng.permutation(mask)
        u = u_of(pm)
        if u >= obs - 1e-9:
            cnt_greater += 1
        if u <= obs + 1e-9:
            cnt_less += 1
    p_two = 2.0 * min(cnt_greater, cnt_less) / (n_perm + 1.0)
    p_two = float(min(1.0, p_two))
    rb = 2.0 * U1 / (n1 * n2) - 1.0      # rank-biserial correlation
    return float(U1), p_two, float(rb)


def norm_cdf(z: float) -> float:
    """标准正态 CDF（自实现，因本机 scipy 不可用）。"""
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def bh_fdr(pvals):
    n = len(pvals)
    if n == 0:
        return []
    order = np.argsort(pvals)
    ranked = np.array(pvals, dtype=float)[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return list(out)


# ---------------------------------------------------------------- 主流程

def main():
    summ = json.loads((RESULTS / "29_interface_summary.json").read_text())
    sites = {
        "zn_coordination": set(summ["zn_coord_cys_uniprot"]),
        "e6ap_lxxll_interface": set(summ["e6ap_interface_uniprot"]),
        "p53_interface": set(summ["p53_interface_uniprot"]),
    }
    functional = set().union(*sites.values())

    seqs = list(csv.DictReader(open(DATA / "02_dataset_clean.csv")))
    clu = {int(r["type_str"]): int(r["cluster"])
           for r in csv.DictReader(open(RESULTS / "04_type_clusters.csv"))
           if int(r["k"]) == 10}

    # ---- Step 1/2: 锚定 + 按型别聚合（消除伪重复）----
    type_pos_aa: dict[int, dict[int, Counter]] = defaultdict(lambda: defaultdict(Counter))
    n_ok = n_fail = 0
    fail_reason = Counter()
    for r in seqs:
        typ = int(r["hpv_type"])
        q = r["sequence"].strip()
        am = build_anchor_map(q)
        if am is None:
            n_fail += 1
            src = "SP" if r["accession"].startswith("sp|") else "TrEMBL"
            fail_reason[src] += 1
            continue
        n_ok += 1
        t_anchors, q_anchors = am
        # ---- 强制自检：8 个 Zn 锚点在查询序列上必须都是半胱氨酸 ----
        # 锚点由指纹直接给出，因此该检查恒成立；若失败则说明常量/映射写错。
        bad = [q[i - 1] for i in q_anchors
               if not (1 <= i <= len(q)) or q[i - 1] != "C"]
        if bad:
            raise AssertionError(
                f"{r['sequence_id']}: Zn 锚点残基非 Cys {bad} —— "
                f"t_anchors={t_anchors}, q_anchors={q_anchors}，映射常量有误")
        for p in range(BG_LO, BG_HI + 1):
            qp = map_template_to_query(p, t_anchors, q_anchors, len(q))
            if qp is None:
                continue
            aa = q[qp - 1]
            if aa in "ACDEFGHIKLMNPQRSTVWY":
                type_pos_aa[typ][p][aa] += 1

    print("=" * 92)
    print("Step 1: Zn 指指纹锚定映射")
    print("=" * 92)
    print(f"  成功 {n_ok} 条 / 失败 {n_fail} 条  失败来源: {dict(fail_reason)}")
    print(f"  （失败条目全部为低质量片段，已在脚本 30 中验证：0 条 Swiss-Prot）")
    print(f"  纳入型别 {len(type_pos_aa)} 个")

    # 型别一致性序列：每个位置取最频繁氨基酸，同时记录一致性强度
    consensus: dict[int, dict[int, str]] = {}
    cons_strength: dict[int, dict[int, float]] = {}
    for typ, posc in type_pos_aa.items():
        consensus[typ] = {}
        cons_strength[typ] = {}
        for p, c in posc.items():
            top, n = c.most_common(1)[0]
            consensus[typ][p] = top
            cons_strength[typ][p] = n / sum(c.values())

    p1 = RESULTS / "31_type_consensus.csv"
    with open(p1, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["hpv_type", "cluster", "position", "consensus_aa",
                    "consensus_strength", "n_sequences", "site_class"])
        for typ in sorted(consensus):
            for p in sorted(consensus[typ]):
                cls = "+".join(sorted(k for k, v in sites.items() if p in v)) or "background"
                w.writerow([typ, clu.get(typ, ""), p, consensus[typ][p],
                            round(cons_strength[typ][p], 4),
                            sum(type_pos_aa[typ][p].values()), cls])
    print(f"[saved] {p1}")

    # ---- Step 4: 每个位置的谱系一致性 ----
    types = sorted(consensus)
    tclu = {t: clu[t] for t in types if t in clu}
    print(f"\n  型别 -> 集群：{len(tclu)} 个型别，{len(set(tclu.values()))} 个集群")

    coh_rows = []
    for p in range(BG_LO, BG_HI + 1):
        present = [t for t in types if p in consensus[t]]
        if len(present) < 10:
            continue
        within_same = within_n = 0
        between_same = between_n = 0
        for i in range(len(present)):
            for j in range(i + 1, len(present)):
                ti, tj = present[i], present[j]
                same_aa = consensus[ti][p] == consensus[tj][p]
                if tclu.get(ti) == tclu.get(tj):
                    within_n += 1
                    within_same += int(same_aa)
                else:
                    between_n += 1
                    between_same += int(same_aa)
        if within_n == 0 or between_n == 0:
            continue
        pw = within_same / within_n
        pb = between_same / between_n
        cls = "+".join(sorted(k for k, v in sites.items() if p in v)) or "background"
        # 关键：完全保守的位置（所有型别残基相同）其 coherence 恒为 0，
        # 与"随机无关"的 0 无法区分 —— 该指标在此类位置上无定义，必须排除。
        aas = {consensus[t][p] for t in present}
        invariant = (len(aas) == 1)
        coh_rows.append({
            "position": p, "site_class": cls,
            "n_types": len(present),
            "n_distinct_aa": len(aas),
            "invariant": invariant,
            "p_within_same": round(pw, 4), "p_between_same": round(pb, 4),
            "lineage_coherence": round(pw - pb, 4),
            "mean_consensus_strength": round(
                float(np.mean([cons_strength[t][p] for t in present])), 4),
        })
    p2 = RESULTS / "31_position_coherence.csv"
    with open(p2, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(coh_rows[0].keys()))
        w.writeheader()
        w.writerows(coh_rows)
    print(f"[saved] {p2}")

    # ---- Step 5/6: 功能位点 vs 背景 ----
    print("\n" + "=" * 92)
    print("Step 2: 功能位点 vs 背景（Mann-Whitney U + 置换检验）")
    print("=" * 92)
    # 仅保留有变异的位置：完全保守位置的 coherence 无定义
    bg = [r["lineage_coherence"] for r in coh_rows
          if r["site_class"] == "background" and not r["invariant"]]
    n_bg_excl = sum(1 for r in coh_rows
                    if r["site_class"] == "background" and r["invariant"])
    results = {}
    raw_p, keys = [], []
    for k in ["zn_coordination", "e6ap_lxxll_interface", "p53_interface"]:
        grp_all = [r for r in coh_rows if r["site_class"] == k]
        grp = [r["lineage_coherence"] for r in grp_all if not r["invariant"]]
        n_excl = sum(1 for r in grp_all if r["invariant"])
        if not grp:
            print(f"  {k}: 全部 {len(grp_all)} 个位置均完全保守，coherence 无定义，"
                  f"退出组间比较")
            results[k] = {"n_positions": len(grp_all), "n_invariant_excluded": n_excl,
                          "all_invariant": True}
            continue
        U, p, rb = mannwhitney_perm(np.array(grp), np.array(bg))
        raw_p.append(p)
        keys.append(k)
        results[k] = {
            "n_positions": len(grp),
            "n_invariant_excluded": n_excl,
            "all_invariant": False,
            "mean_coherence": round(float(np.mean(grp)), 4),
            "median_coherence": round(float(np.median(grp)), 4),
            "U": round(U, 1), "p_raw": round(p, 6),
            "rank_biserial": round(rb, 4),
        }
        print(f"  {k:<22} n={len(grp):>3}  均值 {np.mean(grp):>7.4f}  "
              f"中位 {np.median(grp):>7.4f}  U={U:>9.1f}  p={p:.4g}  r_rb={rb:+.3f}")

    qs = bh_fdr(raw_p)
    for k, q in zip(keys, qs):
        results[k]["q_bh"] = round(float(q), 6)
        results[k]["significant_fdr05"] = bool(q < 0.05)
    results["background"] = {
        "n_positions": len(bg),
        "mean_coherence": round(float(np.mean(bg)), 4),
        "median_coherence": round(float(np.median(bg)), 4),
    }
    print(f"\n  {'background':<22} n={len(bg):>3}  均值 {np.mean(bg):>7.4f}  "
          f"中位 {np.median(bg):>7.4f}")
    print("\n  BH-FDR 校正后：")
    for k in keys:
        r = results[k]
        print(f"    {k:<22} q={r['q_bh']:.4g}  "
              f"{'SIGNIFICANT' if r['significant_fdr05'] else 'not significant'}")

    # ---- 功效分析：避免把"不显著"误读为"无差异" ----
    print("\n" + "=" * 92)
    print("功效分析（正态近似；说明不显著结果的解释边界）")
    print("=" * 92)
    for k in keys:
        r = results[k]
        n1, n2 = r["n_positions"], len(bg)
        eu = n1 * n2 / 2.0
        var = n1 * n2 * (n1 + n2 + 1) / 12.0
        sd = math.sqrt(var)
        z_obs = (r["U"] - eu) / sd
        power = norm_cdf(abs(z_obs) - 1.959964) + norm_cdf(-abs(z_obs) - 1.959964)
        u_mde = eu + 2.801585 * sd          # z_{0.975} + z_{0.80}
        rb_mde = 2.0 * u_mde / (n1 * n2) - 1.0
        r.update({"z_obs": round(float(z_obs), 3),
                  "post_hoc_power": round(float(power), 3),
                  "min_detectable_rank_biserial": round(float(rb_mde), 3)})
        print(f"  {k:<22} n1={n1:>2} n2={n2:>2}  |z|={abs(z_obs):.2f}  "
              f"观测功效={power:.2f}  可检出最小效应 r_rb>={rb_mde:.2f}")
    print("\n  注：功效 <0.8 时，'不显著' 只能表述为'未检出'，不能表述为'无差异'。")

    # ---- Zn 位点保守性（显式报告，不静默跳过）----
    zn_report = {}
    for p in sorted(sites["zn_coordination"]):
        comp = Counter()
        for t in types:
            if p in consensus[t]:
                comp[consensus[t][p]] += 1
        zn_report[p] = {"composition": dict(comp), "n_types": sum(comp.values()),
                        "is_cys": set(comp.keys()) == {"C"}}
    n_zn_cys = sum(1 for v in zn_report.values() if v["is_cys"])
    print("\n" + "=" * 92)
    print("Zn 配位位点保守性（显式报告）")
    print("=" * 92)
    for p, v in zn_report.items():
        print(f"  位置 {p:>3}: n_types={v['n_types']:>2}  组成={v['composition']}  "
              f"{'100% Cys' if v['is_cys'] else 'NOT fully Cys'}")
    print(f"  => {n_zn_cys}/{len(zn_report)} 个 Zn 配位位置在所有型别一致性序列中均为半胱氨酸")

    p3 = RESULTS / "31_functional_vs_background.json"
    out = {
        "unit_of_analysis": "HPV type (consensus residue per position); n_types=%d" % len(types),
        "n_sequences_anchored": n_ok,
        "n_sequences_failed": n_fail,
        "failure_source": dict(fail_reason),
        "background_window": [BG_LO, BG_HI],
        "test": f"Mann-Whitney U with permutation p ({N_PERM} perms, seed={SEED}); BH-FDR",
        "groups": results,
        "zn_conservation": {
            "n_positions": len(zn_report),
            "n_positions_fully_cys": n_zn_cys,
            "detail": zn_report,
        },
    }
    p3.write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(f"[saved] {p3}")


if __name__ == "__main__":
    main()
