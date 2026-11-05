#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
30_lineage_structure_variation.py — 结构约束下的谱系变异分析

目的：稿件 Part B 声称 lineage → pocket geometry / binding。直接比较同源模型的结合能
会被建模误差淹没（跨型同一性低至 26%）。本脚本改用一条零建模误差的路径：
把脚本 29 从实验结构算出的"结构关键位置集"锚定映射回全部 E6 序列，统计各谱系
（Part A 的 k=10 集群）在这些功能位点上的氨基酸组成差异。

映射方法（关键技术选择）：
  不用全局序列比对 —— 脚本 27/28 已证明 PairwiseAligner 在 <50% 同一性下会产生
  比对滑移，把真实存在的配位 Cys 错报为缺失。改用 **Zn 指指纹锚定**：
    - 每条序列用 C-X-X-C ...(25-35)... C-X-X-C 指纹定位 8 个配位 Cys（脚本 28 方法）
    - 这 8 个 Cys 是与模板结构一一对应的硬约束锚点
    - 锚点之间做分段线性插值，得到全序列的 序列位置 <-> 模板(UniProt)位置 映射
  这样得到的对应关系由结构约束保证，不依赖序列相似性。

三类结构关键位置（来自 results/29_interface_summary.json）：
  A. Zn 配位 Cys      （双锌指骨架）
  B. E6AP LxxLL 界面  （4GIZ 实测）
  C. p53 界面         （4XR8 实测）

统计（预注册）：
  - 每个位置 × 每个集群 的氨基酸组成与 Shannon 熵
  - 每个位置跨集群的差异：卡方检验（列联表 集群 × 氨基酸类别），BH-FDR 校正
  - 保守性：某位置在某集群内的优势氨基酸频率

输出：
  results/30_position_cluster_composition.csv
  results/30_position_conservation.csv
  results/30_lineage_variation_summary.json
"""
from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

# 说明：本机 scipy 二进制损坏（scipy.sparse.linalg._propack 加载失败），无法使用
# scipy.stats.chi2_contingency。改用置换检验：列联表普遍存在稀疏格（多个集群仅
# 1-3 条序列），卡方渐近近似本就不可靠；置换检验给出精确条件 p 值，且无分布假设。

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
RESULTS = REBUILD / "results"
DATA = REBUILD / "data"

# 模板（4GIZ chain C）的配位 Cys，结构编号；以及到 UniProt P03126 的偏移
TEMPLATE_ZN2_STRUCT = [30, 33, 63, 66]      # 序列靠前的 Zn 位点
TEMPLATE_ZN1_STRUCT = [103, 106, 136, 139]  # 序列靠后的 Zn 位点
GIZ_TO_UNP = 7
TEMPLATE_ZN2_UNP = [n + GIZ_TO_UNP for n in TEMPLATE_ZN2_STRUCT]
TEMPLATE_ZN1_UNP = [n + GIZ_TO_UNP for n in TEMPLATE_ZN1_STRUCT]

# 代表性理化分类，用于降维后的卡方检验（避免 20×10 稀疏列联表）
AA_CLASS = {
    "A": "small_hydrophobic", "V": "small_hydrophobic", "L": "small_hydrophobic",
    "I": "small_hydrophobic", "M": "small_hydrophobic",
    "F": "aromatic", "W": "aromatic", "Y": "aromatic",
    "S": "polar", "T": "polar", "N": "polar", "Q": "polar", "C": "polar",
    "K": "positive", "R": "positive", "H": "positive",
    "D": "negative", "E": "negative",
    "G": "glycine", "P": "proline",
}
CLASSES = ["small_hydrophobic", "aromatic", "polar", "positive", "negative",
           "glycine", "proline"]

# 预注册：置换次数与随机种子（每个随机步骤前立即固定，保证可复现）
N_PERM = 1999
SEED = 42


def cys_positions(seq: str) -> list[int]:
    return [i + 1 for i, ch in enumerate(seq) if ch == "C"]


def find_zn_fingerprints(seq: str, loop_min: int = 25, loop_max: int = 35) -> list[list[int]]:
    """返回按序列顺序排列的 Zn 位点，每个位点为 4 个配位 Cys 位置。"""
    cys = cys_positions(seq)
    pairs = []
    for i in range(len(cys)):
        for j in range(i + 1, len(cys)):
            if cys[j] - cys[i] == 3:
                pairs.append((cys[i], cys[j]))
    sites, used = [], set()
    for i in range(len(pairs) - 1):
        p1, p2 = pairs[i], pairs[i + 1]
        loop = p2[0] - p1[1] - 1
        if loop_min <= loop <= loop_max:
            quad = [p1[0], p1[1], p2[0], p2[1]]
            if not (set(quad) & used):      # 避免同一 Cys 被两个位点共用
                sites.append(quad)
                used |= set(quad)
    return sites


def build_anchor_map(query: str) -> tuple[list[int], list[int]] | None:
    """用 Zn 指纹建立 模板(UniProt)锚点 -> 查询序列位置 的映射。

    返回 (t_anchors, q_anchors)；若指纹不完整则返回 None。
    """
    sites = find_zn_fingerprints(query)
    if len(sites) < 2:
        return None
    # 按序列位置排序，取前两个完整位点：靠前的对应模板 Zn2，靠后的对应 Zn1
    sites = sorted(sites, key=lambda s: s[0])[:2]
    q_anchors = sites[0] + sites[1]
    t_anchors = TEMPLATE_ZN2_UNP + TEMPLATE_ZN1_UNP
    if len(q_anchors) != 8:
        return None
    return t_anchors, q_anchors


def map_template_to_query(t: int, t_anchors: list[int], q_anchors: list[int],
                          qlen: int) -> int | None:
    """把模板(UniProt)位置 t 映射到查询序列位置。锚点间线性插值，锚点外平移外推。"""
    q = None
    if t <= t_anchors[0]:
        q = q_anchors[0] - (t_anchors[0] - t)
    elif t >= t_anchors[-1]:
        q = q_anchors[-1] + (t - t_anchors[-1])
    else:
        for i in range(len(t_anchors) - 1):
            if t_anchors[i] <= t <= t_anchors[i + 1]:
                frac = (t - t_anchors[i]) / (t_anchors[i + 1] - t_anchors[i])
                q = q_anchors[i] + frac * (q_anchors[i + 1] - q_anchors[i])
                break
    if q is None:
        return None
    q = int(round(q))
    return q if 1 <= q <= qlen else None


def shannon_entropy(counts: Counter) -> float:
    tot = sum(counts.values())
    if tot == 0:
        return float("nan")
    return -sum((c / tot) * math.log2(c / tot) for c in counts.values() if c > 0)


def chi2_stat_fast(lab: np.ndarray, cls: np.ndarray,
                   n_l: int, n_c: int) -> float:
    """列联表卡方统计量（向量化）。lab/cls 为整数编码的观测标签。"""
    idx = lab * n_c + cls
    counts = np.bincount(idx, minlength=n_l * n_c).reshape(n_l, n_c).astype(float)
    row = counts.sum(axis=1, keepdims=True)
    col = counts.sum(axis=0, keepdims=True)
    tot = counts.sum()
    if tot == 0:
        return 0.0
    exp = row @ col / tot
    with np.errstate(divide="ignore", invalid="ignore"):
        safe = np.where(exp > 0, exp, 1.0)
        stat = np.where(exp > 0, (counts - exp) ** 2 / safe, 0.0)
    return float(stat.sum())


def permutation_test(lab: np.ndarray, cls: np.ndarray,
                     n_perm: int = 1999, seed: int = 42) -> tuple[float, float, int]:
    """置换检验：打乱集群标签，重算卡方，得到条件零分布。

    返回 (observed_chi2, p_value, dof)。固定随机种子以保证可复现。
    """
    n_l = int(lab.max()) + 1
    n_c = int(cls.max()) + 1
    obs = chi2_stat_fast(lab, cls, n_l, n_c)
    rng = np.random.default_rng(seed)   # 每个随机步骤前立即设种子
    cnt = 0
    for _ in range(n_perm):
        perm = rng.permutation(lab)
        if chi2_stat_fast(perm, cls, n_l, n_c) >= obs - 1e-12:
            cnt += 1
    p = (1 + cnt) / (1 + n_perm)
    dof = (n_l - 1) * (n_c - 1)
    return obs, p, dof


def bh_fdr(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg 校正。"""
    n = len(pvals)
    if n == 0:
        return []
    order = np.argsort(pvals)
    ranked = np.array(pvals)[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return list(out)


def main() -> None:
    # ---------- 载入结构关键位置 ----------
    summ = json.loads((RESULTS / "29_interface_summary.json").read_text())
    sites = {
        "zn_coordination": summ["zn_coord_cys_uniprot"],
        "e6ap_lxxll_interface": summ["e6ap_interface_uniprot"],
        "p53_interface": summ["p53_interface_uniprot"],
    }
    print("=" * 92)
    print("结构关键位置集（UniProt P03126 编号）")
    print("=" * 92)
    for k, v in sites.items():
        print(f"  {k:22s} n={len(v):2d}  {v}")

    # ---------- 载入序列与集群 ----------
    seqs = list(csv.DictReader(open(DATA / "02_dataset_clean.csv")))
    clu = {int(r["type_str"]): int(r["cluster"])
           for r in csv.DictReader(open(RESULTS / "04_type_clusters.csv"))
           if int(r["k"]) == 10}
    print(f"\n序列数 {len(seqs)}；型别->集群映射 {len(clu)} 个型别")

    # ---------- 对每条序列做锚定映射并提取关键位置残基 ----------
    # position -> cluster -> Counter(aa)
    pos_cluster_aa: dict[int, dict[int, Counter]] = defaultdict(lambda: defaultdict(Counter))
    all_positions = sorted({p for v in sites.values() for p in v})
    n_ok = n_fail = 0
    fail_by_type: Counter = Counter()

    for r in seqs:
        typ = int(r["hpv_type"])
        cl = clu.get(typ)
        if cl is None:
            continue
        q = r["sequence"].strip()
        am = build_anchor_map(q)
        if am is None:
            n_fail += 1
            fail_by_type[typ] += 1
            continue
        n_ok += 1
        t_anchors, q_anchors = am
        for p in all_positions:
            qp = map_template_to_query(p, t_anchors, q_anchors, len(q))
            if qp is None:
                continue
            pos_cluster_aa[p][cl][q[qp - 1]] += 1

    print(f"\n锚定映射成功 {n_ok} / 失败 {n_fail} 条")
    if fail_by_type:
        print(f"  失败型别分布（前 10）: {fail_by_type.most_common(10)}")

    # ---------- 输出 1：位置 × 集群 氨基酸组成 ----------
    comp_rows = []
    for p in all_positions:
        which = [k for k, v in sites.items() if p in v]
        for cl in sorted(pos_cluster_aa[p].keys()):
            c = pos_cluster_aa[p][cl]
            tot = sum(c.values())
            top_aa, top_n = c.most_common(1)[0]
            comp_rows.append({
                "uniprot_pos": p,
                "site_class": "+".join(sorted(which)),
                "cluster": cl,
                "n_seq": tot,
                "top_aa": top_aa,
                "top_freq": round(top_n / tot, 4),
                "entropy_bits": round(shannon_entropy(c), 4),
                "composition": ";".join(f"{a}:{n}" for a, n in c.most_common()),
            })
    p1 = RESULTS / "30_position_cluster_composition.csv"
    with open(p1, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(comp_rows[0].keys()))
        w.writeheader()
        w.writerows(comp_rows)
    print(f"[saved] {p1}")

    # ---------- 输出 2：位置层面的跨集群差异检验 ----------
    cons_rows = []
    raw_p = []
    for p in all_positions:
        which = [k for k, v in sites.items() if p in v]
        clusters = sorted(pos_cluster_aa[p].keys())
        if len(clusters) < 2:
            continue
        # 展开为逐条观测的 (集群, 理化类别) 标签，供置换检验使用
        labels, classes = [], []
        for cl in clusters:
            c = pos_cluster_aa[p][cl]
            for aa, n in c.items():
                klass = AA_CLASS.get(aa, "other")
                labels.extend([cl] * n)
                classes.extend([klass] * n)
        if len(labels) == 0:
            continue
        # 类别重新编码为连续整数，并去掉只含一个水平的退化情形
        uniq_cls = sorted(set(classes))
        if len(uniq_cls) < 2 or len(set(labels)) < 2:
            continue
        cls_idx = {c: i for i, c in enumerate(uniq_cls)}
        lab_idx = {c: i for i, c in enumerate(sorted(set(labels)))}
        lab_arr = np.array([lab_idx[x] for x in labels], dtype=int)
        cls_arr = np.array([cls_idx[x] for x in classes], dtype=int)
        chi2, pv, dof = permutation_test(lab_arr, cls_arr, n_perm=N_PERM, seed=SEED)
        allc = Counter()
        for cl in clusters:
            allc.update(pos_cluster_aa[p][cl])
        tot = sum(allc.values())
        top_aa, top_n = allc.most_common(1)[0]
        raw_p.append((p, pv))
        cons_rows.append({
            "uniprot_pos": p,
            "site_class": "+".join(sorted(which)),
            "n_clusters": len(clusters),
            "n_seq": int(tot),
            "global_top_aa": top_aa,
            "global_top_freq": round(top_n / tot, 4),
            "global_entropy_bits": round(shannon_entropy(allc), 4),
            "chi2": round(float(chi2), 3),
            "dof": int(dof),
            "p_raw": float(pv),
            "test": f"permutation({N_PERM}), seed={SEED}",
        })
    # BH-FDR
    qs = bh_fdr([r["p_raw"] for r in cons_rows])
    for r, q in zip(cons_rows, qs):
        r["q_bh"] = round(float(q), 6)
        r["significant_fdr05"] = bool(q < 0.05)

    p2 = RESULTS / "30_position_conservation.csv"
    with open(p2, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(cons_rows[0].keys()))
        w.writeheader()
        w.writerows(cons_rows)
    print(f"[saved] {p2}")

    # ---------- 汇总 ----------
    n_sig = sum(1 for r in cons_rows if r["significant_fdr05"])
    summary = {
        "n_sequences_anchored": n_ok,
        "n_sequences_failed": n_fail,
        "n_positions_tested": len(cons_rows),
        "n_significant_fdr05": n_sig,
        "mapping_method": "Zn-finger CXXC fingerprint anchored piecewise-linear mapping",
        "test": "chi-square contingency (cluster x physicochemical class), BH-FDR",
        "site_definitions": {k: v for k, v in sites.items()},
        "by_site_class": {},
    }
    for k in sites:
        sub = [r for r in cons_rows if k in r["site_class"]]
        if not sub:
            continue
        ent = [r["global_entropy_bits"] for r in sub]
        summary["by_site_class"][k] = {
            "n_positions": len(sub),
            "n_significant_fdr05": sum(1 for r in sub if r["significant_fdr05"]),
            "mean_entropy_bits": round(float(np.mean(ent)), 4),
            "max_entropy_bits": round(float(np.max(ent)), 4),
            "mean_top_freq": round(float(np.mean([r["global_top_freq"] for r in sub])), 4),
        }
    p3 = RESULTS / "30_lineage_variation_summary.json"
    p3.write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"[saved] {p3}")

    # ---------- 打印结论 ----------
    print("\n" + "=" * 92)
    print("结果汇总")
    print("=" * 92)
    print(f"  锚定成功 {n_ok} 条序列，检验 {len(cons_rows)} 个位置，"
          f"BH-FDR<0.05 显著 {n_sig} 个")
    print(f"\n  {'位点类别':<24} {'位置数':>6} {'显著':>5} {'平均熵(bits)':>13} "
          f"{'最大熵':>8} {'平均优势频率':>13}")
    for k, v in summary["by_site_class"].items():
        print(f"  {k:<24} {v['n_positions']:>6} {v['n_significant_fdr05']:>5} "
              f"{v['mean_entropy_bits']:>13.3f} {v['max_entropy_bits']:>8.3f} "
              f"{v['mean_top_freq']:>13.3f}")

    print("\n  跨集群差异最显著的前 12 个位置：")
    top = sorted(cons_rows, key=lambda r: (r["q_bh"], -r["global_entropy_bits"]))[:12]
    print(f"    {'位置':>5} {'类别':<26} {'熵':>6} {'优势AA':>7} {'频率':>7} "
          f"{'chi2':>8} {'q(BH)':>10} 显著")
    for r in top:
        print(f"    {r['uniprot_pos']:>5} {r['site_class']:<26} "
              f"{r['global_entropy_bits']:>6.2f} {r['global_top_aa']:>7} "
              f"{r['global_top_freq']:>7.3f} {r['chi2']:>8.1f} {r['q_bh']:>10.2e} "
              f"{'YES' if r['significant_fdr05'] else 'no'}")


if __name__ == "__main__":
    main()
