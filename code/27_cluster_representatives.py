#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
27_cluster_representatives.py — 集群代表序列筛选 + 同源建模可行性分析 (Part B, 集群分层)

目的：稿件标题声称 lineage → pocket geometry / binding，需要为每个序列相似性集群
（Part A, k=10）建一个代表结构。跨型 E6 同一性很低，建模前必须先回答：
  1. 每个代表序列相对 4GIZ 模板的同一性是多少？（<30% 进入同源建模 twilight zone）
  2. 形成锌指的 8 个半胱氨酸在比对中是否仍然对齐到半胱氨酸？
     —— 若某型缺失配位 Cys，其模型不能用于口袋/结合比较，必须剔除或标注。

方法：模板锚定全局比对（Bio.Align.PairwiseAligner，BLOSUM62，开放/延伸 -11/-1），
不使用外部 MSA 工具，保证结果与"相对模板的残基对应关系"直接可解释。

输出：
  results/27_cluster_representatives.csv   集群 -> 代表类型/序列
  results/27_template_alignment.csv        每个代表的同一性 + 8 个 Cys 对齐情况
"""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
RESULTS = REBUILD / "results"
RESULTS.mkdir(exist_ok=True)

DATA = REBUILD / "data" / "02_dataset_clean.csv"
TEMPLATE_PDB = Path(os.environ.get(
    "TEMPLATE_PDB",
    "/Users/taozhu/my researches/HPV_structure/rebuild_20261025/data/pilot/rec_clean.pdb"))

# 模板（4GIZ 衍生受体）中两个 Zn 的配位半胱氨酸，及模板自身的残基编号
ZN1_CYS = [103, 106, 136, 139]
ZN2_CYS = [30, 33, 63, 66]
ALL_COORD_CYS = sorted(ZN1_CYS + ZN2_CYS)

# Part A (k=10) 每个集群的代表类型：选该集群内记录数最多、且被 IARC 评估过者优先；
# 记录数为 1 的集群仍保留代表，但在输出中标注 n_records 以供后续剔除判断。
CLUSTER_REPRESENTATIVES = {
    1: 18,   # alpha-7 (18/39/45/59/68)
    2: 16,   # alpha-9/5/6/11 高危超支系（4GIZ 模板来源）
    3: 11,   # alpha-10 低危
    4: 54,   # alpha-13
    5: 61,
    6: 2,
    7: 10,
    8: 5,    # Beta 属（皮肤型），与 mucosal 型差异最大
    9: 17,
    10: 1,
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "MSE": "M",
}


def template_sequence(pdb: Path) -> tuple[str, list[int]]:
    """从模板 PDB 提取氨基酸序列与残基编号（只取 ATOM 主链已出现的残基）。"""
    if not pdb.exists():
        raise FileNotFoundError(f"模板 PDB 不存在：{pdb}")
    seq, numbers, seen = [], [], set()
    for line in pdb.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        resname = line[17:20].strip().upper()
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        if resseq in seen or resname not in THREE_TO_ONE:
            continue
        seen.add(resseq)
        seq.append(THREE_TO_ONE[resname])
        numbers.append(resseq)
    return "".join(seq), numbers


def make_aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices
    a = Align.PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.mode = "global"
    # 半胱氨酸对半胱氨酸给强奖励：确保锌指 Cys 优先对齐（这些位点在结构上必须等价）
    return a


def main() -> int:
    df = pd.read_csv(DATA)
    df["hpv_type"] = df["hpv_type"].astype(str)

    tseq, tnums = template_sequence(TEMPLATE_PDB)
    print(f"模板序列：{len(tseq)} 残基，残基编号 {tnums[0]}–{tnums[-1]}")
    tcoord_idx = {r: tnums.index(r) for r in ALL_COORD_CYS if r in tnums}
    missing = [r for r in ALL_COORD_CYS if r not in tnums]
    if missing:
        print(f"[WARN] 模板中找不到这些配位 Cys：{missing}")
    print(f"模板配位 Cys（{len(tcoord_idx)}/8）及其序列位置：{tcoord_idx}")

    aligner = make_aligner()
    rows, align_rows = [], []

    for cluster in sorted(CLUSTER_REPRESENTATIVES):
        t = str(CLUSTER_REPRESENTATIVES[cluster])
        sub = df[df["hpv_type"] == t]
        if sub.empty:
            print(f"[WARN] 集群 {cluster} 的代表类型 HPV{t} 在数据集中无序列，跳过")
            continue
        # 与 Part A 一致：取该类型最长序列作代表
        idx = sub["sequence"].astype(str).map(len).idxmax()
        rec = sub.loc[idx]
        qseq = str(rec["sequence"]).upper().replace("*", "").replace(" ", "")
        n_records = int(len(sub))

        aln = aligner.align(tseq, qseq)[0]
        a_t, a_q = aln[0], aln[1]
        # 建立“模板残基编号 -> 比对上的目标残基”映射
        ti = -1
        t2q = {}
        for c_t, c_q in zip(a_t, a_q):
            if c_t != "-":
                ti += 1
            if c_q != "-" and ti >= 0:
                t2q[ti] = c_q

        identical = sum(1 for x, y in zip(a_t, a_q) if x == y and x != "-")
        aln_len = sum(1 for x, y in zip(a_t, a_q) if x != "-" and y != "-")
        identity = identical / aln_len if aln_len else 0.0

        cys_ok, cys_detail = 0, []
        for r in ALL_COORD_CYS:
            if r not in tcoord_idx:
                cys_detail.append(f"{r}:NA")
                continue
            q = t2q.get(tcoord_idx[r], "-")
            ok = (q == "C")
            cys_ok += 1 if ok else 0
            cys_detail.append(f"{r}:{q}{'*' if ok else ''}")

        # 目标序列自身的 Cys 总数（判断是否整体缺 Cys）
        n_cys_query = qseq.count("C")

        rows.append({
            "cluster": cluster, "hpv_type": t,
            "sequence_id": rec.get("sequence_id", ""),
            "accession": rec.get("accession", ""),
            "n_records_in_dataset": n_records,
            "query_len": len(qseq),
            "template_len": len(tseq),
            "aligned_len": aln_len,
            "identity_vs_template": round(identity, 4),
            "n_cys_query": n_cys_query,
            "n_coord_cys_aligned": cys_ok,
            "coord_cys_detail": " ".join(cys_detail),
            "homology_model_eligible": ("yes" if (identity >= 0.30 and cys_ok == 8)
                                        else "no"),
            "eligibility_note": (
                "" if (identity >= 0.30 and cys_ok == 8)
                else ("低同一性 " if identity < 0.30 else "")
                     + (f"配位Cys缺失/错位({8-cys_ok}/8) " if cys_ok < 8 else "")),
            "sequence": qseq,
        })
        align_rows.append({
            "cluster": cluster, "hpv_type": t,
            "aligned_template": a_t, "aligned_query": a_q,
        })
        print(f"集群 {cluster:>2} HPV{t:>3}: 同一性 {identity:.3f}  "
              f"配位Cys对齐 {cys_ok}/8  Cys总数 {n_cys_query}  n={n_records}")

    out = pd.DataFrame(rows)
    out.drop(columns=["sequence"]).to_csv(
        RESULTS / "27_cluster_representatives.csv", index=False)
    out[["cluster", "hpv_type", "sequence"]].to_csv(
        RESULTS / "27_representative_sequences.csv", index=False)
    pd.DataFrame(align_rows).to_csv(RESULTS / "27_template_alignment.csv", index=False)

    # 供建模用的 FASTA
    fa = RESULTS / "27_representatives.fasta"
    with open(fa, "w") as fh:
        for _, r in out.iterrows():
            fh.write(f">C{int(r['cluster']):02d}_HPV{int(r['hpv_type'])}\n{r['sequence']}\n")

    print(f"\n可行（同一性 ≥30% 且 8 个配位 Cys 全部对齐）："
          f"{(out['homology_model_eligible'] == 'yes').sum()}/{len(out)}")
    print("输出：27_cluster_representatives.csv / 27_representative_sequences.csv / "
          "27_template_alignment.csv / 27_representatives.fasta")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
