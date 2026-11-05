#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
29_define_interface_positions.py — 从实验结构定义 E6 的"结构关键位置集"

目的：稿件 Part B 声称 lineage → pocket geometry / binding。要检验这个 claim，
必须先把 4GIZ(HPV16 E6) 上"哪些残基构成功能位点"从实验结构中算出来，而不是引用
记忆中的残基列表。本脚本纯用坐标计算，产出可核验的位置集。

三类位置：
  A. Zn 配位 Cys（双锌指骨架）—— 由脚本 28 的 C-X-X-C 指纹独立给出，此处用结构坐标复核
  B. E6AP LxxLL 结合界面 —— 4GIZ 中 E6(chain C) 与 MBP-E6AP 融合蛋白 C 端 LxxLL motif 的接触残基
  C. p53 结合界面      —— 4XR8 中 E6(chain F/H) 与 p53 DBD(chain C/D) 的接触残基

同时核查工程化突变（4GIZ: F47R + C80S/C97S/C111S/C140S）是否落在界面内 ——
若落在界面内，所有基于该模板的对接/打分都必须显式声明这一 caveat。

判定规则（预注册）：任一残基的重原子与配体/伙伴链任一重原子距离 < CUTOFF (5.0 Å)
即定义为界面残基。Zn 配位判定：Cys SG 与 ZN 距离 < 3.0 Å。

输出：
  results/29_interface_positions.csv   每个位置的归类（可多选）
  results/29_interface_summary.json    汇总 + 突变核查
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
RESULTS = REBUILD / "results"
RESULTS.mkdir(exist_ok=True)

STRUCT = Path("/Users/taozhu/my researches/HPV_structure/protein_structures")
GIZ = STRUCT / "4GIZ.pdb"
XR8 = STRUCT / "4XR8.pdb"

CUTOFF = 5.0        # 界面判定阈值 (Å)
ZN_CUTOFF = 3.0     # Zn 配位判定阈值 (Å)

# 4GIZ 中 HPV16 E6 的工程化突变（结构编号 -> 野生型），来自 SEQADV 记录
GIZ_ENGINEERED = {47: "F47R", 80: "C80S", 97: "C97S", 111: "C111S", 140: "C140S"}
# 结构编号 -> UniProt(P03126) 编号 的偏移（由 SEQADV 交叉验证：结构 47↔UNP 54 等，偏移恒为 7）
GIZ_TO_UNP_OFFSET = 7

# E6AP LxxLL motif 核心序列（用于定位融合蛋白中的肽段）
LXXLL_MOTIF = "LTLQELL"


def heavy_atoms(res):
    return np.array([a.coord for a in res if a.element != "H"], dtype=float)


def min_dist(a: np.ndarray, b: np.ndarray) -> float:
    """两组坐标之间的最小距离（Å）。"""
    if len(a) == 0 or len(b) == 0:
        return float("inf")
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)
    return float(d.min())


def seq_of(chain) -> tuple[str, list[int]]:
    res = [r for r in chain if r.id[0] == " "]
    s = "".join(protein_letters_3to1.get(r.get_resname(), "X") for r in res)
    return s, [r.id[1] for r in res]


def locate_motif(seq: str, nums: list[int], motif: str) -> list[int]:
    """返回 motif 在链中的残基编号列表。"""
    idx = seq.find(motif)
    if idx < 0:
        return []
    return nums[idx: idx + len(motif)]


def main() -> None:
    parser = PDBParser(QUIET=True)
    out_rows = []

    # ---------------- A. Zn 配位（4GIZ chain C，含 2 个 ZN） ----------------
    print("=" * 92)
    print("A. Zn 配位复核（4GIZ chain C vs ZN 201/202）")
    print("=" * 92)
    giz = parser.get_structure("4giz", str(GIZ))[0]
    e6 = giz["C"]
    zn_coords = {}
    for res in giz["C"]:
        if res.get_resname() == "ZN":
            zn_coords[res.id[1]] = np.array(
                [a.coord for a in res if a.element != "H"], dtype=float)
    print(f"  chain C 中 Zn 数: {len(zn_coords)} 编号={sorted(zn_coords)}")

    zn_cys = {}
    for zn_id, zc in sorted(zn_coords.items()):
        hits = []
        for res in e6:
            if res.id[0] != " ":
                continue
            sg = np.array([a.coord for a in res if a.get_name() == "SG"], dtype=float)
            if len(sg) == 0:
                continue
            d = min_dist(sg, zc)
            if d < ZN_CUTOFF:
                hits.append((res.id[1], round(d, 3)))
        hits.sort()
        zn_cys[zn_id] = hits
        print(f"  ZN {zn_id}: 配位 Cys = {[h[0] for h in hits]}  距离 = {[h[1] for h in hits]}")

    all_zn_cys = sorted({c for hits in zn_cys.values() for c, _ in hits})
    print(f"  => 结构坐标确认的配位 Cys ({len(all_zn_cys)}): {all_zn_cys}")

    # ---------------- B. E6AP LxxLL 界面（4GIZ chain C vs chain A 的 LxxLL） ----------------
    print("\n" + "=" * 92)
    print("B. E6AP LxxLL 结合界面（4GIZ chain C vs chain A LxxLL motif）")
    print("=" * 92)
    chA = giz["A"]
    seqA, numsA = seq_of(chA)
    motif_nums = locate_motif(seqA, numsA, LXXLL_MOTIF)
    print(f"  chain A 长度 {len(seqA)}，LxxLL motif '{LXXLL_MOTIF}' 位于残基 {motif_nums}")

    # 肽段取 motif 及其两侧各 3 个残基，覆盖完整结合肽
    lo, hi = min(motif_nums) - 3, max(motif_nums) + 3
    pep_res = [r for r in chA if r.id[0] == " " and lo <= r.id[1] <= hi]
    pep_xyz = np.vstack([heavy_atoms(r) for r in pep_res if len(heavy_atoms(r))])
    print(f"  结合肽段取 {lo}-{hi}（含 motif 两侧各 3 残基），共 {len(pep_res)} 残基")

    e6ap_iface = []
    for res in e6:
        if res.id[0] != " ":
            continue
        rc = heavy_atoms(res)
        if len(rc) == 0:
            continue
        d = min_dist(rc, pep_xyz)
        if d < CUTOFF:
            e6ap_iface.append((res.id[1], round(d, 2), res.get_resname()))
    e6ap_iface.sort()
    print(f"  => E6AP 界面残基 ({len(e6ap_iface)}): "
          f"{[f'{n}{protein_letters_3to1.get(rn, chr(88))}' for n, d, rn in e6ap_iface]}")

    # ---------------- C. p53 界面（4XR8 E6 vs p53 DBD） ----------------
    print("\n" + "=" * 92)
    print("C. p53 结合界面（4XR8 E6 chain F/H vs p53 chain C/D）")
    print("=" * 92)
    xr8 = parser.get_structure("4xr8", str(XR8))[0]
    # p53 DBD 特征序列片段，用于确认哪条链是 p53
    p53_marker = "MFCQLAKTCPVQLWVDSTPPP"
    e6_chains, p53_chains = [], []
    for chain in xr8:                      # Model 迭代返回 Chain 对象，需取 .id
        cid = chain.id
        s, _ = seq_of(chain)
        if p53_marker in s:
            p53_chains.append(cid)
        elif "FQDPQERPRKLPQLCTEL" in s:
            e6_chains.append(cid)
    print(f"  p53 链: {p53_chains}   E6 链: {e6_chains}")

    # 所有 p53 重原子合并，用于突变位点的距离核查
    p53_all_xyz = None
    _stack = []
    for pid in p53_chains:
        for r in xr8[pid]:
            if r.id[0] == " " and len(heavy_atoms(r)):
                _stack.append(heavy_atoms(r))
    if _stack:
        p53_all_xyz = np.vstack(_stack)

    p53_iface = {}
    for cid in e6_chains:
        e6c = xr8[cid]
        best = []
        for pid in p53_chains:
            pc = xr8[pid]
            p_xyz = np.vstack([heavy_atoms(r) for r in pc
                               if r.id[0] == " " and len(heavy_atoms(r))])
            hits = []
            for res in e6c:
                if res.id[0] != " ":
                    continue
                rc = heavy_atoms(res)
                if len(rc) == 0:
                    continue
                d = min_dist(rc, p_xyz)
                if d < CUTOFF:
                    hits.append((res.id[1], round(d, 2), res.get_resname()))
            hits.sort()
            print(f"  E6 chain {cid} vs p53 chain {pid}: {len(hits)} 个界面残基")
            if len(hits) > len(best):
                best = hits
        p53_iface[cid] = best
        print(f"  => E6 chain {cid} p53 界面: "
              f"{[f'{n}{protein_letters_3to1.get(rn, chr(88))}' for n, d, rn in best]}")

    # 4XR8 与 4GIZ 的 E6 编号体系一致（偏移 0），已由 SEQADV 交叉验证：
    # 4XR8 SER F80/H80 <-> UNP CYS 87，与 4GIZ SER C80 同一位置。
    # 注意 4XR8 只带 4 个 C->S（80/97/111/140），不含 4GIZ 的 F47R。
    XR8_TO_GIZ = 0
    p53_iface_giz = sorted({n + XR8_TO_GIZ for hits in p53_iface.values()
                            for n, d, rn in hits})
    print(f"\n  => p53 界面换算到 4GIZ 编号体系 ({len(p53_iface_giz)}): {p53_iface_giz}")

    # 构建 4XR8 E6 的 编号->坐标 表（两个拷贝 F/H 取距离更小者）
    xr8_e6_coord = {}
    for cid in e6_chains:
        for res in xr8[cid]:
            if res.id[0] != " ":
                continue
            rc = heavy_atoms(res)
            if len(rc) == 0:
                continue
            num = res.id[1] + XR8_TO_GIZ
            if p53_all_xyz is None:
                continue
            d = min_dist(rc, p53_all_xyz)
            if num not in xr8_e6_coord or d < xr8_e6_coord[num][1]:
                xr8_e6_coord[num] = (rc, d)
    xr8_e6_coord = {k: v[0] for k, v in xr8_e6_coord.items()}
    print(f"  => 已构建 4XR8 E6 坐标表（{len(xr8_e6_coord)} 残基）用于突变距离核查")

    # ---------------- D. 工程化突变是否落在界面内 ----------------
    print("\n" + "=" * 92)
    print("D. 工程化突变核查（4GIZ: F47R + C80S/C97S/C111S/C140S）")
    print("=" * 92)
    e6ap_set = {n for n, d, rn in e6ap_iface}
    p53_set = set(p53_iface_giz)
    zn_set = set(all_zn_cys)
    mut_rows = []
    for num, label in sorted(GIZ_ENGINEERED.items()):
        in_e6ap = num in e6ap_set
        in_p53 = num in p53_set
        in_zn = num in zn_set
        if in_e6ap or in_p53:
            risk = "HIGH — 位于功能界面，必须声明"
        elif in_zn:
            risk = "CRITICAL — 位于 Zn 配位位点"
        else:
            risk = "LOW — 不在已定义的功能位点内"
        # 实测最小距离，避免 5 Å 二值判定掩盖临界情形
        d_e6ap = d_p53 = float("nan")
        for res in e6:
            if res.id[0] == " " and res.id[1] == num:
                rc = heavy_atoms(res)
                d_e6ap = round(float(min_dist(rc, pep_xyz)), 2)
                break
        # p53 距离必须用 4XR8 的 E6 坐标（4GIZ 与 4XR8 是不同的晶体坐标系）
        if p53_all_xyz is not None and xr8_e6_coord:
            rc = xr8_e6_coord.get(num)
            if rc is not None:
                d_p53 = round(float(min_dist(rc, p53_all_xyz)), 2)
        print(f"  {label:>7} (结构 {num:>3}, UniProt {num + GIZ_TO_UNP_OFFSET:>3}): "
              f"E6AP界面={in_e6ap}(d={d_e6ap})  p53界面={in_p53}(d={d_p53})  "
              f"Zn配位={in_zn}  → {risk}")
        mut_rows.append({
            "mutation": label, "struct_num": num,
            "uniprot_num": num + GIZ_TO_UNP_OFFSET,
            "in_e6ap_interface": in_e6ap, "in_p53_interface": in_p53,
            "in_zn_site": in_zn, "risk": risk,
            "min_dist_to_lxxll_A": d_e6ap, "min_dist_to_p53_A": d_p53,
        })

    # ---------------- 输出位置集 ----------------
    all_pos = sorted(zn_set | e6ap_set | p53_set)
    for n in all_pos:
        rn = None
        for res in e6:
            if res.id[0] == " " and res.id[1] == n:
                rn = res.get_resname()
                break
        out_rows.append({
            "struct_num": n,
            "uniprot_num": n + GIZ_TO_UNP_OFFSET,
            "resname_4giz": protein_letters_3to1.get(rn, "X") if rn else "",
            "zn_coordination": n in zn_set,
            "e6ap_lxxll_interface": n in e6ap_set,
            "p53_interface": n in p53_set,
            "engineered_mutation": GIZ_ENGINEERED.get(n, ""),
            "n_functional_sites": int(n in zn_set) + int(n in e6ap_set) + int(n in p53_set),
        })

    outp = RESULTS / "29_interface_positions.csv"
    import csv
    with open(outp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"\n[saved] {outp}")

    summary = {
        "cutoff_interface_A": CUTOFF,
        "cutoff_zn_A": ZN_CUTOFF,
        "zn_coord_cys_struct": all_zn_cys,
        "zn_coord_cys_uniprot": [n + GIZ_TO_UNP_OFFSET for n in all_zn_cys],
        "e6ap_interface_struct": sorted(e6ap_set),
        "e6ap_interface_uniprot": sorted(n + GIZ_TO_UNP_OFFSET for n in e6ap_set),
        "p53_interface_struct": p53_iface_giz,
        "p53_interface_uniprot": sorted(n + GIZ_TO_UNP_OFFSET for n in p53_iface_giz),
        "n_zn": len(zn_set), "n_e6ap": len(e6ap_set), "n_p53": len(p53_iface_giz),
        "n_total_positions": len(all_pos),
        "engineered_mutations": mut_rows,
        "giz_to_uniprot_offset": GIZ_TO_UNP_OFFSET,
        "structure_notes": (
            "4GIZ = HPV16 E6 (chain C/D, 142 res, numbering 2-143) + MBP-E6AP fusion "
            "(chain A/B); E6AP LxxLL motif 'LTLQELL' sits at the C-terminus of the "
            "fusion chain. 4XR8 = HPV16 E6 (chain F/H) + p53 DBD (chain C/D) + "
            "MBP-E6AP fusion (chain A/B). 4GIZ carries 5 engineered mutations "
            "(F47R, C80S, C97S, C111S, C140S) relative to UniProt P03126."
        ),
    }
    sp = RESULTS / "29_interface_summary.json"
    sp.write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"[saved] {sp}")

    # 交叉验证：Zn 配位 Cys 是否与脚本 28 的序列指纹一致
    fp = {37, 40, 70, 73, 110, 113, 143, 146}  # 28 脚本指纹给出的 UniProt 位置
    got = set(summary["zn_coord_cys_uniprot"])
    print("\n" + "=" * 92)
    print("交叉验证：结构坐标 vs 序列指纹（脚本 28）")
    print("=" * 92)
    print(f"  结构坐标给出 (UniProt): {sorted(got)}")
    print(f"  序列指纹给出 (UniProt): {sorted(fp)}")
    print(f"  一致: {got == fp}" + ("" if got == fp else f"  差异: {got ^ fp}"))


if __name__ == "__main__":
    main()
