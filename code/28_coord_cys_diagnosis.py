#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
28_coord_cys_diagnosis.py — 诊断"配位 Cys 缺失"是真实缺失还是比对滑移

背景：27_cluster_representatives.py 报告 HPV18 (cluster 1) 与 HPV10 (cluster 7)
各缺失 1 个配位 Cys（对齐位置 139 分别是 E 和 W）。但 HPV18 是公认高危型，
其 E6 必须具备完整双锌指，否则无法结合 E6AP/p53 —— 生物学上不成立。
因此需要先判断：这是真实的残基替换，还是低同一性下的比对滑移。

E6 锌指的结构指纹（Hofmann 1996 / Nomine 2006 / Zanier 2013）：
    C-X-X-C ...(约 29 残基)... C-X-X-C        → 一个 Zn 位点（4 配位）
    全长为两个这样的单元，共 4 对 = 8 个配位 Cys
    HPV16 E6 (P03126) 配位 Cys：30/33/63/66/103/106/136/139

本脚本不依赖任何 MSA 工具，直接对每条代表序列做：
  1. 列出全部 Cys 位置
  2. 用 C-X-X-C 指纹枚举候选配位 Cys（pair 内间隔恰为 2，pair 间间隔 25-35）
  3. 与模板配位 Cys 的模板锚定比对结果对照，判定"缺失"的性质
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
RESULTS = REBUILD / "results"

# 模板（4GIZ 衍生受体，残基编号 2-143）中两个 Zn 的配位 Cys
ZN1_CYS = [103, 106, 136, 139]
ZN2_CYS = [30, 33, 63, 66]
ALL_COORD_CYS = sorted(ZN1_CYS + ZN2_CYS)

# 代表序列相对模板的残基编号偏移（受体 = 模板序列 2..143，故模板编号 = 序列编号 + 2）
# 该偏移由 27 脚本的模板锚定比对内部确定，这里用于把比对位置换算回序列位置
SEQ_OFFSET = 2


def cys_positions(seq: str) -> list[int]:
    return [i + 1 for i, ch in enumerate(seq) if ch == "C"]


def find_zn_fingerprints(seq: str, loop_min: int = 25, loop_max: int = 35) -> list[dict]:
    """按 C-X-X-C ...(25-35)... C-X-X-C 指纹枚举候选 Zn 位点。

    返回每个候选位点的 4 个 Cys 位置与 loop 长度。
    """
    cys = cys_positions(seq)
    # 找所有 C-X-X-C pair（位置差恰为 3 => 中间隔 2 个残基）
    # 注意：不能只检查相邻 Cys —— 若中间夹了一个非配位 Cys（如 C-X-C-X-X-C），
    # 相邻检查会漏掉真正的配位对（HPV11 的 C64/C67 中间夹着 C66 即属此类）
    pairs = []
    for i in range(len(cys)):
        for j in range(i + 1, len(cys)):
            if cys[j] - cys[i] == 3:
                pairs.append((cys[i], cys[j]))
    sites = []
    for i in range(len(pairs) - 1):
        p1, p2 = pairs[i], pairs[i + 1]
        loop = p2[0] - p1[1] - 1  # 两个 pair 之间的残基数
        if loop_min <= loop <= loop_max:
            sites.append({
                "cys": [p1[0], p1[1], p2[0], p2[1]],
                "loop": loop,
                "pattern": f"C{p1[0]}-XX-C{p1[1]}-({loop})-C{p2[0]}-XX-C{p2[1]}",
            })
    return sites


def main() -> None:
    seqs = list(csv.DictReader(open(RESULTS / "27_representative_sequences.csv")))
    # identity / 比对得到的配位 Cys 对齐情况来自 27 的可行性输出
    aln = {int(r["cluster"]): r for r in
           csv.DictReader(open(RESULTS / "27_cluster_representatives.csv"))}

    out_rows = []
    print("=" * 100)
    print("代表序列 Zn 指指纹诊断（C-X-X-C ...(25-35)... C-X-X-C）")
    print("=" * 100)

    for r in seqs:
        cl = int(r["cluster"])
        typ = r["hpv_type"]
        seq = r["sequence"].strip()
        cys = cys_positions(seq)
        sites = find_zn_fingerprints(seq)

        # 模板锚定比对给出的配位 Cys 对齐情况
        a = aln.get(cl, {})
        detail = a.get("coord_cys_detail", "")
        n_aligned = a.get("n_coord_cys_aligned", "")
        identity = a.get("identity_vs_template", "")

        # 收集指纹命中的所有 Cys
        fingerprint_cys = sorted({c for s in sites for c in s["cys"]})

        print(f"\n--- cluster {cl} / HPV{typ}  (len={len(seq)}, n_cys={len(cys)}, "
              f"identity={identity}) ---")
        print(f"  Cys 位置: {cys}")
        if sites:
            for s in sites:
                print(f"  指纹位点: {s['pattern']}")
            print(f"  指纹覆盖的 Cys: {fingerprint_cys}  (共 {len(fingerprint_cys)} 个)")
        else:
            print("  指纹位点: 未检出")
        print(f"  模板锚定比对: {n_aligned}/8  {detail}")

        # 判定：指纹是否给出完整的 8 个配位 Cys
        if len(fingerprint_cys) >= 8:
            verdict = "COMPLETE_FINGERPRINT"
            note = (f"序列自身含 {len(fingerprint_cys)} 个符合 Zn 指指纹的 Cys"
                    f"（>=8），比对报缺属滑移")
        elif len(fingerprint_cys) >= 4:
            verdict = "PARTIAL_FINGERPRINT"
            note = (f"仅检出 {len(fingerprint_cys)} 个指纹 Cys，可能仅一个 Zn 位点完整"
                    f"，或 loop 长度超出 25-35 窗口")
        else:
            verdict = "NO_FINGERPRINT"
            note = "未检出成对 Zn 指指纹"

        # 未进入指纹的"孤立 Cys"
        orphan = [c for c in cys if c not in fingerprint_cys]
        print(f"  非指纹 Cys: {orphan}")
        print(f"  判定: {verdict} — {note}")

        out_rows.append({
            "cluster": cl,
            "hpv_type": typ,
            "seq_len": len(seq),
            "n_cys": len(cys),
            "cys_positions": " ".join(map(str, cys)),
            "identity_vs_template": identity,
            "n_coord_cys_aligned": n_aligned,
            "alignment_detail": detail,
            "n_fingerprint_sites": len(sites),
            "fingerprint_patterns": " | ".join(s["pattern"] for s in sites),
            "fingerprint_cys": " ".join(map(str, fingerprint_cys)),
            "n_fingerprint_cys": len(fingerprint_cys),
            "orphan_cys": " ".join(map(str, orphan)),
            "verdict": verdict,
            "note": note,
        })

    outp = RESULTS / "28_coord_cys_diagnosis.csv"
    with open(outp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"\n[saved] {outp}")

    # 汇总
    print("\n" + "=" * 100)
    print("汇总")
    print("=" * 100)
    print(f"{'cluster':>7} {'HPV':>5} {'ident':>6} {'nCys':>5} {'指纹Cys':>7} {'比对':>5}  {'判定'}")
    for r in out_rows:
        print(f"{r['cluster']:>7} {r['hpv_type']:>5} "
              f"{float(r['identity_vs_template']):>6.3f} {r['n_cys']:>5} "
              f"{r['n_fingerprint_cys']:>7} {str(r['n_coord_cys_aligned'])+'/8':>6}  "
              f"{r['verdict']}")


if __name__ == "__main__":
    main()
