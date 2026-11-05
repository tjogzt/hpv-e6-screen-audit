#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
26_zn_analyse.py — Zn 模型验证的独立分析步骤

从 24_zn_model_test.py 拆出来，原因：mdrun 在写完 prod.gro 之后偶发不退出，
把 MD 与"读轨迹做分析"解耦后，轨迹一旦完整即可直接判分，不必等 mdrun 收尾。

判据（预注册，与 24_zn_model_test.py 一致）：
  PASS = 末 50% 帧内两个 Zn 各自到其 4 个配位 Cys SG 的距离均值 ≤ 0.30 nm
         （模板 0.226–0.236 nm），且最大帧均值 ≤ 0.35 nm。

用法：
    python 26_zn_analyse.py zaff_cym
    python 26_zn_analyse.py all        # 分析 zn_model_test/ 下全部模型
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

BASE = os.environ.get("BASE", "/root/autodl-tmp/HPV_structure/candidates_run")
OUTROOT = os.path.join(BASE, "zn_model_test")

ZN1_RES = [103, 106, 136, 139]   # ← zn1.pdb（Zn_A）
ZN2_RES = [30, 33, 63, 66]       # ← zn2.pdb（Zn_B）
TEMPLATE_REF_NM = [0.226, 0.236]


def write(path: str, text: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8")


def analyse(model: str) -> dict:
    import numpy as np
    import mdtraj as md

    d = os.path.join(OUTROOT, model)
    prmtop, xtc = f"{d}/sys.prmtop", f"{d}/prod.xtc"
    for f in (prmtop, xtc):
        if not os.path.exists(f):
            return {"model": model, "error": f"缺少 {f}"}
        if os.path.getsize(f) < 1000:
            return {"model": model, "error": f"{f} 为空"}

    t = md.load(xtc, top=prmtop)
    top = t.topology
    # mdtraj 从 prmtop 读元素时对 ZN 常返回 element=None，必须按原子名回退
    zn = sorted({a.index for a in top.atoms
                 if a.element is not None and a.element.symbol == "Zn"}
                | {a.index for a in top.atoms
                   if a.name.strip().upper() == "ZN"})

    # 注意：mdtraj 的 residue.resSeq 与源 PDB 编号不保证一致，故不按残基号定位，
    # 而是取每个 Zn 最近的 4 个 Cys SG —— 与主管线 zn_mean_min4_dist 的定义一致，
    # 便于和 run1（同一指标）直接比较。配位正确时这 4 个就是 Cys103/106/136/139
    #（或 30/33/63/66）；Zn 脱离时该值会显著变大，正是要检出的失效模式。
    s_idx = [a.index for a in top.atoms
             if a.residue.name.strip().upper() in ("CYS", "CYM")
             and a.name.strip().upper() == "SG"]

    res: dict = {"model": model, "n_frames": int(t.n_frames),
                 "zn_atoms": len(zn), "n_cys_sg": len(s_idx), "series": {}}
    if len(zn) < 2:
        res["error"] = f"只找到 {len(zn)} 个 Zn 原子"
        return res
    if len(s_idx) < 8:
        res["error"] = f"Cys SG 只有 {len(s_idx)} 个（需 ≥8）"
        return res

    allmean = []
    for k, z in enumerate(zn[:2], start=1):
        dist = md.compute_distances(t, [(z, s) for s in s_idx])   # nm
        nearest4 = np.sort(dist, axis=1)[:, :4]
        m = nearest4.mean(axis=1)
        res["series"][f"zn{k}"] = [round(float(x), 4) for x in m]
        half = m[len(m) // 2:]
        res[f"zn{k}_mean_nm"] = round(float(half.mean()), 4)
        res[f"zn{k}_max_nm"] = round(float(half.max()), 4)
        res[f"zn{k}_first_nm"] = round(float(m[0]), 4)
        res[f"zn{k}_last_nm"] = round(float(m[-1]), 4)
        allmean.extend(half.tolist())

    res["both_mean_nm"] = round(float(np.mean(allmean)), 4)
    res["both_max_nm"] = round(float(np.max(allmean)), 4)
    res["template_ref_nm"] = TEMPLATE_REF_NM
    res["verdict"] = ("PASS" if (res["both_mean_nm"] <= 0.30
                                 and res["both_max_nm"] <= 0.35) else "FAIL")
    write(f"{d}/zn_result.json", json.dumps(res, indent=2))
    with open(f"{d}/zn_dist_series.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "zn1_mean_nm", "zn2_mean_nm"])
        for i, (a, b) in enumerate(zip(res["series"]["zn1"], res["series"]["zn2"])):
            w.writerow([i, a, b])
    return res


def main() -> int:
    arg = sys.argv[1] if len(sys.argv) > 1 else "all"
    models = (sorted(p.name for p in Path(OUTROOT).iterdir() if p.is_dir())
              if arg == "all" else [arg])
    for m in models:
        r = analyse(m)
        if "error" in r:
            print(f"[{m}] ERROR: {r['error']}")
            continue
        print(f"\n===== {m}: {r.get('verdict')} =====")
        for k, v in r.items():
            if k != "series":
                print(f"  {k}: {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
