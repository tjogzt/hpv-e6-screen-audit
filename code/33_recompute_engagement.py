#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
33_recompute_engagement.py — 从坐标重算起始 pose 与受体的最小重原子距离（替换前/后）

为什么需要：results/20_start_pose_engagement_audit.csv 只保留了**替换后**（selected pose）
的距离，三个被替换配体的 MODEL 1 原始距离未留在 CSV 中，而稿件 §3.6 / Figure 3A 引用了
"替换前 9.16 / 10.05 / 9.24 Å"。为避免引用无据，本脚本直接从坐标重算并交叉核对。

重要陷阱（本机已踩过）：**不能用 Bio.PDB.PDBParser 解析这些配体 PDB**。Open Babel 输出的
配体原子名全部为 "C"/"O"（每个原子同名），而 Biopython 的 Residue 以原子名建索引，
会把 28 个原子折叠成 2 个（实测 CHEMBL717：28 → 2），距离随之失真。
本脚本统一使用**逐行解析**，对受体与配体一视同仁。

输入：
  data/pilot/rec_clean.pdb                  —— 管线实际使用的受体（4GIZ 衍生，残基 2–143）
  data/poses/<LIG>_pose1.pdb                —— 当前使用的（selected）pose
  <项目>/docking_results_full/4GIZ/<LIG>_out.pdbqt —— Vina 原始输出，含 9 个 MODEL

输出：
  results/33_engagement_recomputed.csv
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REBUILD = HERE.parent
RESULTS = REBUILD / "results"
PROJ = REBUILD.parent

RECEPTOR = REBUILD / "data" / "pilot" / "rec_clean.pdb"
POSE_DIR = REBUILD / "data" / "poses"
DOCK = PROJ / "docking_results_full" / "4GIZ"
AUDIT = RESULTS / "20_start_pose_engagement_audit.csv"

HYDROGEN_NAMES = {"H", "HD", "HS", "HA", "HB", "HG", "HH", "HZ", "HD1", "HD2"}


def parse_heavy_atoms(path: Path) -> np.ndarray:
    """逐行解析 PDB 的重原子坐标（忽略元素列缺失的情况，按原子名判氢）。"""
    xyz = []
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        # 原子名去掉前导数字（如 "1HB" -> "HB"）
        stripped = name.lstrip("0123456789")
        if stripped in HYDROGEN_NAMES or name in HYDROGEN_NAMES:
            continue
        try:
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except ValueError:
            continue
    return np.array(xyz, dtype=float)


def parse_models_pdbqt(path: Path) -> list[np.ndarray]:
    """逐行解析 pdbqt 中每个 MODEL 的重原子坐标（按 MODEL 顺序）。"""
    models, cur = [], []
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            stripped = name.lstrip("0123456789")
            if stripped in HYDROGEN_NAMES or name in HYDROGEN_NAMES:
                continue
            cur.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        elif line.startswith("ENDMDL"):
            models.append(np.array(cur, dtype=float))
            cur = []
    return models


def min_dist(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)
    return float(d.min())


def main() -> None:
    stopifnot = lambda c, m: (_ for _ in ()).throw(AssertionError(m)) if not c else None
    stopifnot(RECEPTOR.exists(), f"缺少受体: {RECEPTOR}")

    rec = parse_heavy_atoms(RECEPTOR)
    print(f"受体 {RECEPTOR.name}: {len(rec)} 个重原子")

    audit = {r["ligand"]: r for r in csv.DictReader(open(AUDIT))}
    rows = []
    for lig in sorted(audit):
        a = audit[lig]
        base = lig.replace("_pose1", "")
        pdb_file = POSE_DIR / f"{base}_pose1.pdb"
        pdbqt = DOCK / f"{base}_out.pdbqt"
        replaced = a["pose_replaced"].strip() == "yes"
        sel = int(a["selected_model"])

        models = parse_models_pdbqt(pdbqt) if pdbqt.exists() else []
        d_m1 = min_dist(models[0], rec) if models else float("nan")
        d_sel_model = min_dist(models[sel - 1], rec) if len(models) >= sel else float("nan")
        pdb_xyz = parse_heavy_atoms(pdb_file) if pdb_file.exists() else np.array([])
        d_pdb = min_dist(pdb_xyz, rec) if len(pdb_xyz) else float("nan")

        rows.append({
            "ligand": lig,
            "n_atoms_pose1_pdb": len(pdb_xyz),
            "n_atoms_model1": len(models[0]) if models else 0,
            "d_model1_A": round(d_m1, 2),
            "d_selected_model_A": round(d_sel_model, 2),
            "d_pose1_pdb_A": round(d_pdb, 2),
            "selected_model": sel,
            "pose_replaced": replaced,
            "audit_csv_dist_A": float(a["min_heavy_dist_to_receptor_A"]),
        })

    print(f"\n{'ligand':<24} {'nat(pdb/m1)':>11} {'pdb文件':>7} {'MODEL1':>7} "
          f"{'selected':>8} {'CSV':>6} {'pdb==CSV':>8} {'替换':>5}")
    for r in rows:
        same = abs(r["d_pose1_pdb_A"] - r["audit_csv_dist_A"]) < 0.05
        print(f"{r['ligand']:<24} {str(r['n_atoms_pose1_pdb'])+'/'+str(r['n_atoms_model1']):>11} "
              f"{r['d_pose1_pdb_A']:>7.2f} {r['d_model1_A']:>7.2f} "
              f"{r['d_selected_model_A']:>8.2f} {r['audit_csv_dist_A']:>6.2f} "
              f"{'yes' if same else 'NO':>8} "
              f"{('MODEL%d' % r['selected_model']) if r['pose_replaced'] else '-':>5}")

    # 一致性判定：pose1.pdb 是否与 CSV 完全吻合（即 CSV 就是用该文件测的）
    agree = [r for r in rows if abs(r["d_pose1_pdb_A"] - r["audit_csv_dist_A"]) < 0.05]
    print(f"\npose1.pdb 距离与 CSV 一致的: {len(agree)}/{len(rows)}")

    # 三个被替换配体的 MODEL 1 原始距离
    repl = sorted([r for r in rows if r["pose_replaced"]],
                  key=lambda x: -x["d_model1_A"])
    print("\n=== 三个被替换配体（MODEL 1 = 替换前）===")
    for r in repl:
        print(f"  {r['ligand']:<24} 前 {r['d_model1_A']:.2f} Å  ->  "
              f"后(MODEL{r['selected_model']}) {r['d_selected_model_A']:.2f} Å")

    # 与稿件现有措辞对照
    ms = sorted([9.16, 10.05, 9.24])
    got = sorted([r["d_model1_A"] for r in repl])
    print(f"\n稿件 §3.6 现写: {ms}")
    print(f"本脚本实测   : {got}")
    print(f"一致: {ms == got}")

    out = RESULTS / "33_engagement_recomputed.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n[saved] {out}")


if __name__ == "__main__":
    main()
