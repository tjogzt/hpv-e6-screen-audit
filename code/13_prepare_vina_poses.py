#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""13_prepare_vina_poses.py — 为全部 14 个 ChEMBL 配体提取 Vina 最优姿态（MODEL 1）。

为什么必须有这一步
------------------
07 号脚本已证明：旧管线上传的起始结构是 `docking_results_full/4GIZ/<LIG>.pdb`（未对接的
原始输入），与 Vina 参考姿态质心偏移 72.0 Å —— 这是旧 MM-PBSA 数据作废的主因之一。
但 07 只处理了 CHEMBL46917 一个配体；其余 13 个配体若直接沿用原始 pdb，重跑会把同一个
错误再犯 13 遍。本脚本把 14 个配体全部转成"已对接的 MODEL 1 姿态"，并逐配体记录：

  - Vina 打分（REMARK VINA RESULT）
  - 与未对接原始 pdb 的质心偏移（Å）—— 证明这批姿态确实位于受体口袋坐标系
  - 原子数（pdbqt 姿态为极性氢模型，原子数通常少于原始 pdb）

输出
----
  data/poses/<LIG>_pose1.pdb        —— 上传远端作为 MD 起始结构
  results/13_vina_poses.csv         —— 逐配体审计表
  results/13_vina_poses.json        —— 汇总（含缺失/失败清单）
"""
from __future__ import annotations

import csv
import json
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REBUILD = Path(__file__).resolve().parents[1]
DOCK_DIR = ROOT / "docking_results_full" / "4GIZ"
POSES = REBUILD / "data" / "poses"
RESULTS = REBUILD / "results"
POSES.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)

OBABEL = shutil.which("obabel")
if not OBABEL:
    raise SystemExit("未找到 obabel（openbabel）。请先安装：brew install open-babel")


def extract_model1(pdbqt: Path, dst: Path) -> int:
    """抽取 MODEL 1（Vina 最优姿态），返回模型中姿态总数。"""
    lines, inside, n_models = [], False, 0
    for line in pdbqt.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            n_models += 1
            inside = (n_models == 1)
            continue
        if line.startswith("ENDMDL"):
            inside = False
            continue
        if inside and line.startswith(("ATOM", "HETATM")):
            lines.append(line)
    dst.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return n_models


def coords(pdb: Path) -> np.ndarray:
    xyz = []
    for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
    return np.asarray(xyz, dtype=float)


def main() -> int:
    pdbqts = sorted(DOCK_DIR.glob("*_out.pdbqt"))
    if not pdbqts:
        raise SystemExit(f"未找到任何 *_out.pdbqt：{DOCK_DIR}")
    print(f"发现 {len(pdbqts)} 个 Vina 输出文件")

    rows, failures = [], []
    tmp = Path(tempfile.mkdtemp())
    t0 = time.time()

    for pdbqt in pdbqts:
        name = pdbqt.name.replace("_out.pdbqt", "")
        undocked = DOCK_DIR / f"{name}.pdb"
        model1 = tmp / f"{name}_model1.pdbqt"
        dst = POSES / f"{name}_pose1.pdb"
        try:
            n_models = extract_model1(pdbqt, model1)
            if n_models == 0:
                raise RuntimeError("pdbqt 中无 MODEL 段")
            subprocess.run([OBABEL, "-ipdbqt", str(model1), "-opdb", "-O", str(dst)],
                           check=True, capture_output=True, text=True)
            scores = [float(m.group(1)) for m in re.finditer(
                r"REMARK VINA RESULT:\s+(-?\d+\.\d+)",
                pdbqt.read_text(encoding="utf-8", errors="replace"))]
            xyz_pose = coords(dst)
            xyz_und = coords(undocked) if undocked.exists() else np.empty((0, 3))
            if len(xyz_pose) == 0:
                raise RuntimeError("姿态转换后无原子")
            shift = (float(np.linalg.norm(xyz_pose.mean(0) - xyz_und.mean(0)))
                     if len(xyz_und) else None)
            rows.append({
                "ligand": name,
                "n_poses_in_pdbqt": n_models,
                "vina_best_kcal_mol": scores[0] if scores else "",
                "pose_atoms": int(len(xyz_pose)),
                "undocked_atoms": int(len(xyz_und)),
                "centroid_shift_A": "" if shift is None else round(shift, 3),
                "pose1_pdb": str(dst.relative_to(ROOT)),
            })
            print(f"  [OK] {name}: 姿态 {n_models} 个，最佳 {rows[-1]['vina_best_kcal_mol']} kcal/mol，"
                  f"原子 {len(xyz_pose)}，质心偏移 {rows[-1]['centroid_shift_A']} Å")
        except Exception as exc:  # noqa: BLE001
            failures.append({"ligand": name, "error": str(exc)[:200]})
            print(f"  [FAIL] {name}: {exc}")

    with open(RESULTS / "13_vina_poses.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["ligand", "n_poses_in_pdbqt", "vina_best_kcal_mol",
                                           "pose_atoms", "undocked_atoms",
                                           "centroid_shift_A", "pose1_pdb"])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    shifts = [r["centroid_shift_A"] for r in rows if r["centroid_shift_A"] != ""]
    summary = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "n_ligands": len(rows),
        "n_failures": len(failures),
        "failures": failures,
        "centroid_shift_A_median": round(float(np.median(shifts)), 3) if shifts else None,
        "centroid_shift_A_max": round(float(np.max(shifts)), 3) if shifts else None,
        "note": ("全部配体现使用 Vina MODEL 1 姿态作为 MD 起始结构；"
                 "质心偏移列为该姿态与未对接原始 pdb 的距离，用于审计"
                 "（旧管线即因使用未对接 pdb 而产生 72 Å 偏移）。"
                 "pdbqt 姿态为极性氢模型，原子数少于原始 pdb 属正常。"),
        "upload_target": "远端 $BASE/ligands_vina/*.pdb（与 LIGDIR 一致）",
    }
    with open(RESULTS / "13_vina_poses.json", "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(f"\n完成 {len(rows)} 个，失败 {len(failures)} 个，"
          f"质心偏移中位 {summary['centroid_shift_A_median']} Å / 最大 {summary['centroid_shift_A_max']} Å")
    print(f"输出：{RESULTS / '13_vina_poses.csv'}")
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
