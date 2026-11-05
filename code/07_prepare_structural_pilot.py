#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""07_prepare_structural_pilot.py

Prepare the corrected structural pilot inputs (P0 of the structural line) and
document, with numbers, why the previous 14-candidate run cannot be used.

What it does locally (no GPU / no MD):
  1. extracts MODEL 1 of the Vina docked pose (CHEMBL46917_out.pdbqt) and converts
     it to PDB with obabel -> data/pilot/CHEMBL46917_pose1.pdb
  2. compares the docked pose with the undocked docking_results_full/4GIZ/CHEMBL46917.pdb
     that was uploaded as the previous MD starting structure (centroid shift + RMSD)
  3. inventories receptor / Zn templates (rec_clean.pdb, zn1.pdb, zn2.pdb, ZN.mol2)
  4. writes the pilot input checklist (results/07_structural_pilot_prep.json)

It does NOT run MD or MM-PBSA; those need the GPU host.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REBUILD = Path(__file__).resolve().parents[1]
OUT = REBUILD / "results"
PILOT = REBUILD / "data" / "pilot"
OUT.mkdir(parents=True, exist_ok=True)
PILOT.mkdir(parents=True, exist_ok=True)

LIG = "CHEMBL46917"
DOCK_DIR = ROOT / "docking_results_full" / "4GIZ"
POSE_PDBQT = DOCK_DIR / f"{LIG}_out.pdbqt"
UNDOCKED_PDB = DOCK_DIR / f"{LIG}.pdb"
OBABEL = shutil.which("obabel")


def extract_model1(pdbqt: Path, dst: Path) -> int:
    lines, inside = [], False
    n_models = 0
    for line in pdbqt.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            n_models += 1
            inside = (n_models == 1)
            continue
        if line.startswith("ENDMDL"):
            inside = False
            continue
        if inside and (line.startswith(("ATOM", "HETATM"))):
            lines.append(line)
    dst.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return n_models


def coords(pdb: Path):
    xyz, elements = [], []
    for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            elements.append(line[76:78].strip() or line[12:16].strip()[:1])
    return np.asarray(xyz), elements


def kabsch_rmsd(a: np.ndarray, b: np.ndarray):
    """RMSD after optimal superposition (Kabsch) on the first min(len) atoms."""
    n = min(len(a), len(b))
    if n == 0:
        return None
    a, b = a[:n], b[:n]
    ac, bc = a - a.mean(0), b - b.mean(0)
    h = ac.T @ bc
    u, s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1, 1, d]) @ u.T
    ar = ac @ rot.T
    return float(np.sqrt(((ar - bc) ** 2).sum(axis=1).mean()))


def zn_inventory():
    inv = {}
    for name in ("rec_clean.pdb", "zn1.pdb", "zn2.pdb", ".workbuddy/_remote_rec_clean.pdb",
                 ".workbuddy/4GIZ_E6_single_clean.pdb", ".workbuddy/4GIZ_full.pdb"):
        p = ROOT / name
        if not p.exists():
            inv[name] = {"exists": False}
            continue
        txt = p.read_text(encoding="utf-8", errors="replace")
        zn = len([l for l in txt.splitlines()
                  if l.startswith(("ATOM", "HETATM")) and l[12:16].strip().upper() == "ZN"])
        cys = len([l for l in txt.splitlines()
                   if l.startswith(("ATOM", "HETATM")) and l[17:20].strip().upper() == "CYS"])
        inv[name] = {"exists": True, "zn_atoms": zn, "cys_residues": cys}
    for cand in ("ZN.mol2", ".workbuddy/ZN.mol2"):
        mol2 = ROOT / cand
        inv[cand] = {"exists": mol2.exists(),
                     "bytes": mol2.stat().st_size if mol2.exists() else 0}
    return inv


def main():
    if not POSE_PDBQT.exists():
        raise SystemExit(f"missing docked pose: {POSE_PDBQT}")
    if not OBABEL:
        raise SystemExit("obabel not found; install openbabel to generate the pilot pose")

    tmp = Path(tempfile.mkdtemp())
    model1 = tmp / f"{LIG}_model1.pdbqt"
    n_models = extract_model1(POSE_PDBQT, model1)
    pose_pdb = PILOT / f"{LIG}_pose1.pdb"
    subprocess.run([OBABEL, "-ipdbqt", str(model1), "-opdb", "-O", str(pose_pdb)],
                   check=True, capture_output=True, text=True)

    vina_scores = [float(m.group(1)) for m in
                   re.finditer(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)",
                               POSE_PDBQT.read_text(encoding="utf-8", errors="replace"))]

    xyz_pose, el_pose = coords(pose_pdb)
    xyz_und, el_und = coords(UNDOCKED_PDB) if UNDOCKED_PDB.exists() else (np.empty((0, 3)), [])
    n = min(len(xyz_pose), len(xyz_und))
    centroid_shift = (float(np.linalg.norm(xyz_pose.mean(0) - xyz_und.mean(0)))
                      if n else None)
    rmsd = kabsch_rmsd(xyz_pose, xyz_und) if n else None

    summary = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "ligand": LIG,
        "docked_poses_in_out_pdbqt": n_models,
        "vina_scores_kcal_mol": vina_scores[:5],
        "pose_pdb_written": str(pose_pdb.relative_to(ROOT)),
        "docked_pose_atoms": int(len(xyz_pose)),
        "undocked_uploaded_atoms": int(len(xyz_und)),
        "atom_count_match": bool(len(xyz_pose) == len(xyz_und)),
        "centroid_shift_A": None if centroid_shift is None else round(centroid_shift, 3),
        "rmsd_after_superposition_A": None if rmsd is None else round(rmsd, 3),
        "interpretation": (
            "The previous MM-PBSA run used docking_results_full/4GIZ/<LIG>.pdb (the undocked "
            "input, not a Vina pose) as the MD starting structure. A 72 A centroid shift shows the "
            "uploaded structure is not even in the receptor frame of the docked pose. Caveat: the "
            "PDBQT pose carries polar hydrogens only (16 atoms) while the uploaded PDB has 26, so "
            "the RMSD below is computed on the first min(N) atoms in file order, not on a verified "
            "atom mapping; the centroid shift, not the RMSD, is the decisive evidence."
        ),
        "zn_and_receptor_inventory": zn_inventory(),
        "pilot_checklist": [
            "1 receptor conformation only (rec_clean.pdb, fixed protonation)",
            "1-2 Vina poses for 1 candidate (data/pilot/<LIG>_pose1.pdb)",
            "1 Zn model decision recorded up-front: non-bonded (Li/Merz) vs bonded/dummy; "
            "if the scientific target is chelation, a fixed tetra-coordinated model is not evidence",
            "fix LEaP atom order: receptor -> Zn ions -> ligand, and record it in the run log",
            "verify ligand charge/protonation before tleap (antechamber/resp)",
            "acceptance gate: both Zn sites stay 4-coordinated and the initial contact persists "
            "in a short run; otherwise stop and fix the model, do not extend sampling",
        ],
    }
    with open(OUT / "07_structural_pilot_prep.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
