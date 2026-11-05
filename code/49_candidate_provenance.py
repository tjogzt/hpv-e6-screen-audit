#!/usr/bin/env python3
"""
49_candidate_provenance.py — R1 (Round-1 revision): candidate-set provenance table.

Reconstructs, from surviving primary artifacts only (no hand-transcribed numbers),
the provenance of the 14 ChEMBL candidate ligands used in Part B:

  1. ligands/01_candidate_ligand_library.csv  — the 74-compound library built by
     scripts/02_fetch_ligand_library.py (ChEMBL web API, 20 indication keywords,
     canonical-SMILES dedup, Lipinski+Veber, Morgan FP Tanimoto < 0.8 diversity).
  2. docking_progress.json (2025-10-07)       — proves the MD candidates are the
     FIRST entries of that library file, in file order.
  3. docking_results_full/4GIZ/<ID>_out.pdbqt — Vina best-mode affinity per ligand.

Outputs:
  results/49_candidate_provenance_table.csv  — 14 rows: ID, name, SMILES,
      properties, Lipinski/Veber, best Vina affinity, library row index.
  results/49_candidate_provenance.json       — library-level provenance record
      (queries, filters, dates read from file mtimes, selection evidence).

Run: /opt/anaconda3/bin/python3 code/49_candidate_provenance.py
"""
import json
import os
import re
from datetime import datetime
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]          # .../rebuild_20261025
PROJ = ROOT.parent                                  # .../HPV_structure
LIB = PROJ / "ligands" / "01_candidate_ligand_library.csv"
PROGRESS = PROJ / "results" / "molecular_docking" / "docking_progress.json"
FETCH_SCRIPT = PROJ / "scripts" / "02_fetch_ligand_library.py"
DOCK_DIR = PROJ / "docking_results_full" / "4GIZ"
OUT_CSV = ROOT / "results" / "49_candidate_provenance_table.csv"
OUT_JSON = ROOT / "results" / "49_candidate_provenance.json"

MD14 = [
    "CHEMBL283639", "CHEMBL29609", "CHEMBL46917", "CHEMBL1129", "CHEMBL1090",
    "CHEMBL184", "CHEMBL118218", "CHEMBL320775", "CHEMBL4650485", "CHEMBL119",
    "CHEMBL235668", "CHEMBL3701238", "CHEMBL848", "CHEMBL717",
]


def best_vina_affinity(pdbqt_path: Path) -> float:
    """First REMARK VINA RESULT line = best mode affinity (kcal/mol)."""
    with open(pdbqt_path) as fh:
        for line in fh:
            if line.startswith("REMARK VINA RESULT"):
                return float(line.split()[3])
    raise ValueError(f"no VINA RESULT in {pdbqt_path}")


def parse_queries_from_fetch_script() -> list:
    """Extract the QUERIES list literally from the archived fetch script."""
    text = FETCH_SCRIPT.read_text(encoding="utf-8")
    m = re.search(r"QUERIES\s*=\s*\[(.*?)\]", text, re.S)
    if not m:
        raise ValueError("QUERIES list not found in 02_fetch_ligand_library.py")
    return re.findall(r"'([^']+)'", m.group(1))


def mtime(p: Path) -> str:
    return datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d")


def main() -> None:
    lib = pd.read_csv(LIB)
    assert len(lib) == 74, f"expected 74-row library, got {len(lib)}"
    # The 14 candidates must be exactly rows 0..13 (file order).
    first14 = lib["id"].head(14).tolist()
    assert first14 == MD14, f"library head mismatch: {first14}"

    # Selection evidence: docking_progress.json entries must match rows 0..n.
    prog = json.loads(PROGRESS.read_text())
    prog_ids = [e["ligand_id"] for e in prog]
    assert prog_ids == lib["id"].head(len(prog_ids)).tolist(), "progress order mismatch"
    assert set(prog_ids) <= set(MD14)

    rows = []
    for i, cid in enumerate(MD14):
        r = lib[lib["id"] == cid].iloc[0]
        out_pdbqt = DOCK_DIR / f"{cid}_out.pdbqt"
        assert out_pdbqt.exists(), f"missing {out_pdbqt}"
        rows.append({
            "library_row": i,
            "chembl_id": cid,
            "pref_name": r["name"],
            "smiles": r["smiles"],
            "MW": round(float(r["MW"]), 2),
            "LogP": round(float(r["LogP"]), 2),
            "HBD": int(r["HBD"]),
            "HBA": int(r["HBA"]),
            "RotBonds": int(r["RotBonds"]),
            "TPSA": round(float(r["TPSA"]), 2),
            "LipinskiPass": bool(r["LipinskiPass"]),
            "VeberPass": bool(r["VeberPass"]),
            "best_vina_4GIZ_kcal_mol": best_vina_affinity(out_pdbqt),
        })
    table = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(exist_ok=True)
    table.to_csv(OUT_CSV, index=False)

    record = {
        "library_file": str(LIB.relative_to(PROJ)),
        "library_rows": int(len(lib)),
        "library_built_by": "scripts/02_fetch_ligand_library.py",
        "library_script_mtime": mtime(FETCH_SCRIPT),
        "library_file_mtime": mtime(LIB),
        "chembl_api": "https://www.ebi.ac.uk/chembl/api/data/molecule/search",
        "chembl_release": "not recorded at access time",
        "queries": parse_queries_from_fetch_script(),
        "n_queries": 20,
        "filters": {
            "dedup": "canonical SMILES",
            "lipinski": "MW<=500, LogP<=5, HBD<=5, HBA<=10",
            "veber": "RotBonds<=10, TPSA<=140",
            "diversity": "Morgan FP r=2 2048-bit, greedy, max pairwise Tanimoto < 0.8",
        },
        "selection_rule": (
            "first 14 entries of the library file, in file order "
            "(positional truncation; no score-, activity- or structure-based "
            "rule recorded)"
        ),
        "selection_evidence": {
            "file": "results/molecular_docking/docking_progress.json",
            "n_entries": len(prog_ids),
            "entries_match_library_rows": "0..%d" % (len(prog_ids) - 1),
            "timestamp_first": prog[0]["timestamp"],
        },
        "vina_scores_source": "docking_results_full/4GIZ/<ID>_out.pdbqt (REMARK VINA RESULT)",
        "vina_range": [
            float(table.best_vina_4GIZ_kcal_mol.min()),
            float(table.best_vina_4GIZ_kcal_mol.max()),
        ],
        "table": str(OUT_CSV.relative_to(ROOT)),
        "generated": datetime.now().isoformat(timespec="seconds"),
    }
    OUT_JSON.write_text(json.dumps(record, indent=2))

    print(f"[49] wrote {OUT_CSV} ({len(table)} rows)")
    print(f"[49] wrote {OUT_JSON}")
    print(f"[49] Vina range: {record['vina_range'][0]:.3f} .. {record['vina_range'][1]:.3f}")
    print("[49] queries:", ", ".join(record["queries"]))


if __name__ == "__main__":
    main()
