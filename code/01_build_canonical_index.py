#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""01_build_canonical_index.py

Build the canonical evidence index and the quarantine/freeze manifest for the
HPV E6 project rebuild. Read-only: it never moves, edits or deletes project files.

Outputs (rebuild_20261025/results/):
  01_canonical_index.csv     assets allowed to support scientific claims (with tier)
  01_control_manifest.csv    quarantined / frozen / pending files with rule + reason
  01_index_summary.json      counts and sizes per tier
"""

from __future__ import annotations

import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Rules. (glob, tier, reason, action, evidence)
# tier: canonical_input | canonical_result | debug_only | quarantine | freeze | pending_review
# ---------------------------------------------------------------------------

QUARANTINE_RULES = [
    (
        "md_simulation_setup/_QUARANTINE_SYNTHETIC_2026-07-29/**/*",
        "quarantine",
        "Synthetic results isolated by the project itself",
        "Never cite; keep for provenance audit only",
        "md_simulation_setup/PROVENANCE_AND_VALIDATION.md",
    ),
    (
        "results/md_analysis/_QUARANTINE_SYNTHETIC_2026-07-29/**/*",
        "quarantine",
        "Synthetic results isolated by the project itself",
        "Never cite; keep for provenance audit only",
        "md_simulation_setup/PROVENANCE_AND_VALIDATION.md",
    ),
    (
        "results/massive_molecular_docking/19_ultra_massive*",
        "quarantine",
        "54,032 rows but only 1,000 compound IDs; MW/logP/TPSA conflict across duplicates",
        "Exclude from evidence chain; regenerate library with unique IDs",
        "results/massive_molecular_docking/19_ultra_massive_docking_results.csv",
    ),
    (
        "results/massive_molecular_docking/19_top_100_candidates.csv",
        "quarantine",
        "Derived from the invalid ultra-massive table",
        "Exclude until parent table is regenerated",
        "results/massive_molecular_docking/19_ultra_massive_docking_results.csv",
    ),
    (
        "md_simulation_setup/fig5_binding_mode.png",
        "quarantine",
        "Per-residue contributions are hard-coded in generate_figures.py:75-88",
        "Withdraw; regenerate only from canonical decomposition output",
        "md_simulation_setup/generate_figures.py:75-88",
    ),
    (
        "md_simulation_setup/generate_figures.py",
        "freeze",
        "Source of the hard-coded Fig.5 contributions",
        "Replace with canonical-result-driven figure script",
        "md_simulation_setup/generate_figures.py:75-88",
    ),
    (
        "md_simulation_setup/manuscript_draft.md",
        "freeze",
        "Contains historic 545/571/313 energies, three-site claims, PCA/FEL without valid evidence",
        "Do not submit or extend; rewrite after canonical results exist",
        "md_simulation_setup/PROVENANCE_AND_VALIDATION.md",
    ),
    (
        "md_simulation_setup/manuscript.tex",
        "freeze",
        "References Fig.1-6 that are not backed by canonical results",
        "Recompile only after figures are regenerated",
        "md_simulation_setup/manuscript.pdf",
    ),
    (
        "md_simulation_setup/manuscript.pdf",
        "freeze",
        "Only 3 pages, 0 figures; not a complete submission",
        "Recompile after canonical figures exist",
        "md_simulation_setup/manuscript.pdf",
    ),
    (
        "results/**/*.png",
        "pending_review",
        "Legacy analysis figures; provenance and version consistency unverified",
        "Redraw only from canonical result tables",
        "systematic review section 3.4",
    ),
    (
        "results/**/*.pdf",
        "pending_review",
        "Legacy analysis figures/reports; provenance unverified",
        "Redraw only from canonical result tables",
        "systematic review section 3.4",
    ),
]

CANONICAL_RULES = [
    (
        "results/10_扩展真实HPV_E6数据集.csv",
        "canonical_input",
        "Raw 878-row ML dataset; labels are type-level proxies, not verified variant phenotypes",
        "Use as cleaning input only; never publish raw-label performance",
        "878 rows, 50 types, 12 physicochemical features",
    ),
    (
        "rec_clean.pdb",
        "canonical_input",
        "Receptor template for docking/MD rebuild",
        "Reuse as the single receptor conformation",
        "local asset present",
    ),
    ("zn1.pdb", "canonical_input", "Zn site 1 template", "Reuse for Zn model rebuild", "local asset present"),
    ("zn2.pdb", "canonical_input", "Zn site 2 template", "Reuse for Zn model rebuild", "local asset present"),
    ("ZN.mol2", "canonical_input", "Zn parameters for LEaP", "Reuse; verify against non-bonded/bonded model choice", "local asset present"),
    (".workbuddy/4GIZ_E6_single_clean.pdb", "canonical_input", "Cleaned 4GIZ E6 structure", "Reuse as receptor source", "local asset present"),
    (".workbuddy/4GIZ_full.pdb", "canonical_input", "Full 4GIZ complex", "Reference only", "local asset present"),
    (
        "docking_results_full/4GIZ/CHEMBL*.pdb",
        "canonical_input",
        "14 ChEMBL candidate ligands (PDB)",
        "Use only with verified charge/protonation; prefer *_out.pdbqt poses",
        "14 ligands",
    ),
    (
        "docking_results_full/4GIZ/CHEMBL*_out.pdbqt",
        "canonical_input",
        "Vina docked poses for the 14 candidates",
        "MUST be the MD/MM-PBSA start pose (previous run used ligand.pdb instead)",
        "Vina output",
    ),
    (
        "candidates_run_final_20260818.tar.gz",
        "debug_only",
        "14 candidate PB endpoint results; reproducible arithmetic but invalid start pose and Zn",
        "Debug/negative example only; do not rank binders from it",
        "0/714 frames keep both Zn sites tetra-coordinated",
    ),
    (
        "candidates_run_backup_20260817.tar.gz",
        "debug_only",
        "Earlier failed run history",
        "Version comparison only; never mix with final archive",
        "archive backup",
    ),
    (
        ".workbuddy/SUMMARY_final.csv",
        "debug_only",
        "Manual recomputation matching the archive components",
        "Use to regression-test the parser, not as a scientific ranking",
        "component sum VDWAALS+EEL+EPB+ENPOLAR+EDISPER",
    ),
    (
        "models/hpv_e6_model_metadata.json",
        "debug_only",
        "Current SVC metrics; arithmetically consistent but validation is biased",
        "Reference baseline only; supersede by grouped nested validation",
        "Accuracy 0.976082 / AUC 0.995488",
    ),
]


def stat_row(path: Path):
    st = path.stat()
    return {
        "size_bytes": st.st_size,
        "mtime": __import__("datetime").datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds"),
    }


def collect():
    canonical, controlled = [], []
    seen = set()

    for pattern, tier, reason, action, evidence in CANONICAL_RULES:
        for p in sorted(ROOT.glob(pattern)):
            if p.is_dir() or str(p) in seen:
                continue
            seen.add(str(p))
            canonical.append(
                {
                    "path": p.relative_to(ROOT).as_posix(),
                    "tier": tier,
                    "role": reason,
                    "caveat_or_usage": action,
                    "evidence": evidence,
                    **stat_row(p),
                }
            )

    for pattern, tier, reason, action, evidence in QUARANTINE_RULES:
        for p in sorted(ROOT.glob(pattern)):
            if p.is_dir() or str(p) in seen:
                continue
            seen.add(str(p))
            controlled.append(
                {
                    "path": p.relative_to(ROOT).as_posix(),
                    "tier": tier,
                    "rule": pattern,
                    "reason": reason,
                    "action": action,
                    "evidence": evidence,
                    **stat_row(p),
                }
            )

    return canonical, controlled


def write_csv(rows, path, fieldnames):
    import csv

    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def main():
    canonical, controlled = collect()

    write_csv(
        canonical,
        OUT / "01_canonical_index.csv",
        ["path", "tier", "role", "caveat_or_usage", "evidence", "size_bytes", "mtime"],
    )
    write_csv(
        controlled,
        OUT / "01_control_manifest.csv",
        ["path", "tier", "rule", "reason", "action", "evidence", "size_bytes", "mtime"],
    )

    summary = {
        "root": str(ROOT),
        "canonical_files": len(canonical),
        "controlled_files": len(controlled),
        "by_tier_canonical": {},
        "by_tier_controlled": {},
    }
    for r in canonical:
        summary["by_tier_canonical"][r["tier"]] = summary["by_tier_canonical"].get(r["tier"], 0) + 1
    for r in controlled:
        summary["by_tier_controlled"][r["tier"]] = summary["by_tier_controlled"].get(r["tier"], 0) + 1

    with open(OUT / "01_index_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(f"canonical files: {len(canonical)}")
    print(f"controlled files: {len(controlled)}")
    print(json.dumps(summary["by_tier_controlled"], ensure_ascii=False))
    print(f"written: {OUT}")


if __name__ == "__main__":
    main()
