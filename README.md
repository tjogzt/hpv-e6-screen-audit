# Auditing a single-structure virtual screen of the HPV16 E6 zinc-finger pocket

Code, data and results archive for the manuscript:

> *Auditing a single-structure virtual screen of the HPV16 E6 zinc-finger pocket:
> a negative verdict, control-calibrated reporting, and the lineage confound in
> E6 sequence analysis* (revised manuscript, Round 1)

**Contents:** an audited negative docking-to-MD screen of 14 ChEMBL candidates
against the HPV16 E6 E6AP-interface pocket (PDB 4GIZ-derived receptor, bonded
ZAFF/CYM zinc model, DMSO + E6AP LxxLL controls, 3 × 20 ns replicas per complex,
MM-PBSA re-scoring, exploratory GaMD check), plus a lineage-confound analysis of
878 HPV E6 sequence records (45 types, type-level statistics, type-grouped CV,
selection-free dual-zinc conservation scan).

## Repository layout

| Directory | Contents |
|---|---|
| `code/` | All analysis and pipeline scripts (numbered in execution order; `48_revision_numbers.py` regenerates every number quoted in the manuscript from the CSVs) |
| `manuscript/` | Revised manuscript (`merged_manuscript_v2.md`) and figures (`figures/Fig*.pdf`) |
| `results/` | Per-analysis outputs (CSV/JSON/MD), batch-2 per-run summaries (`batch2/`), GaMD analysis outputs (`gamd_v2/`) |
| `templates/` | Receptor/zinc structural templates, the machine-readable **reconstructed ZAFF frcmod**, the non-bonded fallback frcmod, and a runnable tleap input (`build_system.leap.in`) that rebuilds the solvated system |
| `ligands/` | The 74-compound ChEMBL library CSV, the 14-candidate provenance table (SMILES, properties, Vina scores) |
| `docking/` | Library-build and batch-docking scripts, the candidates' Vina inputs/outputs for both receptors (`4GIZ_poses/` used; `4XR8_poses/` excluded and documented), and `docking_progress.json` — the decisive provenance evidence for the 14-compound subset |
| `data/partA/` | The 878-record E6 sequence dataset (English-named copy `10_extended_real_HPV_E6_dataset.csv`; original filename preserved alongside, see `FILENAME_MAP.md`) |

## Reproducing the headline numbers

Every number in the manuscript's results sections is generated from the archived
CSVs at write time — nothing is transcribed by hand:

```bash
python3 code/48_revision_numbers.py     # all Part A / revision numbers -> results/48_revision_numbers.json
python3 code/49_candidate_provenance.py # candidate-set provenance table
python3 code/50_tier_floor_bootstrap.py # control-tier threshold bootstrap
```

Part A statistics: `code/02_clean_ml_taxonomy.py` → `08_clade_signature_descriptive.py`
(both units of analysis, both CV constructions) → `46_k_robustness.py` →
`47_selection_free_zinc_scan.py`.

Part B: the MD/MM-PBSA pipeline ran on a GPU instance (see `SOFTWARE_VERSIONS.md`,
`TRAJECTORIES.md`); per-run outputs are archived in `results/batch2/` and aggregated
by `code/36_pull_full_screen.py` → `results/36_screen_full.csv`. The solvated receptor
system can be rebuilt from `templates/` with the archived tleap input.

## Documentation files

- `SOFTWARE_VERSIONS.md` — software stack with versions (and honest "not recorded" entries)
- `DATA_PROVENANCE.md` — sources, access dates and query strings for every dataset
- `TRAJECTORIES.md` — what trajectory data survives, what was lost, and how to regenerate
- `FILENAME_MAP.md` — mapping of legacy Chinese filenames to English names
- `MANIFEST.sha256` — SHA-256 checksums of every archived file

## License

- **Code** (`code/`, `docking/`): MIT — see `LICENSE`
- **Data, results, manuscript, templates**: CC BY 4.0 — see `LICENSE-CC-BY-4.0.txt`

## Citation

[DOI will appear here after Zenodo archiving of the v1.0.0 release.]

## Contact

Correspondence: [to be completed by the authors]
