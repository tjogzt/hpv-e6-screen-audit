# Filename map (legacy → English)

The working project accumulated some non-English (Chinese) filenames during
development. For this archive, one original filename is preserved alongside
its English-named byte-identical copy so that (a) the archive is portable
across systems with poor Unicode support and (b) the link back to the
original working files is not lost.

| Archived English name | Original (preserved) name | Bytes | Note |
|---|---|---|---|
| `data/partA/10_extended_real_HPV_E6_dataset.csv` | `data/partA/10_扩展真实HPV_E6数据集.csv` | identical (SHA-256 verified in `MANIFEST.sha256`) | the 878-record HPV E6 sequence dataset used throughout Part A |

All other archived files kept their original (already ASCII) working names;
no renaming was applied, so script-internal cross-references (e.g. paths
quoted in `code/*.py` docstrings) remain traceable to the original project
layout under `rebuild_20261025/`.

## Layout mapping (original project → this archive)

| Original location (working project) | Archive location |
|---|---|
| `rebuild_20261025/code/*.py` / `*.R` | `code/` |
| `rebuild_20261025/manuscript/merged_manuscript_v2.md` | `manuscript/` |
| `rebuild_20261025/results/Fig*.pdf` | `manuscript/figures/` (and byte-identical copies kept in `results/`) |
| `rebuild_20261025/results/*.csv|json|md` | `results/` |
| `rebuild_20261025/results/batch2/` | `results/batch2/` |
| `rebuild_20261025/results/gamd_v2/` | `results/gamd_v2/` |
| `rebuild_20261025/ligands/01_candidate_ligand_library.csv` | `ligands/` |
| `rebuild_20261025/code/02_fetch_ligand_library.py` (ligand-library stage) | `docking/02_fetch_ligand_library.py` |
| `docking_results_full/4GIZ/` (14 candidate pose sets) | `docking/4GIZ_poses/` |
| `docking_results_full/4XR8/` (excluded-receptor control) | `docking/4XR8_poses/` |
| `results/molecular_docking/docking_progress.json` | `docking/docking_progress.json` |
| `10_扩展真实HPV_E6数据集.csv` (project root area) | `data/partA/` (both names) |
| `.workbuddy/` receptor/zinc templates (`rec_clean.pdb`, `zn1.pdb`, `zn2.pdb`, `ZN.mol2`, `peptide_LxxLL.pdb`, `zinc_nonbond_fallback.frcmod`) | `templates/` |
| ZAFF parameters reconstructed from `code/24_zn_model_test.py` | `templates/zaff_reconstructed.frcmod` (machine-readable, header documents reconstruction) |

Files that were never pulled back from the compute instances (trajectories,
`.tpr`/`.cpt`, batch-1 per-run tables, smoke test) are not in this archive —
see `TRAJECTORIES.md`.
