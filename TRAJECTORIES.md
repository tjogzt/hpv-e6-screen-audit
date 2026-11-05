# Trajectory data: what survives, what was lost, how to regenerate

## Status: the raw trajectories are LOST

All molecular-dynamics trajectories (conventional MD and GaMD) were produced
on rented cloud GPU instances whose local disks have since been released.
No `.xtc`, `.gro` (production frames), `.cpt`, or `.tpr` files survive
locally. What survives is **everything the instance pushed back** — the
per-run analysis summaries and per-frame metric tables listed below — plus
the complete input-side state (structures, parameters, run scripts, run
configurations) needed to regenerate the trajectories from scratch.

**The manuscript's Round-1 revision states this explicitly** (§3.9, §6):
the GaMD reading is reported as a qualitative, non-reproducible-from-archive
stability check, and the loss is listed among the limitations.

## What IS archived

### Conventional MD + MM-PBSA screen (batch 2: 14 candidates + DMSO + LxxLL, 3 × 20 ns each)

- Per-run MM-PBSA summary tables: `results/batch2/SUMMARY_v3_<ligand>_rep<1-3>.csv`
  (48 runs + controls) and the combined `results/batch2/SUMMARY_v3_all.csv`
- Aggregated screen: `results/36_screen_full.csv` /
  `results/36_screen_byligand.csv` / `results/36_screen_report.md`
  (built by `code/36_pull_full_screen.py`)
- Exploratory rescue re-scoring: `results/39_rescue_deltaH.csv`
- Zn-model acceptance record: `results/zn_model_selection.md`

### GaMD check (CHEMBL848 / CHEMBL1129 / LxxLL, 3 × 100 ns each)

Pulled-back per-replica outputs (NOT trajectories — per-frame derived
metrics):

- `results/gamd_v2/<system>_rep<1-3>__frames.csv` — per-frame
  interface-superposed and global-frame ligand RMSD, interface drift,
  minimum ligand–receptor distance, contact count
- `results/gamd_v2/<system>_rep<1-3>__rescontact.csv` — per-frame residue
  contact bitmaps
- `results/gamd_v2/<system>_rep<1-3>__resocc.csv` — per-residue contact
  occupancies
- `results/gamd_v2/<system>_rep<1-3>__gamd-para.dat` / `__gamd-reweight.dat`
  — GaMD boost parameters and reweighting data as written by the fork
- `results/gamd_v2/remote_logs/<system>_rep<k>__gamd.in` / `__mdout.mdp` /
  `__analyse.log` — the exact run configurations and analysis logs
- Aggregates: `results/gamd_v2/gamd_summary.csv`,
  `results/gamd_v2/reweight_check.csv`, `results/gamd_v2/contact_occupancy.csv`,
  `results/gamd_v2/gamd_report.md`, `results/ess_scan_s025.csv`

### Also lost (and not archived)

- Batch-1 per-run outputs (non-bonded Zn²⁺ protocol): only the aggregated
  audit trail survives (`results/21_screen_report.md`,
  `results/21_screen_byligand.csv`)
- The 2 ns smoke-test run and the σ₀ = 6.0 kcal/mol failed-boost GaMD
  attempt: survive only as textual records in the manuscript and the
  boost-diagnostics discussion

## How to regenerate

1. **Rebuild the solvated system** from `templates/` with the archived tleap
   input `templates/build_system.leap.in` (receptor `rec_cym.pdb`, zinc
   ions `zn1_chainC.pdb`/`zn2_chainC.pdb`, reconstructed ZAFF parameters
   `templates/zaff_reconstructed.frcmod` — see its header for the
   reconstruction statement and acceptance basis).
2. **Re-run docking** with `docking/batch_docking_75.py` (Vina inputs
   already in `docking/4GIZ_poses/`) or reuse the archived best poses.
3. **Re-run the MD/MM-PBSA pipeline** with `code/12_pipeline_v3.py` /
   `code/16_pipeline_v3_amber.py` on a GPU host (see `SOFTWARE_VERSIONS.md`
   for the required GROMACS + AmberTools stack).
4. **Re-run GaMD** with `code/40_gamd_build.py` (build the GROMACS
   2025.4+GaMD v2.1 fork) and `code/41_gamd_run.py` (deployment, mdp
   templates and the exact `gamd.in` generation); analysis via
   `code/42_gamd_analyse.py` / `code/43_gamd_reweight.py`.

Note: exact regeneration of the *same numbers* is not guaranteed (GPU
non-determinism, unrecorded conda versions, and the batch-2 start structures
having been the docking-pose systems built on the remote host). The archived
per-frame CSVs, not a re-run, are the canonical record of what happened.
