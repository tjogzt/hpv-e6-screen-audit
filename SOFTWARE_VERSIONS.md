# Software versions

This file records the software stack used to produce the archived results,
with an honest "not recorded" wherever the version was not logged at run time.

## Local analysis environment (macOS, arm64)

Interpreter: `/opt/anaconda3/bin/python3.12` (Python 3.12.7)

| Package | Version | Used by |
|---|---|---|
| rdkit | 2026.03.3 | library building, property computation, fingerprints (`docking/02_fetch_ligand_library.py`, Part A) |
| scikit-learn | 1.6.1 | grouped/nested CV, permutation importance (`code/03_*.py`) |
| scipy | 1.13.1 | statistics (KW tests, bootstrap) |
| numpy | 2.2.6 | throughout |
| pandas | 2.3.3 | throughout |
| requests | (not recorded) | ChEMBL web API fetch |

R: 4.5.0 (2025-04-11) — figures (`code/32_*.R`, `34_*.R`, `35_*.R`, `44_*.R`).
Required R packages: ggplot2, patchwork (versions not recorded).

Statistical software for Part A type-level tests: scipy (Python), as above.

## Remote compute environment (MD / MM-PBSA / GaMD)

The MD pipeline ran on rented GPU instances accessed over SSH; the host
details below were recorded in the deployment scripts.

- **GROMACS 2025.4 + GaMD v2.1 fork** — built by `code/40_gamd_build.py`
  into `/usr/local/gromacs-gamd` on the remote host. The stock GROMACS
  2025.4 (`/opt/gmx2025.4`) does **not** contain GaMD. Fork constraints
  (documented in `code/41_gamd_run.py`): `integrator = md`, `gamd = yes`,
  `nstcalcenergy = 1`, C-rescale NPT only, GPU offload for non-bonded/PME/
  bonded/update.
- **GaMD replicas:** one replica per RTX A4000 GPU (4 GPUs, round-robin;
  9 jobs). Dual-boost (`igamd = 3`), σ₀ = 0.25 kcal/mol for both potential
  and dihedral boost.
- **AmberTools (MMPBSA.py) conda stack on the remote host:** Python 3.10
  with `ambertools`, `openbabel`, `parmed`, `mdtraj` from conda-forge.
  **Exact versions were not pinned at install time and are not recorded.**
- **AutoDock Vina:** used for the docking screen
  (`docking/batch_docking_75.py`). **Version not recorded.**
- **Batch-1/batch-2 conventional MD (20 ns × 3 replicas per complex):**
  GROMACS on the remote GPU instance (same host family as above; exact
  build for the cMD stage not separately recorded — see `TRAJECTORIES.md`).

## Data-source software

- **ChEMBL web API** (`https://www.ebi.ac.uk/chembl/api/`): accessed
  2025-10-02; **ChEMBL release identifier not recorded at access time.**
  See `DATA_PROVENANCE.md` for the query strings.

## Web resources

- UniProt and NCBI Protein access dates and record counts:
  see `DATA_PROVENANCE.md`.

## Notes

- The trajectories and the remote instance are gone (see
  `TRAJECTORIES.md`); the versions above are reconstructed from the
  archived deployment scripts and local logs, not from a surviving
  environment. Anything not recoverable from those artefacts is marked
  "not recorded" rather than guessed.
