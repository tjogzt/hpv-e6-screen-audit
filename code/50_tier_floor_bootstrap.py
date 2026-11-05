#!/usr/bin/env python3
"""
50_tier_floor_bootstrap.py — S4 (Round-1 revision): sampling uncertainty of the
control-calibrated tier thresholds.

The tier floor was defined as min of the n=3 LxxLL positive-control retention
fractions and the ceiling as max of the n=3 DMSO fractions. Both statistics are
extremes of tiny samples; this script quantifies their instability by
nonparametric bootstrap over the three replica values (resample with
replacement, seed 42, 100 000 resamples) and reports percentile 95% CIs.

Inputs (primary, no hand transcription):
  results/36_screen_full.csv  (batch-2 run table; controls = LxxLL_peptide / DMSO_pose1)

Output:
  results/50_tier_floor_bootstrap.json

Run: /opt/anaconda3/bin/python3 code/50_tier_floor_bootstrap.py
"""
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results" / "36_screen_full.csv"
OUT = ROOT / "results" / "50_tier_floor_bootstrap.json"

N_BOOT = 100_000
SEED = 42


def main() -> None:
    df = pd.read_csv(SRC)
    lxx = df[df.CHEMBL == "LxxLL_peptide"]["kept_frac"].to_numpy(float)
    dmso = df[df.CHEMBL == "DMSO_pose1"]["kept_frac"].to_numpy(float)
    assert len(lxx) == 3 and len(dmso) == 3, "expected 3 replicas per control"

    rng = np.random.default_rng(SEED)
    boot_floor = rng.choice(lxx, size=(N_BOOT, 3), replace=True).min(axis=1)
    boot_ceil = rng.choice(dmso, size=(N_BOOT, 3), replace=True).max(axis=1)

    rec = {
        "source": "results/36_screen_full.csv (batch 2)",
        "lxxll_replicas": [round(x, 6) for x in lxx],
        "dmso_replicas": [round(x, 6) for x in dmso],
        "floor_point": round(float(lxx.min()), 6),
        "ceiling_point": round(float(dmso.max()), 6),
        "bootstrap": {
            "n_resamples": N_BOOT,
            "seed": SEED,
            "floor_ci95": [round(float(np.percentile(boot_floor, 2.5)), 6),
                           round(float(np.percentile(boot_floor, 97.5)), 6)],
            "ceiling_ci95": [round(float(np.percentile(boot_ceil, 2.5)), 6),
                             round(float(np.percentile(boot_ceil, 97.5)), 6)],
            "p_floor_equals_lowest_replica": round(float((boot_floor == lxx.min()).mean()), 4),
            "p_floor_at_or_above_0.5": round(float((boot_floor >= 0.5).mean()), 4),
        },
        "interpretation": (
            "With n = 3 the floor is the minimum of three draws; the bootstrap "
            "distribution is discrete and wide. The 95% CI spans up to the "
            "second-lowest replica, and there is a non-trivial probability that "
            "a re-drawn control trio would set the floor at or above 0.5 — i.e. "
            "candidate-vs-floor comparisons near the boundary are not stable, "
            "and per-replica differences within a ligand are not ranking evidence."
        ),
        "generated": datetime.now().isoformat(timespec="seconds"),
    }
    OUT.write_text(json.dumps(rec, indent=2))
    print(json.dumps(rec["bootstrap"], indent=2))
    print(f"[50] wrote {OUT}")


if __name__ == "__main__":
    main()
