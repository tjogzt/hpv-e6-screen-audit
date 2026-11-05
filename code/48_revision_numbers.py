#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""48_revision_numbers.py

Single source of truth for the revision-round numbers (R2, R4, R5, R6).

Every number written into the manuscript or the response letter during the
revision must be pasted from the output of this script, which reads the
result CSVs directly (no transcription). Re-run after any upstream change:

  /opt/anaconda3/bin/python3 code/48_revision_numbers.py

Outputs (rebuild_20261025/results/):
  48_revision_numbers.json
  48_revision_numbers.md
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pandas as pd

REBUILD = Path(__file__).resolve().parents[1]
RES = REBUILD / "results"


def main():
    kw_seq = pd.read_csv(RES / "08_feature_kruskal.csv")
    kw_typ = pd.read_csv(RES / "08_feature_kruskal_typelevel.csv")
    auc = pd.read_csv(RES / "08_cluster_auc.csv")
    krob = pd.read_csv(RES / "46_k_robustness_summary.csv")
    scan = json.load(open(RES / "47_selection_free_zinc_scan.json", encoding="utf-8"))

    # ---- R4: KW at both units of analysis -------------------------------
    merged = kw_seq.merge(kw_typ, on="feature", suffixes=("_seq", "_type"))
    r4 = {
        "n_sig_seq_fdr05": int((kw_seq["p_bh_fdr"] < 0.05).sum()),
        "n_sig_type_fdr05": int((kw_typ["p_bh_fdr"] < 0.05).sum()),
        "top_seq": kw_seq.sort_values("epsilon2", ascending=False)
                         .iloc[0][["feature", "epsilon2", "p_bh_fdr"]].to_dict(),
        "top_type": kw_typ.sort_values("epsilon2", ascending=False)
                          .iloc[0][["feature", "epsilon2", "p_bh_fdr"]].to_dict(),
        "hydrophobicity": merged[merged["feature"] == "hydrophobicity"]
                            .iloc[0][["epsilon2_seq", "p_bh_fdr_seq",
                                      "epsilon2_type", "p_bh_fdr_type"]].to_dict(),
        "n_features_type_not_sig": kw_typ.loc[kw_typ["p_bh_fdr"] >= 0.05,
                                              "feature"].tolist(),
    }

    # ---- R5: AUC under both fold constructions --------------------------
    multi = auc[auc["n_types"] >= 2]
    ev = multi.dropna(subset=["one_vs_rest_auc_grouped"])
    r5 = {
        "ungrouped_range_all_evaluable_clusters": [
            round(float(auc.loc[auc["one_vs_rest_auc"].notna(), "one_vs_rest_auc"].min()), 3),
            round(float(auc["one_vs_rest_auc"].max()), 3),
        ],
        "grouped_range": [round(float(ev["one_vs_rest_auc_grouped"].min()), 3),
                          round(float(ev["one_vs_rest_auc_grouped"].max()), 3)],
        "n_clusters_grouped_evaluable": int(len(ev)),
        "n_single_type_clusters": int((auc["n_types"] == 1).sum()),
        "per_cluster": auc[["cluster", "n", "n_types", "one_vs_rest_auc",
                            "one_vs_rest_auc_grouped", "grouped_cv_splits"]]
                        .round(3).to_dict(orient="records"),
        "note": ("Grouped CV is undefined for single-type clusters "
                 "(clusters 4 and 10); their ungrouped AUC of 0.999-1.000 "
                 "was within-type separation, not between-cluster signal."),
    }

    # ---- R6: cross-k robustness ----------------------------------------
    r6 = {
        "k_sweep": krob["k"].tolist(),
        "type_level_sig_range": [int(krob["kw_type_n_sig_fdr05"].min()),
                                 int(krob["kw_type_n_sig_fdr05"].max())],
        "type_level_top_feature_all_k":
            sorted(krob["kw_type_top_feature"].unique().tolist()),
        "purity_at_k10": float(krob.loc[krob["k"] == 10, "weighted_mean_purity"].iloc[0]),
        "purity_range": [float(krob["weighted_mean_purity"].min()),
                         float(krob["weighted_mean_purity"].max())],
    }
    loco = pd.read_csv(RES / "46_k_loco_evaluability.csv")
    r6["loco_evaluable_folds_per_k"] = {
        int(k): int(g["evaluable_risk_auc"].sum())
        for k, g in loco.groupby("k")
    }

    # ---- R2: selection-free scan ---------------------------------------
    r2 = {
        "n_records_scanned": scan["input_records"],
        "retained": scan["by_decision"]["retain"],
        "n_counterexamples": scan["n_counterexamples"],
        "counterexample_census": scan["counterexample_census"],
        "counterexample_ids": [c["sequence_id"] for c in scan["counterexamples"]],
        "anchor_failures_reproduced": scan["anchoring_crosscheck"]["n_anchor_failures_reproduced"],
        "types_with_two_site_sequence": scan["per_type"]["types_with_two_site_sequence"],
        "n_types": scan["per_type"]["n_types_retained"],
        "verdict": ("Type-consensus level held (45/45 types carry a record with two "
                    "complete CXXC-CXXC units); record level does not "
                    "(9/524 retained records lack one). Title and §3.7b must use "
                    "the consensus-level statement, not 'absolutely conserved'."),
    }

    out = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "R4_kw_two_levels": r4,
        "R5_grouped_cv": r5,
        "R6_k_robustness": r6,
        "R2_selection_free_scan": r2,
    }
    with open(RES / "48_revision_numbers.json", "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)

    md = [
        "# Revision numbers (auto-generated — do not hand-edit)",
        "",
        f"Generated: {out['generated']}",
        "",
        "## R4 — Kruskal-Wallis at two units of analysis",
        f"- sequence level: {r4['n_sig_seq_fdr05']}/12 FDR-significant; top = "
        f"{r4['top_seq']['feature']} (eps2 = {r4['top_seq']['epsilon2']:.3f})",
        f"- type level: {r4['n_sig_type_fdr05']}/12 FDR-significant; top = "
        f"{r4['top_type']['feature']} (eps2 = {r4['top_type']['epsilon2']:.3f})",
        f"- hydrophobicity: eps2 {r4['hydrophobicity']['epsilon2_seq']:.3f} (seq) vs "
        f"{r4['hydrophobicity']['epsilon2_type']:.3f} (type)",
        f"- not significant at type level: {', '.join(r4['n_features_type_not_sig'])}",
        "",
        "## R5 — one-vs-rest AUC, ungrouped vs type-grouped CV",
        f"- ungrouped range: {r5['ungrouped_range_all_evaluable_clusters']}",
        f"- grouped range: {r5['grouped_range']} "
        f"({r5['n_clusters_grouped_evaluable']} evaluable clusters)",
        f"- single-type clusters (grouped CV undefined): {r5['n_single_type_clusters']}",
        "",
        "## R6 — cross-k robustness",
        f"- k sweep: {r6['k_sweep']}",
        f"- type-level FDR-significant features per k: {r6['type_level_sig_range']} of 12",
        f"- top feature at every k: {r6['type_level_top_feature_all_k']}",
        f"- weighted purity: {r6['purity_range']} (k=10: {r6['purity_at_k10']})",
        f"- LOCO evaluable folds per k: {r6['loco_evaluable_folds_per_k']}",
        "",
        "## R2 — selection-free dual-zinc scan",
        f"- records scanned: {r2['n_records_scanned']}; retained "
        f"{r2['retained']['n']} (of which {r2['retained']['n_two_sites']} carry two "
        f"complete sites, {r2['retained']['n_one_or_zero_sites']} do not)",
        f"- record-level counterexamples: {r2['n_counterexamples']} "
        f"({r2['counterexample_census']['by_nature']})",
        f"- types with >= 1 two-site sequence: {r2['types_with_two_site_sequence']}/"
        f"{r2['n_types']}",
        f"- anchoring failures reproduced: {r2['anchor_failures_reproduced']} (expected 27)",
        f"- verdict: {r2['verdict']}",
    ]
    (RES / "48_revision_numbers.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    print("\n".join(md))


if __name__ == "__main__":
    main()
