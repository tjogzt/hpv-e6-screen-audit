#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""08_clade_signature_descriptive.py

Descriptive re-frame of the E6 ML line: E6 physicochemical properties as
signatures of HPV similarity clusters (data-derived, k=10), NOT risk prediction.

Analyses (all descriptive / inferential statistics with multiple-testing control):
  1. Cluster composition table (types, n sequences, risk label consistency).
  2. Per-cluster feature profiles (median, IQR) for the 12 physicochemical features.
  3. Kruskal-Wallis test per feature across clusters; BH-FDR correction;
     epsilon-squared effect size — computed at TWO units of analysis:
     sequence level (n = 524; historical) and type level (per-type medians,
     n = 45; pseudoreplication-free, revision R4). Conclusions follow the
     type-level results.
  4. One-vs-rest logistic AUC per cluster with 5-fold CV, scaler in pipeline,
     reported under TWO fold constructions (revision R5):
     ungrouped StratifiedKFold (historical) and StratifiedGroupKFold with
     groups = HPV type (no type shared between train and test). Grouped CV
     is undefined for single-type clusters (reported as NaN by design).
  5. Figure: Fig1-style signature heatmap (row-scaled median profiles) + AUC barplot,
     exported as PDF with 7pt fonts, English labels only (project rules).

Outputs (rebuild_20261025/results/):
  08_cluster_composition.csv
  08_cluster_feature_profiles.csv
  08_feature_kruskal.csv            sequence level (n = 524)
  08_feature_kruskal_typelevel.csv  type level (per-type medians, n = 45)
  08_cluster_auc.csv                both fold constructions side by side
  Fig3A_clade_signature_heatmap.pdf
  Fig3B_cluster_discriminability.pdf
  08_clade_summary.json
"""

from __future__ import annotations

import json
from datetime import datetime
from difflib import SequenceMatcher
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import (StratifiedGroupKFold, StratifiedKFold,
                                     cross_val_predict)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
REBUILD = Path(__file__).resolve().parents[1]
DATA = REBUILD / "data" / "02_dataset_clean.csv"
OUT = REBUILD / "results"
OUT.mkdir(parents=True, exist_ok=True)

FEATURES = [
    "length", "molecular_weight", "aromaticity", "instability_index",
    "helix_propensity", "turn_propensity", "sheet_propensity",
    "hydrophobicity", "isoelectric_point", "net_charge",
    "positive_charged", "negative_charged",
]
SEED = 42
N_CLUSTERS = 10  # resolution with the most granular, still interpretable clusters

plt.rcParams.update({
    "font.size": 7,
    "font.family": "DejaVu Sans",
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi": 300,
})


def representative_sequences(df):
    return {t: max(sub["sequence"].astype(str), key=len)
            for t, sub in df.groupby("type_str")}


def assign_clusters(types, reps, k):
    n = len(types)
    sim = np.eye(n)
    for i, j in combinations(range(n), 2):
        r = SequenceMatcher(None, reps[types[i]], reps[types[j]]).ratio()
        sim[i, j] = sim[j, i] = r
    dist = 1.0 - sim
    np.fill_diagonal(dist, 0.0)
    Z = linkage(dist[np.triu_indices(n, k=1)], method="average")
    labels = fcluster(Z, t=k, criterion="maxclust")
    return {t: int(l) for t, l in zip(types, labels)}


def epsilon_squared_kw(h_stat, n, k):
    return max(0.0, (h_stat - k + 1) / (n - k))


def one_vs_rest_auc(X, y_bin, seed):
    if len(np.unique(y_bin)) < 2 or min(np.bincount(y_bin)) < 5:
        return float("nan")
    pipe = Pipeline([("scaler", StandardScaler()),
                     ("clf", LogisticRegression(max_iter=5000, class_weight="balanced",
                                                random_state=seed))])
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    try:
        prob = cross_val_predict(pipe, X, y_bin, cv=cv, method="predict_proba")[:, 1]
        return float(roc_auc_score(y_bin, prob))
    except ValueError:
        return float("nan")


def one_vs_rest_auc_grouped(X, y_bin, groups, seed):
    """Type-grouped CV: no HPV type appears in both train and test folds.

    n_splits = min(5, n_types_in_cluster); undefined (NaN) for single-type
    clusters, where any type-honest evaluation is impossible by construction.
    """
    pos_types = np.unique(groups[y_bin == 1])
    n_types = len(pos_types)
    if n_types < 2 or min(np.bincount(y_bin)) < 2:
        return float("nan"), 0
    n_splits = min(5, n_types)
    pipe = Pipeline([("scaler", StandardScaler()),
                     ("clf", LogisticRegression(max_iter=5000, class_weight="balanced",
                                                random_state=seed))])
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    try:
        prob = cross_val_predict(pipe, X, y_bin, groups=groups, cv=cv,
                                 method="predict_proba")[:, 1]
        return float(roc_auc_score(y_bin, prob)), n_splits
    except ValueError:
        return float("nan"), n_splits


def fig_heatmap(profiles, path):
    """Row-scaled median feature profiles per cluster (features x clusters)."""
    mat = profiles.set_index("cluster")[FEATURES].T
    z = (mat.subtract(mat.mean(axis=1), axis=0)
            .divide(mat.std(axis=1).replace(0, np.nan), axis=0))
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(z.columns)))
    ax.set_xticklabels([f"C{c}" for c in z.columns])
    ax.set_yticks(range(len(z.index)))
    ax.set_yticklabels(z.index)
    ax.set_xlabel("Similarity cluster")
    ax.set_ylabel("Feature (row z-score)")
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cb.ax.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_auc_bar(auc_df, path):
    """Paired bars: ungrouped vs type-grouped CV (revision R5).

    Grouped CV is undefined for single-type clusters (NaN); those slots are
    left empty and annotated 'n/a (1 type)'.
    """
    d = auc_df.sort_values("one_vs_rest_auc", ascending=True).reset_index(drop=True)
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    ax.barh(y + 0.2, d["one_vs_rest_auc"], color="#9dc3dd", height=0.38,
            label="ungrouped 5-fold CV")
    ax.barh(y - 0.2, d["one_vs_rest_auc_grouped"], color="#215f85", height=0.38,
            label="type-grouped CV")
    for yi, (_, r) in zip(y, d.iterrows()):
        if pd.isna(r["one_vs_rest_auc_grouped"]):
            ax.text(0.02, yi - 0.2, "n/a (1 type)", fontsize=7, va="center",
                    color="#555555")
    ax.axvline(0.5, color="#888888", lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"C{c}" for c in d["cluster"]])
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("One-vs-rest AUC")
    ax.set_ylabel("Similarity cluster")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2,
              frameon=False, borderaxespad=0.0)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    df = pd.read_csv(DATA)
    for c in FEATURES:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=FEATURES).reset_index(drop=True)
    df["type_str"] = df["hpv_type"].astype(str)
    df["y_risk"] = df["risk_level"].astype(str).str.lower().str.startswith("high").astype(int)

    types = sorted(df["type_str"].unique())
    reps = representative_sequences(df)
    cluster_of = assign_clusters(types, reps, N_CLUSTERS)
    df["cluster"] = df["type_str"].map(cluster_of)

    # 1. composition
    comp = (df.groupby(["cluster", "type_str"])
              .agg(n=("sequence_id", "size"), risk=("risk_level", "first"))
              .reset_index())
    risk_mix = df.groupby("cluster")["y_risk"].nunique()
    comp["cluster_risk_homogeneous"] = comp["cluster"].map(risk_mix) == 1
    comp.to_csv(OUT / "08_cluster_composition.csv", index=False)

    # 2. profiles
    profiles = (df.groupby("cluster")[FEATURES]
                  .median()
                  .reset_index())
    profiles.to_csv(OUT / "08_cluster_feature_profiles.csv", index=False)

    # 3. Kruskal-Wallis + BH-FDR + effect size — sequence level (n = 524)
    rows = []
    n = len(df)
    k = df["cluster"].nunique()
    for feat in FEATURES:
        groups = [g[feat].values for _, g in df.groupby("cluster")]
        h, p = kruskal(*groups)
        rows.append({"feature": feat, "H": float(h), "p_raw": float(p),
                     "epsilon2": float(epsilon_squared_kw(h, n, k))})
    kw = pd.DataFrame(rows)
    kw["p_bh_fdr"] = multipletests(kw["p_raw"], method="fdr_bh")[1]
    kw = kw.sort_values("epsilon2", ascending=False).reset_index(drop=True)
    kw.to_csv(OUT / "08_feature_kruskal.csv", index=False)

    # 3b. Same test at type level (per-type medians, n = 45) — the unit the
    #     manuscript declares in §2.5c; free of within-type pseudoreplication
    #     (revision R4). Conclusions in the text follow this table.
    dft = (df.groupby(["cluster", "type_str"])[FEATURES]
             .median()
             .reset_index())
    rows_t = []
    nt = len(dft)
    for feat in FEATURES:
        groups_t = [g[feat].values for _, g in dft.groupby("cluster")]
        h, p = kruskal(*groups_t)
        rows_t.append({"feature": feat, "H": float(h), "p_raw": float(p),
                       "epsilon2": float(epsilon_squared_kw(h, nt, k))})
    kwt = pd.DataFrame(rows_t)
    kwt["p_bh_fdr"] = multipletests(kwt["p_raw"], method="fdr_bh")[1]
    kwt = kwt.sort_values("epsilon2", ascending=False).reset_index(drop=True)
    kwt.to_csv(OUT / "08_feature_kruskal_typelevel.csv", index=False)

    # 4. one-vs-rest AUC per cluster — ungrouped and type-grouped CV (revision R5)
    X = df[FEATURES].values
    groups_all = df["type_str"].values
    auc_rows = []
    for cl in sorted(df["cluster"].unique()):
        yb = (df["cluster"] == cl).astype(int).values
        np.random.seed(SEED)
        auc = one_vs_rest_auc(X, yb, SEED)
        np.random.seed(SEED)
        auc_g, n_splits_g = one_vs_rest_auc_grouped(X, yb, groups_all, SEED)
        auc_rows.append({
            "cluster": int(cl),
            "n": int(yb.sum()),
            "n_types": int(df.loc[df["cluster"] == cl, "type_str"].nunique()),
            "risk_homogeneous": bool(risk_mix[cl] == 1),
            "one_vs_rest_auc": auc,
            "one_vs_rest_auc_grouped": auc_g,
            "grouped_cv_splits": n_splits_g,
        })
    auc_df = pd.DataFrame(auc_rows)
    auc_df.to_csv(OUT / "08_cluster_auc.csv", index=False)

    # 5. figures
    fig_heatmap(profiles, OUT / "Fig3A_clade_signature_heatmap.pdf")
    fig_auc_bar(auc_df, OUT / "Fig3B_cluster_discriminability.pdf")

    summary = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "frame": "Descriptive: physicochemical signatures of data-derived similarity clusters "
                 "(NOT risk prediction; risk labels are confounded with clusters)",
        "n_rows": int(len(df)),
        "n_types": int(len(types)),
        "n_clusters": int(k),
        "clusters_risk_homogeneous": int((risk_mix == 1).sum()),
        "top_features_by_epsilon2": kw.head(5)[["feature", "epsilon2", "p_bh_fdr"]].to_dict(orient="records"),
        "top_features_by_epsilon2_typelevel": kwt.head(5)[["feature", "epsilon2", "p_bh_fdr"]].to_dict(orient="records"),
        "cluster_auc": auc_df.to_dict(orient="records"),
        "figures": ["Fig3A_clade_signature_heatmap.pdf", "Fig3B_cluster_discriminability.pdf"],
        "seed": SEED,
        "revision_notes": {
            "R4": "KW/epsilon2 now reported at sequence level (n=524) and type level "
                  "(per-type medians, n=45); conclusions follow the type-level table.",
            "R5": "One-vs-rest AUC reported under ungrouped StratifiedKFold and "
                  "StratifiedGroupKFold (groups=HPV type, n_splits=min(5, n_types)); "
                  "grouped CV is undefined for single-type clusters (NaN by design).",
        },
    }
    with open(OUT / "08_clade_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(json.dumps({k2: summary[k2] for k2 in
                      ["n_rows", "n_types", "n_clusters", "clusters_risk_homogeneous",
                       "top_features_by_epsilon2"]}, ensure_ascii=False, indent=2))
    print(auc_df.to_string(index=False))


if __name__ == "__main__":
    main()
