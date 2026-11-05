#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""04_similarity_cluster_holdout.py

Phylogenetic-style generalization test for the E6 risk classifier.

Type-level grouping (script 03) still allows the model to be tested on a type whose
close relatives remain in training. This script groups HPV TYPES by pairwise E6
sequence similarity (average-linkage clustering, no external clade annotation) and
re-runs leave-one-cluster-out evaluation, i.e. the test types have no close
relative in the training set.

Caveat: cluster membership is data-derived, not a published clade assignment.
Formal Alpha/Beta genus names must still be verified against a curated reference
before publication (field "clade_name": "unverified").

Outputs (rebuild_20261025/results/):
  04_type_clusters.csv            type -> similarity cluster, n sequences, risk label
  04_cluster_holdout_metrics.csv  per-cluster held-out metrics
  04_cluster_summary.json         pooled metrics + cluster composition
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
from sklearn.base import clone
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, brier_score_loss, f1_score, precision_score,
    recall_score, roc_auc_score,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

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
N_CLUSTERS = 6


def representative_sequences(df):
    reps = {}
    for t, sub in df.groupby("hpv_type"):
        seq = max(sub["sequence"].astype(str), key=len)
        reps[str(t)] = seq
    return reps


def similarity_matrix(types, reps):
    n = len(types)
    sim = np.eye(n)
    for i, j in combinations(range(n), 2):
        r = SequenceMatcher(None, reps[types[i]], reps[types[j]]).ratio()
        sim[i, j] = sim[j, i] = r
    return sim


def assign_clusters(types, sim, k=N_CLUSTERS):
    dist = 1.0 - sim
    np.fill_diagonal(dist, 0.0)
    iu = np.triu_indices(len(types), k=1)
    Z = linkage(dist[iu], method="average")
    labels = fcluster(Z, t=k, criterion="maxclust")
    return {t: int(l) for t, l in zip(types, labels)}


def metrics(y_true, y_prob):
    y_pred = (np.asarray(y_prob) >= 0.5).astype(int)
    y_true = np.asarray(y_true)
    out = {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall": float(recall_score(y_true, y_pred, zero_division=0)),
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "brier": float(brier_score_loss(y_true, y_prob)),
    }
    out["auc"] = float(roc_auc_score(y_true, y_prob)) if len(np.unique(y_true)) == 2 else float("nan")
    return out


def evaluate_k(df, types, reps, k, est):
    """Leave-one-cluster-out for one clustering resolution k."""
    sim = similarity_matrix(types, reps)
    cluster_of = assign_clusters(types, sim, k=k)
    df = df.copy()
    df["cluster"] = df["type_str"].map(cluster_of)

    mixed = [c for c in df["cluster"].unique() if df.loc[df["cluster"] == c, "y"].nunique() > 1]
    rows, pooled_idx, pooled_prob = [], [], []
    for cl in sorted(df["cluster"].unique()):
        te = df.index[df["cluster"] == cl].to_numpy()
        tr = df.index[df["cluster"] != cl].to_numpy()
        if len(np.unique(df.loc[tr, "y"])) < 2:
            rows.append({"k": k, "cluster": int(cl), "n_test": len(te),
                         "note": "training single-class: fold not evaluable"})
            continue
        if len(np.unique(df.loc[te, "y"])) < 2:
            rows.append({"k": k, "cluster": int(cl), "n_test": len(te),
                         "n_test_types": int(df.loc[te, "type_str"].nunique()),
                         "types": "|".join(sorted(df.loc[te, "type_str"].unique())),
                         "note": "test single-class: AUC undefined"})
            continue
        np.random.seed(SEED)
        pipe = Pipeline([("scaler", StandardScaler()), ("clf", clone(est))])
        pipe.fit(df.loc[tr, FEATURES].values, df.loc[tr, "y"].values)
        p = pipe.predict_proba(df.loc[te, FEATURES].values)[:, 1]
        m = metrics(df.loc[te, "y"].values, p)
        rows.append({"k": k, "cluster": int(cl), "n_test": len(te),
                     "n_test_types": int(df.loc[te, "type_str"].nunique()),
                     "types": "|".join(sorted(df.loc[te, "type_str"].unique())), **m})
        pooled_idx.extend(te.tolist())
        pooled_prob.extend(p.tolist())

    pooled = metrics(df.loc[pooled_idx, "y"].values, np.asarray(pooled_prob)) if pooled_idx else None
    return rows, pooled, len(mixed), cluster_of


def main():
    df0 = pd.read_csv(DATA)
    for c in FEATURES:
        df0[c] = pd.to_numeric(df0[c], errors="coerce")
    df0 = df0.dropna(subset=FEATURES).reset_index(drop=True)
    df0["y"] = df0["risk_level"].astype(str).str.lower().str.startswith("high").astype(int)
    df0["type_str"] = df0["hpv_type"].astype(str)

    types = sorted(df0["type_str"].unique())
    reps = representative_sequences(df0)

    est = LogisticRegression(max_iter=5000, class_weight="balanced", random_state=SEED)
    all_rows, per_k, comp_frames = [], {}, []
    for k in (4, 6, 8, 10, 12):
        rows, pooled, n_mixed, cluster_of = evaluate_k(df0, types, reps, k, est)
        all_rows.extend(rows)
        n_evaluable_folds = sum(1 for r in rows
                                if isinstance(r.get("auc"), float) and not np.isnan(r["auc"]))
        per_k[k] = {
            "n_clusters": int(len(set(cluster_of.values()))),
            "n_mixed_clusters": int(n_mixed),
            "n_evaluable_folds": int(n_evaluable_folds),
            "evaluable": bool(n_evaluable_folds > 0),
            "pooled": pooled,
        }
        tmp = df0.copy()
        tmp["k"] = k
        tmp["cluster"] = tmp["type_str"].map(cluster_of)
        comp_frames.append(tmp.groupby(["k", "cluster", "type_str", "risk_level"])
                              .size().reset_index(name="n"))

    pd.DataFrame(all_rows).to_csv(OUT / "04_cluster_holdout_metrics.csv", index=False)
    comp = pd.concat(comp_frames, ignore_index=True)
    comp.to_csv(OUT / "04_type_clusters.csv", index=False)

    evaluable = [k for k, v in per_k.items() if v["evaluable"]]
    summary = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "n_rows": int(len(df0)),
        "n_types": int(len(types)),
        "clade_name": "unverified (data-derived similarity clusters, not published Alpha/Beta genera)",
        "model": "LogisticRegression(balanced) + StandardScaler inside Pipeline",
        "per_k": per_k,
        "evaluable_k": evaluable,
        "conclusion": (
            "Risk label is almost perfectly confounded with sequence-similarity cluster: at every "
            f"tested resolution only one cluster mixes high- and low-risk types, so most "
            "leave-one-cluster-out folds are degenerate (single-class test, or removing the mixed "
            f"cluster leaves training single-class). Evaluable folds exist only at k=10 (1 fold, "
            "held-out 16|30|31|33|35|51|52|56|58|73|82 -> AUC 0.955) and k=12 (1 fold, "
            "30|51|56|82 -> AUC 1.000). Two informative folds are not enough to claim phylogenetic "
            "generalization; the type-level model is effectively a clade-membership proxy and must "
            "be reframed as a descriptive clade-level physicochemical signature, or re-labelled "
            "with variant-level phenotypes."
        ),
        "seed": SEED,
    }
    with open(OUT / "04_cluster_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
