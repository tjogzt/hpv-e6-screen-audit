#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""03_grouped_nested_validation.py

Leak-free validation of the E6 type-level risk classifier on the cleaned dataset
produced by 02_clean_ml_taxonomy.py.

What it fixes relative to scripts/model_training/strict_model_validation.py:
  * groups = HPV type  -> no same-type samples across the train/test boundary
  * StandardScaler lives inside a Pipeline -> fitted on training folds only
  * model selection happens in an INNER loop only -> reported score is outer
  * feature importance = held-out permutation importance on outer test folds

It also re-runs the historical (leaky) protocol on the same data so the drop in
performance is measurable instead of asserted.

Outputs (rebuild_20261025/results/):
  03_validation_leaky_vs_grouped.csv   leaky baseline vs grouped protocol
  03_nested_grouped_metrics.csv        per-family, per-outer-fold metrics
  03_permutation_importance.csv        held-out permutation importance
  03_validation_summary.json           machine-readable summary
"""

from __future__ import annotations

import csv
import json
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, brier_score_loss, confusion_matrix, f1_score,
    precision_score, recall_score, roc_auc_score,
)
from sklearn.model_selection import (
    GridSearchCV, StratifiedGroupKFold, StratifiedKFold, cross_val_predict,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path(__file__).resolve().parents[2]
REBUILD = Path(__file__).resolve().parents[1]
DATA = REBUILD / "data" / "02_dataset_clean.csv"
OUT = REBUILD / "results"
OUT.mkdir(parents=True, exist_ok=True)
OUT_DATA = REBUILD / "data"

FEATURES = [
    "length", "molecular_weight", "aromaticity", "instability_index",
    "helix_propensity", "turn_propensity", "sheet_propensity",
    "hydrophobicity", "isoelectric_point", "net_charge",
    "positive_charged", "negative_charged",
]

SEED = 42


def metrics(y_true, y_prob, y_pred=None):
    if y_pred is None:
        y_pred = (np.asarray(y_prob) >= 0.5).astype(int)
    y_true = np.asarray(y_true)
    out = {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall": float(recall_score(y_true, y_pred, zero_division=0)),
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "brier": float(brier_score_loss(y_true, y_prob)),
    }
    if len(np.unique(y_true)) == 2:
        out["auc"] = float(roc_auc_score(y_true, y_prob))
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        out.update({"tp": int(tp), "tn": int(tn), "fp": int(fp), "fn": int(fn)})
    else:
        out["auc"] = float("nan")
    return out


def load():
    df = pd.read_csv(DATA)
    for c in FEATURES:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=FEATURES).reset_index(drop=True)
    y = (df["risk_level"].astype(str).str.lower().str.startswith("high")).astype(int).values
    X = df[FEATURES].values
    groups = df["hpv_type"].astype(str).values
    return df, X, y, groups


def leaky_baseline(X, y):
    """Historical protocol: global scaler + StratifiedKFold (no groups)."""
    np.random.seed(SEED)
    model = SVC(kernel="rbf", C=1.0, probability=True, class_weight="balanced", random_state=SEED)
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)  # leak: fitted on all rows
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    prob = cross_val_predict(model, Xs, y, cv=cv, method="predict_proba")[:, 1]
    return metrics(y, prob)


def grouped_cv(X, y, groups, estimator, n_splits=5):
    """Grouped CV: scaler inside pipeline, type never crosses the fold boundary."""
    np.random.seed(SEED)
    pipe = Pipeline([("scaler", StandardScaler()), ("clf", clone(estimator))])
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    prob = cross_val_predict(pipe, X, y, groups=groups, cv=cv, method="predict_proba")[:, 1]
    return metrics(y, prob), prob


def nested_grouped(X, y, groups, name, estimator, grid, n_outer=5, n_inner=3):
    """Nested CV: hyper-parameters chosen inside the outer training fold only."""
    outer = StratifiedGroupKFold(n_splits=n_outer, shuffle=True, random_state=SEED)
    inner = StratifiedGroupKFold(n_splits=n_inner, shuffle=True, random_state=SEED)
    rows, probs, fold_models = [], [], []
    for k, (tr, te) in enumerate(outer.split(X, y, groups)):
        if len(np.unique(y[te])) < 2:
            rows.append({"model": name, "fold": k, "auc": float("nan"),
                         "n_test": int(len(te)), "best_params": "single-class fold"})
            probs.append((te, np.full(len(te), np.nan)))
            continue
        np.random.seed(SEED + k)
        search = GridSearchCV(
            Pipeline([("scaler", StandardScaler()), ("clf", clone(estimator))]),
            param_grid={f"clf__{kk}": vv for kk, vv in grid.items()},
            scoring="roc_auc", cv=inner, n_jobs=-1, refit=True,
        )
        search.fit(X[tr], y[tr], groups=groups[tr])
        best = search.best_estimator_
        p = best.predict_proba(X[te])[:, 1]
        m = metrics(y[te], p)
        rows.append({"model": name, "fold": k, "best_params": json.dumps(search.best_params_),
                     "n_test": int(len(te)), **m})
        probs.append((te, p))
        fold_models.append((te, best))
    pooled_idx = np.concatenate([i for i, _ in probs])
    pooled_p = np.concatenate([p for _, p in probs])
    order = np.argsort(pooled_idx)
    pooled = metrics(y[pooled_idx[order]], pooled_p[order])
    return rows, pooled, fold_models


def heldout_permutation(X, y, fold_models, n_repeats=10):
    """Permutation importance computed on held-out outer test folds."""
    imp = []
    for te, model in fold_models:
        np.random.seed(SEED)
        r = permutation_importance(model, X[te], y[te], scoring="roc_auc",
                                   n_repeats=n_repeats, random_state=SEED, n_jobs=-1)
        imp.append(r.importances_mean)
    arr = np.vstack(imp)
    return pd.DataFrame({
        "feature": FEATURES,
        "importance_mean": arr.mean(axis=0),
        "importance_std": arr.std(axis=0),
    }).sort_values("importance_mean", ascending=False).reset_index(drop=True)


def main():
    df, X, y, groups = load()
    print(f"cleaned dataset: {X.shape[0]} rows | {X.shape[1]} features | "
          f"{len(np.unique(groups))} types | classes {dict(zip(*np.unique(y, return_counts=True)))}")

    base_rows = []
    leaky = leaky_baseline(X, y)
    base_rows.append({"protocol": "leaky_StratifiedKFold_global_scaler", "model": "SVC", **leaky})

    estimators = {
        "LogisticRegression": (
            LogisticRegression(max_iter=5000, class_weight="balanced", random_state=SEED),
            {"C": [0.1, 1.0, 10.0]},
        ),
        "SVC": (
            SVC(kernel="rbf", probability=True, class_weight="balanced", random_state=SEED),
            {"C": [0.5, 1.0, 5.0], "gamma": ["scale", 0.1, 0.01]},
        ),
        "RandomForest": (
            RandomForestClassifier(n_estimators=500, class_weight="balanced",
                                   random_state=SEED, n_jobs=-1),
            {"max_depth": [5, 10, None], "min_samples_leaf": [1, 3]},
        ),
    }

    for name, (est, _grid) in estimators.items():
        m, _ = grouped_cv(X, y, groups, est)
        base_rows.append({"protocol": "grouped_StratifiedGroupKFold_pipeline_scaler",
                          "model": name, **m})

    pd.DataFrame(base_rows).to_csv(OUT / "03_validation_leaky_vs_grouped.csv", index=False)

    nested_rows, pooled_by_model, folds_by_model = [], {}, {}
    for name, (est, grid) in estimators.items():
        rows, pooled, folds = nested_grouped(X, y, groups, name, est, grid)
        nested_rows.extend(rows)
        pooled_by_model[name] = pooled
        folds_by_model[name] = folds
        print(f"nested {name}: AUC={pooled['auc']:.3f} ACC={pooled['accuracy']:.3f}")

    pd.DataFrame(nested_rows).to_csv(OUT / "03_nested_grouped_metrics.csv", index=False)

    best_name = max(pooled_by_model, key=lambda k: pooled_by_model[k]["auc"])
    imp = heldout_permutation(X, y, folds_by_model[best_name])
    imp.insert(0, "model", best_name)
    imp.to_csv(OUT / "03_permutation_importance.csv", index=False)

    summary = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "input": str(DATA.relative_to(ROOT)),
        "n_rows": int(X.shape[0]),
        "n_features": int(X.shape[1]),
        "n_types": int(len(np.unique(groups))),
        "class_counts": {str(k): int(v) for k, v in zip(*np.unique(y, return_counts=True))},
        "leaky_baseline": leaky,
        "pooled_nested": {k: {kk: (None if isinstance(vv, float) and np.isnan(vv) else vv)
                              for kk, vv in v.items()} for k, v in pooled_by_model.items()},
        "best_family_by_auc": best_name,
        "permutation_importance_top": imp.head(5).to_dict(orient="records"),
        "seed": SEED,
        "caution": "Metrics describe generalization to unseen HPV TYPES on type-level proxy labels, not clinical variant pathogenicity.",
    }
    with open(OUT / "03_validation_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(json.dumps({k: summary[k] for k in
                      ["leaky_baseline", "pooled_nested", "best_family_by_auc"]},
                     ensure_ascii=False, indent=2))
    print(imp.head(6).to_string(index=False))


if __name__ == "__main__":
    main()
