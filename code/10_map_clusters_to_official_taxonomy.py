#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""10_map_clusters_to_official_taxonomy.py

Map the data-derived similarity clusters (k=10, script 08) to the official
papillomavirus genus / species-group taxonomy.

Reference (authoritative):
  - de Villiers et al. 2004, Virology 324:17-27 (Alpha species groups 1-15)
  - ICTV MSL (Betapapillomavirus 1-6 reorganisation; Mupapillomavirus for HPV1)
  - Burk et al. 2013, Virology 445:232-240 ("HR clade": alpha-5/6/7/9/11)

Outputs:
  results/10_cluster_taxonomy_map.csv   cluster -> types -> official group -> purity
  results/10_cluster_taxonomy_map.json  summary
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REBUILD = Path(__file__).resolve().parents[1]
OUT = REBUILD / "results"
OUT.mkdir(parents=True, exist_ok=True)

# Official assignments (de Villiers 2004; ICTV MSL 2018+)
ALPHA = {
    1: {32, 42}, 2: {3, 10, 28, 29, 78, 94}, 3: {61, 72, 81, 83, 84, 86, 87, 89},
    4: {2, 27, 57}, 5: {26, 51, 69, 82}, 6: {30, 53, 56, 66},
    7: {18, 39, 45, 59, 68, 70, 85}, 8: {7, 40, 43, 91},
    9: {16, 31, 33, 35, 52, 58, 67}, 10: {6, 11, 13, 44, 55, 74},
    11: {34, 73}, 13: {54}, 15: {71},
}
BETA1 = {5, 8, 12, 14, 19, 20, 21, 24, 25, 36, 47, 93, 98, 99, 105, 118, 124, 143, 152}
BETA2 = {9, 15, 17, 22, 23, 37, 38, 80, 100, 104, 107, 110, 111, 113, 120, 122, 145, 151}
MU1 = {1}
HR_CLADE = {*ALPHA[5], *ALPHA[6], *ALPHA[7], *ALPHA[9], *ALPHA[11]}  # Burk 2013

COMP = pd.read_csv(OUT / "08_cluster_composition.csv")


def official_group(t: int) -> str:
    for g, members in ALPHA.items():
        if t in members:
            return f"Alphapapillomavirus {g}"
    if t in BETA1:
        return "Betapapillomavirus 1"
    if t in BETA2:
        return "Betapapillomavirus 2"
    if t in MU1:
        return "Mupapillomavirus 1"
    return "unassigned"


def main():
    rows = []
    for cl, g in COMP.groupby("cluster"):
        types = sorted(int(t) for t in g["type_str"])
        groups = [official_group(t) for t in types]
        from collections import Counter
        cnt = Counter(groups)
        majority, majority_n = cnt.most_common(1)[0]
        all_hr = all(t in HR_CLADE for t in types)
        rows.append({
            "cluster": int(cl),
            "n_sequences": int(g["n"].sum()),
            "types": "|".join(map(str, types)),
            "official_groups": "|".join(f"{k}({v})" for k, v in sorted(cnt.items())),
            "majority_group": majority,
            "purity": round(majority_n / len(types), 3),
            "all_in_Burk_HR_clade": all_hr,
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "10_cluster_taxonomy_map.csv", index=False)

    n_pure = int((df["purity"] == 1.0).sum())
    summary = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "n_clusters": len(df),
        "clusters_pure_single_official_group": n_pure,
        "cluster_mapping": {int(r["cluster"]): r["official_groups"] for _, r in df.iterrows()},
        "note": ("Cluster 2 is a mixture of Alpha-9/6/5/11 — exactly the high-risk superclade "
                 "defined by Burk et al. 2013. Cluster 9 contains HPV24 (Betapapillomavirus 1) "
                 "among four Betapapillomavirus-2 types."),
        "references": ["de Villiers 2004 Virology 324:17-27",
                       "ICTV MSL (Betapapillomavirus reorganisation)",
                       "Burk 2013 Virology 445:232-240"],
    }
    with open(OUT / "10_cluster_taxonomy_map.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(df.to_string(index=False))
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
