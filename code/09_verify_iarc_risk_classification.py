#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""09_verify_iarc_risk_classification.py

Verify the dataset's type-level risk labels against official IARC classifications
for the 45 HPV types retained after taxonomy cleaning (02_dataset_clean.csv).

Reference (authoritative): IARC Scientific Council report quoting Monographs
vol. 100B (events.iarc.who.int, GC66 SC report):
  Group 1 (carcinogenic)      : 16,18,31,33,35,39,45,51,52,56,58,59
  Group 2A (probably)         : 68
  Group 2B (possibly)         : 26,53,66,67,70,73,82
  Group 3 (not classifiable)  : 6,11 + beta/gamma genus types
  HPV 5,8                     : limited evidence for non-melanoma skin cancer
                                (epidermodysplasia verruciformis patients)

Known secondary-source conflict on 68 (2025 Nanjing study lists it as 2B;
PMC12113743 claims reclassification to Group 1) -> official IARC report wins,
conflict recorded in the note column.

Consistency rule:
  dataset high_risk  == IARC Group 1                      -> consistent
  dataset high_risk  vs  Group 2A/2B                      -> flagged (label too strong)
  dataset low_risk   vs  Group 2A/2B                      -> flagged (label too weak)
  dataset low_risk   ==  Group 3 / unclassified / other   -> consistent
  type not evaluated by IARC                              -> consistent_unofficial (dataset label
                                                             cannot be verified; flagged for review)

Outputs:
  results/09_iarc_risk_verification.csv   one row per retained type
  results/09_iarc_risk_verification.json  summary + discrepancy list
"""

from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REBUILD = Path(__file__).resolve().parents[1]
DATA = REBUILD / "data" / "02_dataset_clean.csv"
OUT = REBUILD / "results"
OUT.mkdir(parents=True, exist_ok=True)

IARC_GROUP1 = {16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59}
IARC_GROUP2A = {68}
IARC_GROUP2B = {26, 53, 66, 67, 70, 73, 82}
IARC_GROUP3 = {6, 11}  # explicitly evaluated, not classifiable
NMSC_LIMITED = {5, 8}

SOURCE_IARC = ("IARC Monographs vol.100B via IARC Scientific Council report "
               "(events.iarc.who.int/event/74/attachments/189/1294/GC66_5_SC_Report.pdf)")
CONFLICT_NOTE = ("secondary sources conflict: Nanjing 2025 (PMC12512415) lists 68 as 2B; "
                 "PMC12113743 claims 68 reclassified Group 1; official IARC SC report = 2A")


def iarc_group(t: int):
    if t in IARC_GROUP1:
        return "1"
    if t in IARC_GROUP2A:
        return "2A"
    if t in IARC_GROUP2B:
        return "2B"
    if t in IARC_GROUP3:
        return "3"
    if t in NMSC_LIMITED:
        return "limited (non-melanoma skin cancer, EV patients)"
    return "not evaluated"


def classify(t: int, dataset_risk: str):
    g = iarc_group(t)
    high = dataset_risk == "high_risk"
    if high and g == "1":
        return g, "consistent", ""
    if high and g in ("2A", "2B"):
        return g, "flag_label_too_strong", CONFLICT_NOTE if t == 68 else ""
    if (not high) and g in ("2A", "2B"):
        return g, "flag_label_too_weak", ""
    if not high and g in ("3", "limited (non-melanoma skin cancer, EV patients)"):
        return g, "consistent", ""
    if g == "not evaluated":
        return g, "consistent_unofficial", "IARC has not evaluated this type; dataset label unverifiable"
    return g, "consistent", ""


def main():
    with open(DATA, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    by_type = {}
    for r in rows:
        t = int(r["hpv_type"])
        by_type.setdefault(t, {"risk": r["risk_level"], "n": 0})
        by_type[t]["n"] += 1

    out_rows = []
    for t in sorted(by_type):
        info = by_type[t]
        g, status, note = classify(t, info["risk"])
        out_rows.append({
            "hpv_type": t,
            "n_sequences": info["n"],
            "dataset_risk_level": info["risk"],
            "iarc_group": g,
            "consistency": status,
            "note": note,
            "source": SOURCE_IARC,
        })

    with open(OUT / "09_iarc_risk_verification.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)

    flagged = [r for r in out_rows if r["consistency"].startswith("flag")]
    unverified = [r for r in out_rows if r["consistency"] == "consistent_unofficial"]
    summary = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "n_types": len(out_rows),
        "n_consistent": len([r for r in out_rows if r["consistency"] == "consistent"]),
        "n_flagged": len(flagged),
        "n_unverified": len(unverified),
        "flagged_types": flagged,
        "unverified_types": [r["hpv_type"] for r in unverified],
        "iarc_group1_in_dataset": sorted(IARC_GROUP1 & set(by_type)),
        "source": SOURCE_IARC,
        "caveat": ("'consistent' means the dataset label agrees with IARC Group 1 vs non-Group-1 "
                   "as a binary high/low split; Group 2A/2B types in either dataset class are flagged. "
                   "consistency_unofficial types have no IARC evaluation and keep the dataset label "
                   "only as a working label."),
    }
    with open(OUT / "09_iarc_risk_verification.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(json.dumps({k: v for k, v in summary.items()
                      if k not in ("flagged_types",)}, ensure_ascii=False, indent=2))
    print("flagged:", [(r["hpv_type"], r["dataset_risk_level"], r["iarc_group"])
                       for r in flagged])


if __name__ == "__main__":
    main()
