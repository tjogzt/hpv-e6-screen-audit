#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""02_clean_ml_taxonomy.py

Row-level taxonomy / type-label audit of results/10_扩展真实HPV_E6数据集.csv.

Rules (all decisions are recorded per row, nothing is deleted):
  os_is_human      : description contains OS=Human papillomavirus*
  os_missing       : no OS= field in description
  type_mismatch_desc: HPV number parsed from the OS/description string != hpv_type
  type_mismatch_acc : HPV number parsed from the accession suffix != hpv_type
  decision
      retain  : human OS present AND no type mismatch (either source)
      exclude : non-human OS, or any confirmed type mismatch
      review  : OS field missing (species cannot be confirmed from the file)

Outputs:
  results/02_ml_taxonomy_audit.csv   one row per input sequence with flags + reason
  data/02_dataset_clean.csv          retained rows only (modeling input)
  results/02_taxonomy_summary.json   counts by decision / reason / type
"""

from __future__ import annotations

import csv
import json
import re
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REBUILD = Path(__file__).resolve().parents[1]
SRC = ROOT / "results" / "10_扩展真实HPV_E6数据集.csv"
OUT_RES = REBUILD / "results"
OUT_DATA = REBUILD / "data"
OUT_RES.mkdir(parents=True, exist_ok=True)
OUT_DATA.mkdir(parents=True, exist_ok=True)

csv.field_size_limit(10 ** 9)

OS_RE = re.compile(r"OS=([^=]+?)(?:\s+OX=|\s+GN=|\s+PE=|\s+SV=|$)")
PV_RE = re.compile(r"([A-Za-z ]*papillomavirus)\s*(?:type\s*)?(\d+)", re.IGNORECASE)
ACC_RE = re.compile(r"_HPV(\d+)", re.IGNORECASE)

FEATURES = [
    "length", "molecular_weight", "aromaticity", "instability_index",
    "helix_propensity", "turn_propensity", "sheet_propensity",
    "hydrophobicity", "isoelectric_point", "net_charge",
    "positive_charged", "negative_charged",
]


def parse_row(r):
    desc = r.get("description") or ""
    acc = r.get("accession") or ""
    hpv_type = str(r.get("hpv_type"))

    m_os = OS_RE.search(desc)
    os_raw = m_os.group(1).strip() if m_os else ""
    os_is_human = bool(re.match(r"^human papillomavirus", os_raw, re.IGNORECASE)) if os_raw else False
    os_missing = not os_raw

    m_pv = PV_RE.search(desc)
    desc_genus = m_pv.group(1).strip() if m_pv else ""
    desc_type = m_pv.group(2) if m_pv else ""

    m_acc = ACC_RE.search(acc)
    acc_type = m_acc.group(1) if m_acc else ""

    mismatch_desc = bool(desc_type) and desc_type != hpv_type
    mismatch_acc = bool(acc_type) and acc_type != hpv_type

    reasons = []
    if os_raw and not os_is_human:
        reasons.append("non_human_os")
    if os_missing:
        reasons.append("os_missing")
    if mismatch_desc:
        reasons.append("type_mismatch_description")
    if mismatch_acc:
        reasons.append("type_mismatch_accession")

    if os_raw and not os_is_human:
        decision = "exclude"
    elif mismatch_desc or mismatch_acc:
        decision = "exclude"
    elif os_missing:
        decision = "review"
    else:
        decision = "retain"

    return {
        "sequence_id": r.get("sequence_id"),
        "hpv_type": hpv_type,
        "risk_level": r.get("risk_level"),
        "pathogenicity": r.get("pathogenicity"),
        "source": r.get("source"),
        "accession": acc,
        "os_raw": os_raw,
        "desc_genus": desc_genus,
        "desc_type": desc_type,
        "acc_type": acc_type,
        "os_is_human": os_is_human,
        "os_missing": os_missing,
        "type_mismatch_desc": mismatch_desc,
        "type_mismatch_acc": mismatch_acc,
        "decision": decision,
        "reason": ";".join(reasons) if reasons else "ok",
    }


def main():
    with open(SRC, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))

    audit = [parse_row(r) for r in rows]

    with open(OUT_RES / "02_ml_taxonomy_audit.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(audit[0].keys()))
        w.writeheader()
        w.writerows(audit)

    keep_idx = [i for i, a in enumerate(audit) if a["decision"] == "retain"]
    out_fields = ["sequence_id", "hpv_type", "risk_level", "pathogenicity", "source",
                  "accession", "sequence"] + FEATURES
    with open(OUT_DATA / "02_dataset_clean.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields)
        w.writeheader()
        for i in keep_idx:
            w.writerow({k: rows[i].get(k, "") for k in out_fields})

    by_decision = {}
    by_reason = {}
    for a in audit:
        by_decision[a["decision"]] = by_decision.get(a["decision"], 0) + 1
        for rsn in (a["reason"].split(";") if a["reason"] != "ok" else ["ok"]):
            by_reason[rsn] = by_reason.get(rsn, 0) + 1

    types_kept = {}
    risk_kept = {}
    for i in keep_idx:
        t = str(rows[i]["hpv_type"])
        types_kept[t] = types_kept.get(t, 0) + 1
        risk_kept[str(rows[i]["risk_level"])] = risk_kept.get(str(rows[i]["risk_level"]), 0) + 1

    summary = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "source": str(SRC.relative_to(ROOT)),
        "n_input_rows": len(rows),
        "by_decision": by_decision,
        "by_reason": by_reason,
        "n_abnormal_rows": by_decision.get("exclude", 0) + by_decision.get("review", 0),
        "retained_rows": len(keep_idx),
        "retained_types": len(types_kept),
        "retained_class_balance": risk_kept,
        "retained_rows_per_type": dict(sorted(types_kept.items(), key=lambda kv: -kv[1])),
        "note": "risk_level/pathogenicity remain dataset labels (type-level proxy); external risk classification not verified in this script",
    }
    with open(OUT_RES / "02_taxonomy_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(json.dumps({k: v for k, v in summary.items() if k != "retained_rows_per_type"}, ensure_ascii=False, indent=2))
    print("retained types:", len(types_kept), "| top:", list(summary["retained_rows_per_type"].items())[:10])


if __name__ == "__main__":
    main()
