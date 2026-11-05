#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""05_fix_candidate_energy_parser.py

Corrected MMPBSA.py FINAL_RESULTS.dat parser + regression test against the
14 archived candidate results.

Bug being fixed (see .workbuddy/candidates_pipeline.py):
  * the "existing results" branch calls parse_dg_bind (correct: sums the
    Differences block),
  * the FRESH-RUN branch (lines ~403-411) reads `DELTA TOTAL`, `DELTA G gas`,
    `DELTA G solv` with a regex. In AmberTools 24 / MMPBSA.py 14.0 PB-only mode
    those summary lines are printed as 0.0000, so every fresh run is recorded as
    dG = 0.0000 with status OK.

Parser contract:
  * dG_PB = VDWAALS + EEL + INTERNAL + EPB + ENPOLAR + EDISPER (components only)
  * raises ParseError when the Differences block is missing -> no silent "NA"+OK
  * records std. dev. propagation and the known-zero DELTA TOTAL trap

Regression test:
  * reads candidates_run_final_20260818.tar.gz (14 FINAL_RESULTS.dat)
  * compares with .workbuddy/SUMMARY_final.csv (manual recomputation)
  * also reports what the buggy fresh-run regex would have produced

Outputs: results/05_parser_regression.csv, results/05_parser_regression.json
Exit code 1 if any regression fails.
"""

from __future__ import annotations

import csv
import json
import math
import re
import tarfile
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REBUILD = Path(__file__).resolve().parents[1]
OUT = REBUILD / "results"
OUT.mkdir(parents=True, exist_ok=True)

ARCHIVE = ROOT / "candidates_run_final_20260818.tar.gz"
REFERENCE = ROOT / ".workbuddy" / "SUMMARY_final.csv"

PB_COMPONENTS = ("VDWAALS", "EEL", "INTERNAL", "EPB", "ENPOLAR", "EDISPER")
TOL = 0.01  # kcal/mol


class ParseError(RuntimeError):
    pass


def parse_dg_bind_strict(txt: str):
    """Return dict with dG, components, std devs and trap flags."""
    m = re.search(r"Differences \(Complex - Receptor - Ligand\):(.*?)DELTA TOTAL", txt, re.S)
    if not m:
        raise ParseError("no 'Differences (Complex - Receptor - Ligand)' block found")

    comps, sds = {}, {}
    for line in m.group(1).splitlines():
        mm = re.match(r"\s*([A-Z][A-Z0-9\- ]*?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
        if not mm:
            continue
        name = mm.group(1).strip()
        if name.upper().startswith("DELTA"):
            continue
        comps[name] = float(mm.group(2))
        sds[name] = float(mm.group(3))

    if not comps:
        raise ParseError("Differences block found but no energy components parsed")

    unknown = [c for c in comps if c not in PB_COMPONENTS]
    dg = sum(comps[c] for c in PB_COMPONENTS if c in comps)
    # conservative uncertainty: square root of summed variances of used components
    var = sum(sds[c] ** 2 for c in PB_COMPONENTS if c in comps)
    sd = math.sqrt(var)

    mt = re.search(r"DELTA TOTAL\s+(-?\d+\.\d+)", txt)
    total = float(mt.group(1)) if mt else None
    trap = bool(total == 0.0 and abs(dg) > 1e-6)

    return {
        "dG_PB_kcal_mol": round(dg, 4),
        "sd_kcal_mol": round(sd, 4),
        "n_frames": _frames(txt),
        "components": comps,
        "unknown_components": unknown,
        "delta_total_reported": total,
        "delta_total_zero_trap": trap,
    }


def _frames(txt: str):
    m = re.search(r"using\s+([\d.]+)\s+complex frames", txt)
    return float(m.group(1)) if m else None


def buggy_fresh_run_value(txt: str):
    """Reproduce the regex used by the fresh-run branch of candidates_pipeline.py."""
    m = re.search(r"DELTA TOTAL\s*=\s*([-\d.]+)\s*±\s*([-\d.]+)", txt)
    if m:
        return float(m.group(1))
    m2 = re.search(r"DELTA TOTAL\s+(-?\d+\.\d+)", txt)
    return float(m2.group(1)) if m2 else None


def load_reference():
    ref = {}
    with open(REFERENCE, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            ref[r["CHEMBL"]] = float(r["dG_PB_kcal_mol"])
    return ref


def main():
    ref = load_reference()
    rows, fails = [], []

    with tarfile.open(ARCHIVE) as tar:
        members = sorted([m for m in tar.getmembers() if m.name.endswith("FINAL_RESULTS.dat")],
                         key=lambda x: x.name)
        for mem in members:
            chembl = Path(mem.name).parent.name
            txt = tar.extractfile(mem).read().decode("utf-8", errors="replace")
            try:
                parsed = parse_dg_bind_strict(txt)
            except ParseError as e:
                fails.append(f"{chembl}: {e}")
                rows.append({"CHEMBL": chembl, "status": "PARSE_ERROR", "error": str(e)})
                continue

            expected = ref.get(chembl)
            delta = None if expected is None else abs(parsed["dG_PB_kcal_mol"] - expected)
            ok = expected is not None and delta <= TOL
            if not ok:
                fails.append(f"{chembl}: parsed {parsed['dG_PB_kcal_mol']} vs reference {expected}")

            rows.append({
                "CHEMBL": chembl,
                "dG_PB_kcal_mol": parsed["dG_PB_kcal_mol"],
                "sd_kcal_mol": parsed["sd_kcal_mol"],
                "n_frames": parsed["n_frames"],
                "reference_dG": expected,
                "abs_delta": None if delta is None else round(delta, 4),
                "buggy_fresh_run_value": buggy_fresh_run_value(txt),
                "delta_total_zero_trap": parsed["delta_total_zero_trap"],
                "unknown_components": "|".join(parsed["unknown_components"]),
                "status": "PASS" if ok else "FAIL",
            })

    # negative controls: the parser must fail loudly, not return NA/0
    for label, sample in (("empty", ""), ("no_block", "MMPBSA.py Version=14.0\n")):
        try:
            parse_dg_bind_strict(sample)
            fails.append(f"negative control {label}: expected ParseError")
            rows.append({"CHEMBL": f"NEGATIVE_{label}", "status": "FAIL",
                         "error": "expected ParseError, none raised"})
        except ParseError:
            rows.append({"CHEMBL": f"NEGATIVE_{label}", "status": "PASS",
                         "error": "ParseError raised as required"})

    with open(OUT / "05_parser_regression.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)

    summary = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "archive": str(ARCHIVE.relative_to(ROOT)),
        "reference": str(REFERENCE.relative_to(ROOT)),
        "n_candidates": len([r for r in rows if r["CHEMBL"].startswith("CHEMBL")]),
        "n_pass": len([r for r in rows if r.get("status") == "PASS"]),
        "n_fail": len([r for r in rows if r.get("status") == "FAIL"]),
        "tolerance_kcal_mol": TOL,
        "failures": fails,
        "all_zero_trap": all(r.get("delta_total_zero_trap") for r in rows
                             if r["CHEMBL"].startswith("CHEMBL")),
    }
    with open(OUT / "05_parser_regression.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(json.dumps({k: summary[k] for k in
                      ["n_candidates", "n_pass", "n_fail", "all_zero_trap", "failures"]},
                     ensure_ascii=False, indent=2))
    for r in rows:
        if r["CHEMBL"].startswith("CHEMBL"):
            print(f"{r['CHEMBL']:<14} dG={r['dG_PB_kcal_mol']:>8} ref={r['reference_dG']:>8} "
                  f"buggy={r['buggy_fresh_run_value']} {r['status']}")

    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
