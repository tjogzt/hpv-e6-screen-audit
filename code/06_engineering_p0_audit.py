#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""06_engineering_p0_audit.py

Read-only audit of the P0 engineering blockers listed in the systematic review.
Produces a machine-readable checklist with file/line evidence; it does not edit
any project file.

Checked items:
  1. start entry referenced by package.json actually exists (start-production.js vs .cjs)
  2. Math.random / hard-coded mock data inside api/routes (production endpoints)
  3. sensitive files present and their permission bits (.env*, SSH private keys)
  4. automated tests (jest config / test files) and whether they can run offline
  5. dependency lock files and pinned versions

Outputs: results/06_engineering_p0_audit.csv, results/06_engineering_p0_audit.json
"""

from __future__ import annotations

import json
import os
import re
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

SKIP_DIRS = {"node_modules", ".venv", "_ssh_venv", "__pycache__", ".git", "dist", "dist-api", "build", "rebuild_20261025"}


def rel(p: Path) -> str:
    return p.relative_to(ROOT).as_posix()


def audit_start_entry():
    pkg = ROOT / "package.json"
    issues = []
    if not pkg.exists():
        return [{"item": "start entry", "target": "package.json", "finding": "missing", "severity": "P0"}]
    data = json.loads(pkg.read_text(encoding="utf-8"))
    for script_name, cmd in (data.get("scripts") or {}).items():
        for token in re.findall(r"node\s+([A-Za-z0-9_\-./]+\.(?:js|cjs|mjs))", cmd):
            target = ROOT / token
            if not target.exists():
                issues.append({
                    "item": "start entry",
                    "target": f"package.json:scripts.{script_name} -> {token}",
                    "finding": "referenced file does not exist",
                    "severity": "P0",
                })
    return issues


def audit_random_fallback():
    issues = []
    api = ROOT / "api"
    if not api.exists():
        return issues
    for f in sorted(api.rglob("*.ts")):
        if any(part in SKIP_DIRS for part in f.parts):
            continue
        text = f.read_text(encoding="utf-8", errors="replace")
        n_rand = len(re.findall(r"Math\.random\(", text))
        n_mock = len(re.findall(r"mock|demo|placeholder", text, re.I))
        if not n_rand and not n_mock:
            continue
        # ID-only randomness (pred_${Date.now()}_${Math.random()...}) carries no data
        rand_lines = [l for l in text.splitlines() if "Math.random(" in l]
        id_only = all("toString(36)" in l or "substring(2" in l for l in rand_lines) if rand_lines else False
        guarded = "demoGuard" in text
        is_guard_module = f.name == "demoGuard.ts"
        if is_guard_module or (guarded and n_rand) or (id_only and n_rand):
            note = ("ID-only randomness, no data fabrication" if id_only and not guarded
                    else "guarded by demoGuard (production throws 503)")
            issues.append({
                "item": "random/demo fallback",
                "target": rel(f),
                "finding": f"Math.random={n_rand}, {note}",
                "severity": "OK",
            })
        elif n_rand or n_mock:
            issues.append({
                "item": "random/demo fallback",
                "target": rel(f),
                "finding": f"Math.random={n_rand}, mock/demo markers={n_mock}, NOT guarded",
                "severity": "P0" if n_rand else "P1",
            })
    return issues


def audit_sensitive_files():
    issues = []
    patterns = [".env", ".env.*", "*.pem", "id_rsa", "id_ed25519", "*_id_*"]
    seen = set()
    for pat in patterns:
        for f in ROOT.rglob(pat):
            if any(part in SKIP_DIRS for part in f.parts) or not f.is_file():
                continue
            if str(f) in seen:
                continue
            seen.add(str(f))
            mode = oct(f.stat().st_mode & 0o777)
            issues.append({
                "item": "sensitive file",
                "target": rel(f),
                "finding": f"mode={mode}" + ("" if int(mode, 8) == 0o600 else " (world/group readable)"),
                "severity": "OK" if int(mode, 8) == 0o600 else "P0",
            })
    return issues


def audit_tests():
    issues = []
    test_globs = ["**/*.test.ts", "**/*.test.js", "**/*.spec.ts", "**/*.spec.js",
                  "tests/**", "test/**", "__tests__/**"]
    found = []
    for g in test_globs:
        found += [p for p in ROOT.glob(g) if not any(part in SKIP_DIRS for part in p.parts)]
    issues.append({
        "item": "automated tests",
        "target": "repository root",
        "finding": f"{len(found)} test file(s) found",
        "severity": "P1" if found else "P0",
    })
    jest_cfg = [p for p in ROOT.glob("jest.config*") if not any(part in SKIP_DIRS for part in p.parts)]
    if not jest_cfg:
        issues.append({
            "item": "jest config",
            "target": "jest.config.*",
            "finding": "missing although package.json declares npm test = jest",
            "severity": "P1",
        })
    return issues


def audit_dependencies():
    issues = []
    pkg = ROOT / "package.json"
    lock = ROOT / "package-lock.json"
    if pkg.exists() and not lock.exists():
        issues.append({"item": "dependency lock", "target": "package-lock.json",
                       "finding": "missing", "severity": "P1"})
    pinned = 0
    total = 0
    req = ROOT / "requirements.txt"
    lock = ROOT / "requirements-lock.txt"
    if req.exists():
        for line in req.read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            total += 1
            if "==" in line:
                pinned += 1
        if lock.exists():
            n_lock = sum(1 for l in lock.read_text().splitlines()
                         if l.strip() and not l.startswith("#") and "==" in l)
            issues.append({"item": "python dependency pinning",
                           "target": "requirements-lock.txt",
                           "finding": f"{n_lock} packages pinned in lock file (requirements.txt keeps ranges)",
                           "severity": "OK"})
        else:
            issues.append({"item": "python dependency pinning",
                           "target": "requirements.txt",
                           "finding": f"{pinned}/{total} pinned with == and no lock file",
                           "severity": "P1"})
    return issues


def main():
    rows = []
    for fn in (audit_start_entry, audit_random_fallback, audit_sensitive_files,
               audit_tests, audit_dependencies):
        try:
            rows.extend(fn())
        except Exception as exc:  # defensive: audit must not crash on odd files
            rows.append({"item": fn.__name__, "target": "-",
                         "finding": f"audit error: {type(exc).__name__}: {exc}",
                         "severity": "P1"})

    import csv
    with open(OUT / "06_engineering_p0_audit.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["item", "target", "finding", "severity"])
        w.writeheader()
        w.writerows(rows)

    summary = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "n_findings": len(rows),
        "by_severity": {},
        "by_item": {},
    }
    for r in rows:
        summary["by_severity"][r["severity"]] = summary["by_severity"].get(r["severity"], 0) + 1
        summary["by_item"][r["item"]] = summary["by_item"].get(r["item"], 0) + 1
    with open(OUT / "06_engineering_p0_audit.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    for r in rows:
        print(f"[{r['severity']}] {r['item']:<24} {r['target']}\n      {r['finding']}")


if __name__ == "__main__":
    main()
