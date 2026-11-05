# Revision numbers (auto-generated — do not hand-edit)

Generated: 2026-11-04T10:33:11

## R4 — Kruskal-Wallis at two units of analysis
- sequence level: 12/12 FDR-significant; top = hydrophobicity (eps2 = 0.411)
- type level: 10/12 FDR-significant; top = isoelectric_point (eps2 = 0.854)
- hydrophobicity: eps2 0.411 (seq) vs 0.770 (type)
- not significant at type level: instability_index, aromaticity

## R5 — one-vs-rest AUC, ungrouped vs type-grouped CV
- ungrouped range: [0.958, 1.0]
- grouped range: [0.484, 1.0] (8 evaluable clusters)
- single-type clusters (grouped CV undefined): 2

## R6 — cross-k robustness
- k sweep: [4, 6, 8, 10, 12]
- type-level FDR-significant features per k: [9, 10] of 12
- top feature at every k: ['isoelectric_point']
- weighted purity: [0.4222, 0.8667] (k=10: 0.8)
- LOCO evaluable folds per k: {4: 0, 6: 0, 8: 0, 10: 1, 12: 1}

## R2 — selection-free dual-zinc scan
- records scanned: 878; retained 524 (of which 497 carry two complete sites, 27 do not)
- record-level counterexamples: 9 ({'unknown_residue_artifact': 6, 'single_cysteine_deviation': 3})
- types with >= 1 two-site sequence: 45/45
- anchoring failures reproduced: 27 (expected 27)
- verdict: Type-consensus level held (45/45 types carry a record with two complete CXXC-CXXC units); record level does not (9/524 retained records lack one). Title and §3.7b must use the consensus-level statement, not 'absolutely conserved'.
