# 全量筛选汇总（14 配体 × 3 replica × 20 ns）

- 明细行数 42，配体 14 个
- STABLE（3/3 副本 QC 通过）**0**；MARGINAL（1–2/3）**4**；UNSTABLE（0/3）**10**

| 配体 | 副本 | QC通过 | 判定 | ΔH 均值±SD (kcal/mol) | 保留帧率均值 | Zn合格率(最低) | 状态 |
|---|---|---|---|---|---|---|---|
| CHEMBL1090_pose1 | 3 | 0 | UNSTABLE | NA | 0.042 | 0.276 | QC_LIG_UNSTABLE |
| CHEMBL1129_pose1 | 3 | 1 | MARGINAL | -14.323±NA | 0.552 | 1.0 | OK;QC_LIG_UNSTABLE |
| CHEMBL118218_pose1 | 3 | 0 | UNSTABLE | NA | 0.062 | 1.0 | QC_LIG_UNSTABLE |
| CHEMBL119_pose1 | 3 | 0 | UNSTABLE | NA | 0.254 | 1.0 | QC_LIG_UNSTABLE |
| CHEMBL184_pose1 | 3 | 1 | MARGINAL | -6.73±NA | 0.53 | 1.0 | OK;QC_LIG_UNSTABLE |
| CHEMBL235668_pose1 | 3 | 0 | UNSTABLE | NA | 0.054 | 0.997 | QC_LIG_UNSTABLE |
| CHEMBL283639_pose1 | 3 | 0 | UNSTABLE | NA | 0.095 | 0.99 | QC_LIG_UNSTABLE |
| CHEMBL29609_pose1 | 3 | 0 | UNSTABLE | NA | 0.219 | 1.0 | QC_LIG_UNSTABLE |
| CHEMBL320775_pose1 | 3 | 0 | UNSTABLE | NA | 0.289 | 0.176 | QC_FAIL;QC_LIG_UNSTABLE |
| CHEMBL3701238_pose1 | 3 | 1 | MARGINAL | -20.552±NA | 0.266 | 1.0 | OK;QC_LIG_UNSTABLE |
| CHEMBL4650485_pose1 | 3 | 0 | UNSTABLE | NA | 0.071 | 0.229 | QC_LIG_UNSTABLE |
| CHEMBL46917_pose1 | 3 | 0 | UNSTABLE | NA | 0.061 | 1.0 | QC_LIG_UNSTABLE |
| CHEMBL717_pose1 | 3 | 0 | UNSTABLE | NA | 0.182 | 0.21 | QC_LIG_UNSTABLE |
| CHEMBL848_pose1 | 3 | 2 | MARGINAL | -16.022±2.094 | 0.457 | 1.0 | OK;QC_LIG_UNSTABLE |

## 判据
- 副本级 qc_pass=1：保留帧 ≥500（配体 RMSD<4 Å）且 双Zn四配位合格率 ≥0.95
- 配体级：3/3 通过=STABLE；1–2/3=MARGINAL；0/3=UNSTABLE（不报 ΔH）
