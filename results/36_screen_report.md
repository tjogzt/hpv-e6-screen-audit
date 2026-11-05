# batch 2 全量筛选汇总（对照校准）

- 明细 48 行（应为 48），配体/对照 16 个
- 对照分离带：DMSO ≤0.005（n=3） < LxxLL ≥0.215（n=3）
- 校准阈值：DMSO-like ≤ 0.050；LxxLL-like ≥ 0.215
- 硬判据：STABLE **0** / MARGINAL **2** / UNSTABLE **12**（14 候选中）
- 对照校准：LxxLL-like **3** / intermediate **4** / DMSO-like **7**

| 对象 | 角色 | 副本QC | 硬判据 | 校准档 | 保留帧率中位 | ΔH 均值±SD (kcal/mol) | Zn合格率(最低) |
|---|---|---|---|---|---|---|---|
| CHEMBL1090_pose1 | candidate | 0/3 | UNSTABLE | DMSO-like | 0.007 | NA | 1.0 |
| CHEMBL1129_pose1 | candidate | 1/3 | MARGINAL | LxxLL-like | 0.2717 | -9.241±NA | 1.0 |
| CHEMBL118218_pose1 | candidate | 0/3 | UNSTABLE | DMSO-like | 0.008 | NA | 1.0 |
| CHEMBL119_pose1 | candidate | 0/3 | UNSTABLE | DMSO-like | 0.017 | NA | 1.0 |
| CHEMBL184_pose1 | candidate | 0/3 | UNSTABLE | intermediate | 0.1239 | NA | 1.0 |
| CHEMBL235668_pose1 | candidate | 1/3 | MARGINAL | LxxLL-like | 0.3467 | 4.339±NA | 1.0 |
| CHEMBL283639_pose1 | candidate | 0/3 | UNSTABLE | DMSO-like | 0.02 | NA | 1.0 |
| CHEMBL29609_pose1 | candidate | 0/3 | UNSTABLE | DMSO-like | 0.011 | NA | 1.0 |
| CHEMBL320775_pose1 | candidate | 0/3 | UNSTABLE | intermediate | 0.2068 | NA | 1.0 |
| CHEMBL3701238_pose1 | candidate | 0/3 | UNSTABLE | DMSO-like | 0.011 | NA | 1.0 |
| CHEMBL4650485_pose1 | candidate | 0/3 | UNSTABLE | intermediate | 0.2078 | NA | 1.0 |
| CHEMBL46917_pose1 | candidate | 0/3 | UNSTABLE | DMSO-like | 0.002 | NA | 1.0 |
| CHEMBL717_pose1 | candidate | 0/3 | UNSTABLE | intermediate | 0.1129 | NA | 1.0 |
| CHEMBL848_pose1 | candidate | 0/3 | UNSTABLE | LxxLL-like | 0.3626 | NA | 1.0 |
| DMSO_pose1 | negative_control | 0/3 | UNSTABLE | DMSO-like | 0.004 | NA | 1.0 |
| LxxLL_peptide | positive_control | 2/3 | MARGINAL | LxxLL-like | 0.5704 | -115.381±0.527 | 1.0 |

## 判据
- 副本级 qc_pass=1：保留帧 ≥500（配体 RMSD<4 Å）且双 Zn 四配位合格率 ≥0.95（管线预注册硬判据）
- 对照校准档：由 DMSO（阴性）与 LxxLL 肽（阳性）实测保留帧率界定，避免 500 帧切点的任意性
- ΔH_PB 仅在 qc_pass=1 的副本上取均值；UNSTABLE 不给 ΔH
