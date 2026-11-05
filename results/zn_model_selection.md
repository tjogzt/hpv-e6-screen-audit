# Zn 模型选择依据（2026-10-29）

## 背景
run1（42 run）中 Zn–S(最近4) 均距 0.825–2.372 nm（中位 0.997 nm），而晶体模板为
0.226–0.236 nm —— 非键 Li/Merz Zn²⁺ + 中性 CYS **不能保持四配位**，锌指折叠塌陷。
这是配体普遍脱离口袋的充分解释，必须先修好 Zn 模型才能重跑筛选。

原 ZAFF 参数档（`zinc_zaff.frcmod`）在远端实例上已不可恢复（`find` 全盘无命中），
故按 ZAFF（Peters et al. 2010）流传参数集重建了一份，并用实测判据验收。

## 实验设计（预注册）
同一 apo 受体（142 残基 + 2 Zn²⁺，TIP3P，12 Å 盒，~31.7k 原子），
EM → 1 ns NVT → 1 ns NPT → 2 ns 生产，同一随机种子 751624，三模型各占一张 RTX A4000 并行。

判据：末 50% 帧内两个 Zn 各自到最近 4 个 Cys SG 的距离均值 ≤ 0.30 nm，
且最大帧均值 ≤ 0.35 nm。（指标定义与 run1 的 `zn_mean_min4_dist` 一致，可直接比较。）

脚本：`code/24_zn_model_test.py`（建模 + MD）、`code/26_zn_analyse.py`（独立分析）。
远端产物：`zn_model_test/<model>/zn_result.json`、`zn_dist_series.csv`。

## 结果

| 模型 | 受体 | Zn 参数 | Zn1 均距 | Zn2 均距 | 双 Zn 均值 | 判定 |
|---|---|---|---|---|---|---|
| `zaff_cym` | rec_cym（8 CYS→CYM） | bonded ZAFF，8 条 ZN–SG 键 | 0.2392 nm | 0.2385 nm | **0.2388 nm** | **PASS** |
| `nonbond_cym` | rec_cym | Li/Merz 12-6 非键，CYM(−1) | 0.1982 nm | 0.1987 nm | 0.1984 nm | PASS（但过紧） |
| `nonbond_cys` | rec_clean（中性 CYS） | Li/Merz 12-6 非键 | 1.0127 nm | 1.0534 nm | 1.0331 nm | **FAIL** |

补充观察：
- `nonbond_cys` 的**首帧**就已到 0.807 / 1.005 nm —— Zn 在能量最小化/预平衡阶段就脱离了配位壳，
  不是缓慢漂移。这完整复现了 run1 的失效模式，证明 run1 的破坏是受体建模问题，不是采样问题。
- `zaff_cym` 全程极稳（max 0.2422 nm，波动 < 0.01 nm），末帧 0.2366 / 0.2375 nm。
- `nonbond_cym` 虽通过阈值，但 1.98 Å 比晶体模板短约 0.3 Å —— 纯静电吸引没有键平衡距离，
  会把锌指压得过紧，可能改变口袋几何。

## 选择
**采用 `zaff_cym`（bonded ZAFF + CYM）** 作为 run2 的 Zn 模型，理由：
1. Zn–S 0.2388 nm 落在晶体实测区间 0.226–0.236 nm 的上沿（差 0.003–0.013 nm），是三个模型中最贴近实验结构的；
2. 轨迹内波动最小（max − mean = 0.0034 nm），配位几何稳定；
3. 过紧的 `nonbond_cym` 会系统性压缩锌指，属于"通过阈值但物理不正确"，不作为首选。

## 待声明事项（写入稿件）
- 所用 ZAFF 参数为**重建**参数集，不是官方 ZAFF 文件的逐字副本；接受依据是实测 Zn–S 距离而非出处。
- 系统净电荷变化：8 个 CYS→CYM 使受体净电荷 −8，配合 Zn²⁺ ×2，体系由 tleap `addions` 自动中和。
- `ZN_S_CUTOFF` 已由 2.8（被当作 Å，实际 28 Å）改为 **0.28 nm**，Zn 门在 run2 才真正生效。
