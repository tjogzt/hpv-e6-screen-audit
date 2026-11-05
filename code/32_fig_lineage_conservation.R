#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# 32_fig_lineage_conservation.R — 功能位点谱系保守性图（支持 manuscript §3.7b）
#
# 输出：results/Fig4_lineage_conservation.pdf（Figure 4）
# 面板 A：沿 E6 序列（UniProt 30–150）的 lineage coherence，按位点类别着色；
#         Zn 配位位点（8 个，全部不变，coherence 无定义）在顶部标记为 100% Cys。
# 面板 B：E6AP / p53 / 背景 三组的 coherence 分布（箱线 + 抖动点），
#         标注 Mann-Whitney U 置换检验 p 值与事后功效。
#
# 规范：cairo_pdf 直出目标尺寸；全部文字 7 pt（panel 标签 9 pt 加粗）；全英文；PDF。

suppressPackageStartupMessages({
  library(jsonlite)
})

# ---- 参数集中管理 ----
set.seed(42)  # 抖动点随机种子
HERE <- "/Users/taozhu/my researches/HPV_structure/rebuild_20261025"
IN_CSV  <- file.path(HERE, "results", "31_position_coherence.csv")
IN_JSON <- file.path(HERE, "results", "31_functional_vs_background.json")
OUT_PDF <- file.path(HERE, "results", "Fig4_lineage_conservation.pdf")

stopifnot(file.exists(IN_CSV), file.exists(IN_JSON))

# ---- 读数据 ----
d <- read.csv(IN_CSV, stringsAsFactors = FALSE)
d$site_class <- as.character(d$site_class)
d$lineage_coherence <- as.numeric(d$lineage_coherence)
d$position <- as.numeric(d$position)
d$invariant <- as.logical(d$invariant)

# 防御性：无 NA/NaN/Inf
stopifnot(!any(is.na(d$lineage_coherence)))

# 组别拆分
is_zn  <- d$site_class == "zn_coordination"
is_e6ap <- d$site_class == "e6ap_lxxll_interface"
is_p53  <- d$site_class == "p53_interface"
is_both <- d$site_class == "e6ap_lxxll_interface+p53_interface"
is_bg   <- d$site_class == "background"

# 面板 B 的可比较组（排除 invariant，coherence 对不变位点无定义）
grp_e6ap <- d$lineage_coherence[is_e6ap & !d$invariant]
grp_p53  <- d$lineage_coherence[is_p53 & !d$invariant]
grp_bg   <- d$lineage_coherence[is_bg & !d$invariant]

# 从 JSON 读组间检验结果（p、功效）
js <- read_json(IN_JSON, simplifyVector = TRUE)
e6ap_p   <- js$groups$e6ap_lxxll_interface$p_raw
e6ap_pow <- js$groups$e6ap_lxxll_interface$post_hoc_power
p53_p    <- js$groups$p53_interface$p_raw
p53_pow  <- js$groups$p53_interface$post_hoc_power
bg_med   <- js$groups$background$median_coherence

# ---- 颜色（全局统一，色盲友好）----
col_zn   <- "#B2182B"   # 红：Zn 配位（保守）
col_e6ap <- "#2166AC"   # 蓝：E6AP 界面
col_p53  <- "#1B7837"   # 绿：p53 界面
col_both <- "#762A83"   # 紫：E6AP∩p53 重叠
col_bg   <- "#8C8C8C"   # 灰：背景

# ---- 画布：并排两面板，等高 ----
# 目标单面板约 3.3 in 宽 → 总宽 6.8 in，高 3.3 in
cairo_pdf(OUT_PDF, width = 6.8, height = 3.3, pointsize = 7)
layout(matrix(1:2, nrow = 1), widths = c(1.35, 1))
par(mar = c(3.2, 3.4, 1.6, 0.6), mgp = c(2.0, 0.6, 0),
    cex.lab = 1.0, cex.axis = 1.0, cex.main = 1.0, family = "sans")

## ---------------- Panel A: 沿序列的 coherence ----------------
xlim <- c(30, 150)
ylim <- c(-0.10, 1.02)

# 背景点
plot(d$position[is_bg], d$lineage_coherence[is_bg],
     pch = 16, cex = 0.55, col = adjustcolor(col_bg, 0.55),
     xlim = xlim, ylim = ylim, xaxt = "n",
     xlab = "", ylab = "Lineage coherence", main = "")
# x 轴
axis(1, at = seq(30, 150, 30), cex.axis = 1.0)
mtext("E6 position (UniProt P03126)", side = 1, line = 1.9, cex = 1.0)
# 背景中位数水平线
abline(h = bg_med, col = col_bg, lty = 2, lwd = 0.8)
# 重叠点
if (any(is_both)) {
  points(d$position[is_both], d$lineage_coherence[is_both],
         pch = 17, cex = 0.85, col = col_both)
}
# E6AP / p53
points(d$position[is_e6ap], d$lineage_coherence[is_e6ap],
       pch = 16, cex = 0.85, col = col_e6ap)
points(d$position[is_p53], d$lineage_coherence[is_p53],
       pch = 16, cex = 0.85, col = col_p53)
# Zn 配位位点（invariant）：顶部菱形标记
zn_pos <- d$position[is_zn]
points(zn_pos, rep(0.97, length(zn_pos)),
       pch = 23, bg = col_zn, col = col_zn, cex = 0.9)
# 顶部说明条
text(90, 1.005, "Zn-coordinating Cys: 100% Cys in all 45 types",
     cex = 0.85, col = col_zn)
# 图例
legend("bottomright", bty = "n", cex = 0.82, pt.cex = 0.9,
       pch = c(16, 16, 17, 16, 23),
       pt.bg = c(NA, NA, NA, NA, col_zn),
       col = c(col_e6ap, col_p53, col_both, col_bg, col_zn),
       legend = c("E6AP LxxLL interface", "p53 interface",
                  "both interfaces", "background",
                  "Zn-coordination (invariant)"))
mtext("A", side = 3, line = 0.4, at = xlim[1], font = 2, cex = 1.29, adj = 0)

## ---------------- Panel B: 组间分布 ----------------
par(mar = c(3.2, 3.4, 1.6, 0.4))
grp_list <- list(E6AP = grp_e6ap, p53 = grp_p53, Background = grp_bg)
cols_grp <- c(col_e6ap, col_p53, col_bg)
bp <- boxplot(grp_list, outline = FALSE, col = adjustcolor(cols_grp, 0.35),
              border = cols_grp, lwd = 0.9, xaxt = "n", ylim = c(-0.1, 0.8),
              ylab = "Lineage coherence", xlab = "", main = "",
              boxwex = 0.55, staplewex = 0.4)
axis(1, at = 1:3, labels = c("E6AP", "p53", "Background"), cex.axis = 1.0)
# 抖动点
for (i in 1:3) {
  xjit <- jitter(rep(i, length(grp_list[[i]])), factor = 0.22)
  points(xjit, grp_list[[i]], pch = 16, cex = 0.6,
         col = adjustcolor(cols_grp[i], 0.7))
}
# 中位数水平参考
abline(h = bg_med, col = col_bg, lty = 2, lwd = 0.8)
# 检验标注
text(1, 0.72, sprintf("p=%.2f", e6ap_p), cex = 0.85)
text(2, 0.72, sprintf("p=%.2f", p53_p), cex = 0.85)
text(1, 0.64, sprintf("power=%.2f", e6ap_pow), cex = 0.75, col = "grey35")
text(2, 0.64, sprintf("power=%.2f", p53_pow), cex = 0.75, col = "grey35")
mtext("B", side = 3, line = 0.4, at = 0.5, font = 2, cex = 1.29, adj = 0)

dev.off()

# ---- 校验 ----
info <- file.info(OUT_PDF)
stopifnot(file.exists(OUT_PDF), info$size > 10000)
cat(sprintf("[saved] %s (%.1f KB)\n", OUT_PDF, info$size / 1024))
cat(sprintf("[check] E6AP n=%d p=%.3f power=%.2f | p53 n=%d p=%.3f power=%.2f | bg n=%d med=%.3f\n",
            length(grp_e6ap), e6ap_p, e6ap_pow,
            length(grp_p53), p53_p, p53_pow,
            length(grp_bg), bg_med))
