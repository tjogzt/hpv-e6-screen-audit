#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# 34_fig2_zn_pose.R — Figure 2：单结构筛选的锌位点完整性与 pose 保留
#
# 输出：results/Fig2_zn_pose_retention.pdf
# 面板 A：每个 run 的 pose 保留率（14 配体 × 3 副本），按中位排序；0.499 阈值线。
# 面板 B：每个 run 的 Zn 到最近 4 个 Cys S 的平均距离，晶体模板参考带 0.226–0.236 nm。
#
# 数据源：默认 results/21_screen_summary.csv（batch 1）；batch 2 全量
# （results/36_screen_full.csv，48 run）就绪后自动切换，也可命令行显式传路径。
#
# 规范：cairo_pdf 直出目标尺寸；文字 7 pt（panel 标签 9 pt 加粗）；全英文；PDF。

set.seed(42)  # 抖动点

HERE <- "/Users/taozhu/my researches/HPV_structure/rebuild_20261025"
OUT_PDF <- file.path(HERE, "results", "Fig2_zn_pose_retention.pdf")

# 数据源：命令行第一个参数可显式指定；否则优先用 batch 2 全量文件
# （results/36_screen_full.csv，列名与 batch 1 的 21_screen_summary.csv 兼容，
#  多出 ctrl_tier 等列），但只有在其为完整 48 run 时才自动取用。
pick_src <- function() {
  f36 <- file.path(HERE, "results", "36_screen_full.csv")
  f21 <- file.path(HERE, "results", "21_screen_summary.csv")
  if (file.exists(f36)) {
    n <- nrow(read.csv(f36, stringsAsFactors = FALSE))
    if (n >= 48) return(f36)
    cat(sprintf("[data] batch2 仅 %d/48 run，暂用 batch1\n", n))
  }
  f21
}
argv <- commandArgs(trailingOnly = TRUE)
IN_CSV <- if (length(argv) >= 1) argv[[1]] else pick_src()

stopifnot(file.exists(IN_CSV))
d <- read.csv(IN_CSV, stringsAsFactors = FALSE)
# 防御性：关键列存在且无 NA
stopifnot(all(c("CHEMBL", "replica", "kept_frac", "zn_mean_min4_dist") %in% names(d)))
d$kept_frac <- as.numeric(d$kept_frac)
d$zn_mean_min4_dist <- as.numeric(d$zn_mean_min4_dist)
stopifnot(!any(is.na(d$kept_frac)), !any(is.na(d$zn_mean_min4_dist)))

THRESH <- 0.499            # 预注册阈值：≥500/1001 帧
BAND   <- c(0.226, 0.236)  # 晶体模板 Zn–S 范围 (nm)

# 配体按 pose 保留率中位数排序（两面板共用同一顺序，便于对照）
med_keep <- tapply(d$kept_frac, d$CHEMBL, median)
lig_order <- names(sort(med_keep))
d$CHEMBL <- factor(d$CHEMBL, levels = lig_order)
shortlab <- sub("_pose1$", "", lig_order)

cat(sprintf("[data] %d runs, %d ligands; kept 中位=%.3f; Zn 中位=%.3f 范围 %.3f-%.3f\n",
            nrow(d), length(lig_order), median(d$kept_frac),
            median(d$zn_mean_min4_dist),
            min(d$zn_mean_min4_dist), max(d$zn_mean_min4_dist)))

# 色板
col_pass <- "#B2182B"   # 红：超阈值的 run
col_fail <- "#8C8C8C"   # 灰
col_zn   <- "#2166AC"   # 蓝：Zn 距离
col_band <- "#D1E5F0"

cairo_pdf(OUT_PDF, width = 6.8, height = 3.8, pointsize = 7)
layout(matrix(1:2, nrow = 1), widths = c(1, 1))
par(mar = c(7.0, 3.2, 1.4, 0.5), mgp = c(1.9, 0.5, 0), family = "sans")

## ---------------- Panel A: pose retention ----------------
xj <- jitter(as.numeric(d$CHEMBL), factor = 0.5)
above <- d$kept_frac >= THRESH
plot(xj, d$kept_frac, pch = 16, cex = 0.7, xaxt = "n",
     col = ifelse(above, col_pass, adjustcolor(col_fail, 0.8)),
     xlab = "", ylab = "Pose retention (fraction of 1001 frames)",
     ylim = c(0, 1), xlim = c(0.5, length(lig_order) + 0.5), main = "")
axis(1, at = seq_along(lig_order), labels = shortlab, las = 2, cex.axis = 0.9)
abline(h = THRESH, lty = 2, lwd = 0.8, col = "grey30")
# 阳性对照下限：若含对照数据，从 LxxLL 三副本现算；否则不画
pos_floor <- if ("LxxLL_peptide" %in% d$CHEMBL) {
  min(d$kept_frac[d$CHEMBL == "LxxLL_peptide"])
} else NA
if (!is.na(pos_floor)) {
  abline(h = pos_floor, lty = 3, lwd = 0.9, col = "#B2182B")
  text(1, pos_floor + 0.035, sprintf("positive-control floor %.3f", pos_floor),
       cex = 0.72, col = "#B2182B", adj = 0)
}
# 每个配体的中位短横线
for (i in seq_along(lig_order)) {
  v <- d$kept_frac[d$CHEMBL == lig_order[i]]
  segments(i - 0.28, median(v), i + 0.28, median(v), lwd = 1.6, col = "grey20")
}
mtext("ligand", side = 1, line = 6.1, cex = 0.95)
text(1, THRESH + 0.055, "threshold 0.499", cex = 0.8, col = "grey30", adj = 0)
mtext("A", side = 3, line = 0.3, at = 0.5, font = 2, cex = 1.29, adj = 0)

## ---------------- Panel B: Zn–4S distance ----------------
par(mar = c(7.0, 3.2, 1.4, 0.5))
# 晶体参考带
plot(NA, NA, xlim = c(0.5, length(lig_order) + 0.5),
     ylim = c(0, max(d$zn_mean_min4_dist) * 1.08), xaxt = "n",
     xlab = "", ylab = expression("Zn" ~ "to 4 nearest Cys S (nm)"), main = "")
rect(par("usr")[1], BAND[1], par("usr")[2], BAND[2],
     col = col_band, border = NA)
abline(h = BAND, col = "#2166AC", lty = 3, lwd = 0.8)
points(xj, d$zn_mean_min4_dist, pch = 16, cex = 0.7,
       col = adjustcolor(col_zn, 0.85))
axis(1, at = seq_along(lig_order), labels = shortlab, las = 2, cex.axis = 0.9)
mtext("ligand", side = 1, line = 6.1, cex = 0.95)
text(length(lig_order) - 0.5, BAND[2] + 0.12,
     "crystal 0.226–0.236 nm", cex = 0.8, col = "#2166AC", adj = 1)
mtext("B", side = 3, line = 0.3, at = 0.5, font = 2, cex = 1.29, adj = 0)

dev.off()

info <- file.info(OUT_PDF)
stopifnot(file.exists(OUT_PDF), info$size > 10000)
cat(sprintf("[saved] %s (%.1f KB)\n", OUT_PDF, info$size / 1024))
