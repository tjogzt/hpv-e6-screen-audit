#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# 35_fig3_startpose_deltaH.R — Figure 3：起始姿态审计 + MM-PBSA 重打分
#
# 输出：results/Fig3_startpose_deltaH.pdf
# 面板 A：每个配体的 MODEL 1（替换前）与选定姿态（替换后）到受体的最小重原子距离，
#         0.36 nm 阈值线；三个被替换的配体以箭头连到其 MODEL 2。
# 面板 B：有 ≥1 个 QC 通过副本的配体的逐副本 ΔH_PB；UNSTABLE 配体不出图；
#         被撤回的 CHEMBL46917 pilot 值（−14.58）以删除线标注。
#
# 数据源：
#   面板 A: results/33_engagement_recomputed.csv（从坐标重算，替代只存替换后值的旧审计 CSV）
#   面板 B: results/21_screen_summary.csv（batch 1）
#
# 规范：cairo_pdf 直出目标尺寸；文字 7 pt（panel 标签 9 pt 加粗）；全英文；PDF。

set.seed(42)

HERE <- "/Users/taozhu/my researches/HPV_structure/rebuild_20261025"
IN_ENG <- file.path(HERE, "results", "33_engagement_recomputed.csv")
OUT_PDF <- file.path(HERE, "results", "Fig3_startpose_deltaH.pdf")

# ΔH 面板数据源：命令行第二个参数可显式指定；否则优先 batch 2 全量
# （36_screen_full.csv），只有在其为完整 48 run 时才自动取用。
pick_src <- function() {
  f36 <- file.path(HERE, "results", "36_screen_full.csv")
  f21 <- file.path(HERE, "results", "21_screen_summary.csv")
  if (file.exists(f36)) {
    n <- nrow(read.csv(f36, stringsAsFactors = FALSE))
    if (n >= 48) return(f36)
    cat(sprintf("[data] batch2 仅 %d/48 run，ΔH 面板暂用 batch1\n", n))
  }
  f21
}
argv <- commandArgs(trailingOnly = TRUE)
IN_SCR <- if (length(argv) >= 2) argv[[2]] else pick_src()
stopifnot(file.exists(IN_ENG), file.exists(IN_SCR))

THRESH_A <- 3.6          # engagement 阈值 (Å)
RETRACT  <- -14.58       # 被撤回的 CHEMBL46917 pilot ΔH_PB

## ---------- 面板 A 数据 ----------
eng <- read.csv(IN_ENG, stringsAsFactors = FALSE)
stopifnot(all(c("ligand", "d_model1_A", "d_selected_model_A",
                "pose_replaced") %in% names(eng)))
eng$before <- as.numeric(eng$d_model1_A)
eng$after  <- as.numeric(eng$d_selected_model_A)
eng$rep    <- as.logical(eng$pose_replaced)
stopifnot(!any(is.na(eng$before)), !any(is.na(eng$after)))

## ---------- 面板 B 数据 ----------
scr <- read.csv(IN_SCR, stringsAsFactors = FALSE)
scr$delta_H <- suppressWarnings(as.numeric(scr$delta_H_PB_kcal_mol))
scr$zn_pass_rate <- as.numeric(scr$zn_pass_rate)
hit <- scr[!is.na(scr$delta_H), ]
stopifnot(nrow(hit) >= 4)

# 色板
col_bef <- "#D6604D"   # 替换前（空点）
col_aft <- "#1B7837"   # 替换后（实点）
col_znfail <- "#762A83" # Zn 门未过的 run

cairo_pdf(OUT_PDF, width = 6.8, height = 3.8, pointsize = 7)
layout(matrix(1:3, nrow = 1), widths = c(1.05, 1.0, 0.34))
par(mar = c(7.0, 3.2, 1.4, 0.5), mgp = c(1.9, 0.5, 0), family = "sans")

## ---------------- Panel A: engagement before/after ----------------
ord <- order(eng$before, decreasing = TRUE)
eng <- eng[ord, ]
n <- nrow(eng)
plot(NA, NA, xlim = c(0.5, n + 0.5), ylim = c(0, 11),
     xaxt = "n", xlab = "", ylab = "Min distance to receptor (A)", main = "")
axis(1, at = 1:n, labels = sub("_pose1$", "", eng$ligand), las = 2, cex.axis = 0.85)
abline(h = THRESH_A, lty = 2, lwd = 0.9, col = "grey30")
text(n - 0.3, THRESH_A + 0.45, "engagement 3.6 A", cex = 0.78, adj = 1, col = "grey30")

# 替换前（全部 14 个 MODEL 1）
points(1:n, eng$before, pch = 1, cex = 0.9, col = col_bef, lwd = 1.1)
# 替换后（仅 3 个 replaced），用箭头从 before 指向 after
for (i in which(eng$rep)) {
  points(i, eng$after[i], pch = 16, cex = 1.0, col = col_aft)
  arrows(i, eng$before[i], i, eng$after[i] + 0.35,
         length = 0.06, lwd = 1.0, col = col_aft)
}
# 非替换的 14 个的 after == before（同一点），不重复画
legend("topright", bty = "n", cex = 0.82,
       pch = c(1, 16), pt.cex = 1.0,
       col = c(col_bef, col_aft),
       legend = c("MODEL 1 (before)", "MODEL 2 (after replacement)"))
mtext("ligand", side = 1, line = 6.1, cex = 0.95)
mtext("A", side = 3, line = 0.3, at = 0.5, font = 2, cex = 1.29, adj = 0)

## ---------------- Panel B: per-replica deltaH ----------------
par(mar = c(7.0, 3.4, 1.4, 0.6))
## 补救打分（同一 PB 流程，补算给「过阳性对照下限但硬门未过」的副本）
F_RESC <- file.path(HERE, "results", "39_rescue_deltaH.csv")
if (file.exists(F_RESC)) {
  resc0 <- read.csv(F_RESC, stringsAsFactors = FALSE)
  resc <- resc0[!is.na(suppressWarnings(as.numeric(resc0$delta_H_PB_kcal_mol))), ]
  resc$delta_H_PB_kcal_mol <- as.numeric(resc$delta_H_PB_kcal_mol)
  cat(sprintf("[rescue] %d/%d 补算成功\n", nrow(resc), nrow(resc0)))
} else {
  resc <- data.frame(ligand = character(), delta_H_PB_kcal_mol = numeric(),
                     replica = integer(), kept_frac = numeric())
}

CTL <- c("LxxLL_peptide", "DMSO_pose1")   # 对照单独进面板 C
hit_c <- hit[!hit$CHEMBL %in% CTL, ]
resc_c <- if (nrow(resc) > 0) resc[!resc$ligand %in% CTL, ] else resc
ligs <- sort(unique(c(unique(hit_c$CHEMBL), unique(resc_c$ligand))))
hit_c$xfac <- factor(hit_c$CHEMBL, levels = ligs)
set.seed(42)
xpos <- as.numeric(hit_c$xfac) + jitter(rep(0, nrow(hit_c)), factor = 0.18)
zn_ok <- hit_c$zn_pass_rate >= 0.95
col_resc <- "#1B7837"   # 补救值（探索性）

ylim <- range(c(hit_c$delta_H, RETRACT, resc_c$delta_H_PB_kcal_mol), na.rm = TRUE)
ylim[1] <- ylim[1] - 2.5; ylim[2] <- ylim[2] + 2.5
n_lig <- length(ligs)
xr <- n_lig + 0.7          # 被撤回值占独立 x 位（46917 无有效 ΔH）
plot(NA, NA, xlim = c(0.5, n_lig + 1.2), ylim = ylim,
     xaxt = "n", xlab = "",
     ylab = expression(Delta * H[PB] ~ "(kcal/mol)"), main = "")
axis(1, at = c(seq_along(ligs), xr),
     labels = c(sub("_pose1$", "", ligs), "CHEMBL46917"),
     las = 2, cex.axis = 0.9)
# 通过 Zn 门的实点；未过的用紫色空心（batch 2 实际 0 例，故不画空图例）
points(xpos[zn_ok], hit_c$delta_H[zn_ok], pch = 16, cex = 1.0, col = "#2166AC")
n_znfail <- sum(!zn_ok)
if (n_znfail > 0) {
  points(xpos[!zn_ok], hit_c$delta_H[!zn_ok], pch = 1, cex = 1.0,
         col = col_znfail, lwd = 1.2)
}
# 补救打分探索值（相同 PB 流程，针对硬门未过但保留率 ≥ 阳性对照下限的副本）
if (nrow(resc_c) > 0) {
  rx <- as.numeric(factor(resc_c$ligand, levels = ligs)) +
        jitter(rep(0, nrow(resc_c)), factor = 0.30)
  points(rx, resc_c$delta_H_PB_kcal_mol, pch = 22, cex = 1.1,
         col = col_resc, bg = adjustcolor(col_resc, 0.30), lwd = 1.0)
}
# 被撤回的 CHEMBL46917 pilot 值（独立 x 位 xr，上文已定义）
points(xr, RETRACT, pch = 4, cex = 1.1, col = "grey50", lwd = 1.4)
text(xr, RETRACT + 1.2, "retracted", cex = 0.75, col = "grey50")
segments(xr - 0.22, RETRACT, xr + 0.22, RETRACT, col = "grey50", lwd = 1.2)

mtext("ligand", side = 1, line = 6.1, cex = 0.95)
lgd_pch <- c(16, 4)
lgd_col <- c("#2166AC", "grey50")
lgd_txt <- c("QC-passing replica", "retracted value")
if (n_znfail > 0) {                      # batch 2 无 Zn 门失败，故默认不出现
  lgd_pch <- c(lgd_pch, 1)
  lgd_col <- c(lgd_col, col_znfail)
  lgd_txt <- c(lgd_txt, "Zn gate failed")
}
if (nrow(resc_c) > 0) {
  lgd_pch <- c(lgd_pch, 22)
  lgd_col <- c(lgd_col, col_resc)
  lgd_txt <- c(lgd_txt, "rescued replica (exploratory)")
}
legend("bottomleft", bty = "n", cex = 0.72, pch = lgd_pch, pt.cex = 1.0,
       col = lgd_col, legend = lgd_txt)
mtext("B", side = 3, line = 0.3, at = 0.5, font = 2, cex = 1.29, adj = 0)

## ---------------- Panel C: controls ----------------
# 阳性对照 ≈ −115 kcal/mol 会把候选区压平，故对照单独成面板（不同 y 量程）
par(mar = c(7.0, 3.0, 1.4, 0.5))
pos_rows <- hit[hit$CHEMBL == "LxxLL_peptide" & !is.na(hit$delta_H), ]
pos_resc <- if (nrow(resc) > 0) resc[resc$ligand == "LxxLL_peptide", ] else NULL
y_all <- c(pos_rows$delta_H,
           if (!is.null(pos_resc) && nrow(pos_resc)) pos_resc$delta_H_PB_kcal_mol)
ylimc <- c(min(y_all) - 4, max(y_all) + 4)
plot(NA, NA, xlim = c(0.55, 2.05), ylim = ylimc, xaxt = "n", xlab = "",
     ylab = expression(Delta * H[PB] ~ "(kcal/mol)"), main = "")
segments(par("usr")[1], 0, par("usr")[2], 0, col = "grey80", lty = 2)
points(jitter(rep(1, nrow(pos_rows)), factor = 0.7), pos_rows$delta_H,
       pch = 16, cex = 1.0, col = "#2166AC")
if (!is.null(pos_resc) && nrow(pos_resc)) {
  points(rep(1.18, nrow(pos_resc)), pos_resc$delta_H_PB_kcal_mol, pch = 22,
         cex = 1.1, col = col_resc, bg = adjustcolor(col_resc, 0.30), lwd = 1.0)
}
axis(1, at = c(1.05, 1.75), labels = c("LxxLL", "DMSO"), las = 2, cex.axis = 0.9)
text(1.75, mean(ylimc), "not scored\n(0/3 QC)", cex = 0.68, col = "grey40",
     adj = 0.5)
text(par("usr")[2], max(ylimc), "note: different y-scale", cex = 0.62,
     col = "grey40", adj = 1)
mtext("control", side = 1, line = 6.1, cex = 0.95)
mtext("C", side = 3, line = 0.3, at = 0.5, font = 2, cex = 1.29, adj = 0)

dev.off()
info <- file.info(OUT_PDF)
stopifnot(file.exists(OUT_PDF), info$size > 10000)
cat(sprintf("[saved] %s (%.1f KB)\n", OUT_PDF, info$size / 1024))
cat(sprintf("[check] replaced ligands: %s\n",
            paste(eng$ligand[eng$rep], collapse = ", ")))
