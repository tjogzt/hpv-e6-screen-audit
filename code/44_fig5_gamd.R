#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# 44_fig5_gamd.R — Figure 5：GaMD 姿态保持与界面结合（§3.9）
#
# 设计原则：**只画不依赖 boost 重加权的量**。
# 原版（v1/v2 早期）用 exp(βΔV) 重加权出 PMF，实测 ESS/N 仅 0.2–3.7%、
# 二阶累积量与指数平均最大偏差 1.4 kcal/mol，PMF 定量不可信；且横坐标
# 「全蛋白 CA 叠合后的配体 RMSD」被 E6 结构域漂移污染（界面区相对全蛋白
# 坐标系漂移 ~19 Å）。故 Figure 5 改为三态分解 + 逐残基占据率：
#
#   面板 A：逐副本配体 RMSD 轨迹（**界面残基 CA 叠合**，主指标），4 Å 参考线
#   面板 B：逐残基接触占据率（偏置系综帧比例，副本间 mean ± SD），
#           E6AP LxxLL 界面残基在底部以刻度标出
#   面板 C：三态分解堆叠条 —— 姿态保持 / 位移但仍接触 / 完全脱离（逐副本）
#
# 数据源：
#   results/<TAG>/<job>__frames.csv        （面板 A，逐帧）
#   results/<TAG>/gamd_summary.csv         （面板 C，三态比例）
#   results/<TAG>/contact_occupancy.csv    （面板 B）
#   results/29_interface_positions.csv     （界面残基注释）
#
# 编号映射：构建体取 4GIZ 2–143 位，故 gro resid = 4GIZ struct_num − 1
# （已用锌配位半胱氨酸核对：gro 29/32/62/65/102/105/135/138
#   ↔ struct 30/33/63/66/103/106/136/139）。
#
# 规范：cairo_pdf 直出目标尺寸；文字 7 pt（panel 标签 9 pt 加粗）；全英文；PDF。

set.seed(42)

argv <- commandArgs(trailingOnly = TRUE)
HERE <- "/Users/taozhu/my researches/HPV_structure/rebuild_20261025"
# 第 1 个参数既接受 TAG（如 gamd_v2），也接受目录全路径（48 收尾脚本传的是路径）
arg1 <- if (length(argv) >= 1) argv[[1]] else "gamd_v2"
if (grepl("/", arg1, fixed = TRUE)) {
  GDIR <- arg1
  TAG <- basename(arg1)
} else {
  TAG <- arg1
  GDIR <- file.path(HERE, "results", TAG)
}
OUT_PDF <- if (length(argv) >= 2) argv[[2]] else
  file.path(HERE, "results", sprintf("Fig5_gamd_contact_pmf_%s.pdf", TAG))

F_SUM <- file.path(GDIR, "gamd_summary.csv")
F_OCC <- file.path(GDIR, "contact_occupancy.csv")
F_IF <- file.path(HERE, "results", "29_interface_positions.csv")
for (f in c(F_SUM, F_OCC, F_IF)) stopifnot(file.exists(f))

summ <- read.csv(F_SUM, stringsAsFactors = FALSE)
occ <- read.csv(F_OCC, stringsAsFactors = FALSE)
ifc <- read.csv(F_IF, stringsAsFactors = FALSE)

LIGS <- c("CHEMBL848", "CHEMBL1129", "LxxLL")
LABS <- c(CHEMBL848 = "CHEMBL848", CHEMBL1129 = "CHEMBL1129",
          LxxLL = "LxxLL peptide")
COLS <- c(CHEMBL848 = "#2166AC", CHEMBL1129 = "#1B7837", LxxLL = "#B2182B")
REPS <- 1:3
SYS_PREFIX <- c(CHEMBL848 = "CHEMBL848_pose1", CHEMBL1129 = "CHEMBL1129_pose1",
                LxxLL = "LxxLL_peptide")

# 界面残基（4GIZ 编号）；occ 的 resid 是 gro resid，故 +1 转回结构编号
# 注意：CSV 里该列是 "True"/"False" 字符串，直接 == TRUE 恒为 FALSE
iface <- ifc$struct_num[tolower(as.character(ifc$e6ap_lxxll_interface)) == "true"]
stopifnot(length(iface) >= 12)

occ$struct_num <- occ$resid + 1
occ$occ_mean <- as.numeric(occ$occ_mean)
occ$occ_sd <- suppressWarnings(as.numeric(occ$occ_sd))

# ---- 读逐帧数据（面板 A）----
frames <- list()
for (L in LIGS) {
  for (r in REPS) {
    job <- sprintf("%s_rep%d", SYS_PREFIX[[L]], r)
    fp <- file.path(GDIR, sprintf("%s__frames.csv", job))
    if (!file.exists(fp)) next
    d <- read.csv(fp, stringsAsFactors = FALSE)
    if (!"lig_rmsd_if" %in% names(d)) stop("frames.csv 缺 lig_rmsd_if 列：", fp)
    frames[[job]] <- data.frame(system = L, replica = r,
                                t_ns = d$time_ps / 1000,
                                rmsd = d$lig_rmsd_if)
  }
}
stopifnot(length(frames) == 9)

cairo_pdf(OUT_PDF, width = 7.2, height = 3.55, pointsize = 7)
layout(matrix(1:3, nrow = 1), widths = c(1.22, 1.26, 1.12))
par(mar = c(4.2, 4.0, 3.8, 0.6), mgp = c(1.85, 0.45, 0), family = "sans",
    tcl = -0.25)

## ---------------- Panel A: pose retention over time ----------------
# 纵轴用平方根刻度：既看得清 0–10 Å 的"姿态保持"区，也不裁掉 40–110 Å 的
# 完全解离帧（线性轴下后者会把前者压成一条线；裁切又会丢信息）。
YMAX <- 110
sq <- function(v) sqrt(v)
plot(NA, NA, xlim = c(0, 100), ylim = sq(c(0, YMAX)), xlab = "", ylab = "",
     bty = "l", xaxt = "n", yaxt = "n")
axis(1, at = seq(0, 100, 25), cex.axis = 1.0)
yt <- c(1, 4, 10, 25, 50, 100)
axis(2, at = sq(yt), labels = yt, cex.axis = 1.0, las = 1)
abline(h = sq(4), col = "grey55", lty = 2, lwd = 0.9)
text(2, sq(4) * 0.78, "4 Å", cex = 0.9, col = "grey35", adj = 0)
for (job in names(frames)) {
  d <- frames[[job]]
  L <- d$system[1]
  lines(d$t_ns, sq(d$rmsd), col = adjustcolor(COLS[[L]], 0.62), lwd = 0.7)
}
mtext("Time (ns)", side = 1, line = 2.5, cex = 1.0)
mtext("Ligand RMSD (Å; sqrt scale)", side = 2, line = 2.7, cex = 1.0)
mtext("A", side = 3, line = 0.55, at = par("usr")[1], font = 2, cex = 1.29,
      adj = 0)

## ---------------- Panel B: per-residue contact occupancy ----------------
occ_p <- occ[occ$struct_num >= 2 & occ$struct_num <= 143, ]
xs <- sort(unique(occ_p$struct_num))

plot(NA, NA, xlim = range(xs), ylim = c(0, 1.02), xlab = "", ylab = "",
     bty = "l", xaxt = "n", yaxt = "n")
axis(1, at = seq(0, 150, 25), cex.axis = 1.0)
axis(2, at = seq(0, 1, 0.25), cex.axis = 1.0, las = 1)
# 界面残基在底边以刻度标出（4GIZ 编号）
ifx <- iface[iface >= min(xs) & iface <= max(xs)]
segments(ifx, 0, ifx, -0.052, col = "#762A83", lwd = 1.4, xpd = NA)
# 用点而非折线：逐残基占据率在 0/1 间跳变，折线会产生大量视觉锯齿
for (L in LIGS) {
  d <- occ_p[occ_p$system == L, ]
  points(d$struct_num, d$occ_mean, pch = 20, cex = 0.42,
         col = adjustcolor(COLS[[L]], 0.75))
}
# 界面 24 残基的平均占据率：左上空白区竖排（不与数据重叠）
text(4, 0.945, "interface mean", cex = 0.85, col = "grey30", adj = c(0, 1))
for (k in seq_along(LIGS)) {
  L <- LIGS[k]
  d <- occ_p[occ_p$system == L & occ_p$struct_num %in% iface, ]
  segments(4, 0.845 - 0.085 * (k - 1), 9, 0.845 - 0.085 * (k - 1),
           col = COLS[[L]], lwd = 2.0, xpd = NA)
  text(11, 0.845 - 0.085 * (k - 1), sprintf("%.2f", mean(d$occ_mean)),
       col = COLS[[L]], cex = 0.92, font = 2, adj = c(0, 0.5))
}
# 界面残基刻度（紫色）的说明写在图注里，避免与坐标轴挤在一起
mtext("E6 residue (4GIZ numbering)", side = 1, line = 2.5, cex = 1.0)
mtext("Contact occupancy", side = 2, line = 2.7, cex = 1.0)
mtext("B", side = 3, line = 0.55, at = par("usr")[1], font = 2, cex = 1.29,
      adj = 0)

## ---------------- Panel C: three-state decomposition ----------------
# 每一副本一根水平堆叠条：姿态保持 / 位移但接触 / 脱离
order_jobs <- unlist(lapply(LIGS, function(L)
  sprintf("%s_rep%d", SYS_PREFIX[[L]], REPS)))
summ$job <- as.character(summ$job)
s <- summ[match(order_jobs, summ$job), ]
stopifnot(all(!is.na(s$job)))
COL_ST <- c("#08306B", "#6BAED6", "#D9D9D9")
par(mar = c(4.2, 4.2, 3.8, 0.4))
plot(NA, NA, xlim = c(0, 1), ylim = c(0.35, 9.65), xlab = "", ylab = "",
     bty = "n", xaxt = "n", yaxt = "n")
# 体系之间留 0.42 行空隙，形成三组
GAP <- 0.42
ypos <- numeric(nrow(s))
for (i in seq_len(nrow(s))) {
  grp <- match(s$system[i], LIGS)
  ypos[i] <- 10 - (i + (grp - 1) * GAP)
}
for (i in seq_len(nrow(s))) {
  f <- c(s$frac_retained[i], s$frac_displaced[i], s$frac_detached[i])
  f[is.na(f)] <- 0
  f <- f / sum(f)
  x0 <- 0
  for (k in 1:3) {
    rect(x0, ypos[i] - 0.33, x0 + f[k], ypos[i] + 0.33, col = COL_ST[k],
         border = "white", lwd = 0.6)
    x0 <- x0 + f[k]
  }
}
SHORT <- c(CHEMBL848 = "848", CHEMBL1129 = "1129", LxxLL = "LxxLL")
axis(2, at = ypos, labels = sprintf("%s  r%d", SHORT[s$system], s$replica),
     las = 1, cex.axis = 0.92, tick = FALSE, line = -0.55, col.axis = "grey10")
# 只标 ≥3% 的"姿态保持"段百分比，避免拥挤
for (i in seq_len(nrow(s))) {
  fr <- s$frac_retained[i]
  if (!is.na(fr) && fr >= 0.03)
    text(fr / 2, ypos[i], sprintf("%.0f", 100 * fr), cex = 0.85,
         col = "white", font = 2)
}
axis(1, at = seq(0, 1, 0.5), cex.axis = 1.0)
mtext("Fraction of production frames", side = 1, line = 2.5, cex = 1.0)
legend(x = 0, y = 10.75, xpd = NA, horiz = TRUE, bty = "n", cex = 0.9,
       fill = COL_ST, border = NA, seg.len = 1.0,
       legend = c("pose retained", "displaced", "detached"))
mtext("C", side = 3, line = 0.55, at = par("usr")[1] - 0.22, font = 2,
      cex = 1.29, adj = 0)

# ---- 共享图例（放设备顶端，避免压住任何面板的数据）----
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
legend(x = grconvertX(0.062, "nfc"), y = grconvertY(0.972, "nfc"),
       xpd = NA, horiz = TRUE, bty = "n", cex = 0.95,
       text.col = COLS[LIGS], col = COLS[LIGS], lwd = 1.3, seg.len = 1.3,
       legend = LABS[LIGS])

dev.off()
info <- file.info(OUT_PDF)
stopifnot(file.exists(OUT_PDF), info$size > 10000)

# ---- 关键数字回显（供 §3.9 文本核对，避免手抄）----
cat(sprintf("[saved] %s (%.1f KB)\n", basename(OUT_PDF), info$size / 1024))
for (L in LIGS) {
  d <- summ[summ$system == L, ]
  cat(sprintf("[%s] RMSD_if mean %s | <=4A %s | detached %s | 全局CA %s | 界面漂移 %s\n",
              L,
              paste(sprintf("%.2f", d$mean_rmsd_if), collapse = "/"),
              paste(sprintf("%.3f", d$frac_rmsd_if_le4), collapse = "/"),
              paste(sprintf("%.3f", d$frac_detached), collapse = "/"),
              paste(sprintf("%.2f", d$mean_rmsd_allca), collapse = "/"),
              paste(sprintf("%.2f", d$mean_if_drift), collapse = "/")))
  o <- occ[occ$system == L & occ$struct_num %in% iface, ]
  cat(sprintf("         interface-residue mean occupancy %.3f (n=%d)\n",
              mean(o$occ_mean), nrow(o)))
}
