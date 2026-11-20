#!/usr/bin/env Rscript
# 表型 leading edge 基因热图（RdBu 浅色版）
# 配色参考 PE_PE94_lipid_metabolism_pathways_clustered_GSVA_heatmap：
#   colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))
# 用法: Rscript heatmap_light.r <expr.xls> <output_prefix> <group1> <group2> [...]
#   expr.xls: 行=基因, 列=样本（列名形如 <group>-<数字>），FPKM 表达值
# 输出: <output_prefix>.pdf (cairo_pdf) + <output_prefix>.png (600 dpi)

suppressPackageStartupMessages({
  library(pheatmap)
  library(gtable)
  library(grid)
  library(ragg)
})
library(showtext)

# 注册 Arial；pheatmap 默认文本无 family，需递归改写 gtable 文本字体
# 字体目录: 优先 PP_FONT_DIR(pipeline 自包含 lib/fonts/Arial), 回退 ~/fonts/Arial
font_dir <- Sys.getenv("PP_FONT_DIR")
if (!nzchar(font_dir)) font_dir <- path.expand("~/fonts/Arial")
font_add("Arial",
         regular    = file.path(font_dir, "ARIAL.TTF"),
         bold       = file.path(font_dir, "ARIALBD.TTF"),
         italic     = file.path(font_dir, "ARIALI.TTF"),
         bolditalic = file.path(font_dir, "ARIALBI.TTF"))
showtext_auto()

# 递归设置 gtable 中所有 text grob 的字体族为 Arial
apply_font <- function(g) {
  if (inherits(g, "text")) {
    gp <- g$gp
    if (is.null(gp)) gp <- gpar()
    gp$fontfamily <- "Arial"
    g$gp <- gp
  } else if (inherits(g, "gTree")) {
    g$children <- lapply(g$children, apply_font)
  }
  g
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("用法: heatmap_light.r <expr.xls> <output_prefix> <group1> <group2> [...] [--anno-colors=#hex1,#hex2,...]")
}
# 可选 --anno-colors: 指定分组注释条配色 (默认浅蓝/浅橙), 顺序同分组参数
custom_cols <- NULL
flag_idx <- grep("^--anno-colors=", args)
if (length(flag_idx)) {
  custom_cols <- strsplit(sub("^--anno-colors=", "", args[flag_idx[1]]), ",")[[1]]
  args <- args[-flag_idx]
}
# 可选 --title: 通路名称, 作为热图主标题(顶部), 不进入分组/图例
main_title <- ""
title_idx <- grep("^--title=", args)
if (length(title_idx)) {
  main_title <- sub("^--title=", "", args[title_idx[1]])
  args <- args[-title_idx]
}
# 标题过长折行 (约 60 字符), 避免超出画布
if (nzchar(main_title) && nchar(main_title) > 60) {
  main_title <- paste(strwrap(main_title, width = 60), collapse = "\n")
}

expr <- read.csv(args[1], header = TRUE, sep = "\t", check.names = FALSE, quote = "")
expr <- expr[!duplicated(expr[, 1]), ]
rownames(expr) <- expr[, 1]
expr <- expr[, -1]

# 仅保留指定分组样本列，分组顺序即图中列顺序
groups <- args[-c(1, 2)]
mcol <- unlist(lapply(groups, function(g) grep(paste0("^", g, "-[0-9]+$"), colnames(expr))))
if (length(mcol) == 0) stop("未匹配到分组样本列: ", paste(groups, collapse = ", "))
expr <- expr[, mcol]

# 过滤零表达行
expr <- expr[apply(expr, 1, mean) > 0, , drop = FALSE]
if (nrow(expr) < 2) stop("可用基因行数不足 2: ", args[1])

# 基因名保持原样（本项目为人源基因符号，大写规范；参考脚本的首字母大写仅适用于鼠源命名）

# 列注释：样本分组
sam <- sub("-[^-]+$", "", colnames(expr))
anno_col <- data.frame(sample = factor(sam, levels = groups))
rownames(anno_col) <- colnames(expr)

# 分组注释条配色: 默认浅蓝/浅橙; 可经 --anno-colors 覆盖 (如项目组色 Normal/Model/HY0094_Model)
group_palette <- if (!is.null(custom_cols)) custom_cols else
  c("#6497B1", "#FDAE6B", "#005B96", "#E6550D", "#B3CDE3", "#FDD0A2")
anno_colors <- list(sample = setNames(group_palette[seq_along(groups)], groups))

# 行 z-score，截断到 +/-1（沿用原 heatmap.r 设定，仅更换配色）
mat <- t(scale(t(expr)))
mat[mat >= 1] <- 1
mat[mat <= -1] <- -1

# RdBu 浅色渐变（参考 PE_PE94 GSVA 热图配色，替代原 blue4/red4 深色系）
heat_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

# 构建 gtable 前打开空设备并保持到绘制结束：
# pheatmap 预构建与 gtable 单位换算均需活动设备，否则 grid 自动生成杂散 Rplots.pdf
grDevices::pdf(NULL)
gt <- pheatmap(mat,
               cluster_rows = FALSE, cluster_cols = FALSE,
               cellwidth = 12, cellheight = 8, fontsize = 6,
               main = if (nzchar(main_title)) main_title else NULL,
               fontsize_row = 6,
               annotation_col = anno_col,
               annotation_colors = anno_colors,
               border_color = NA,
               color = heat_colors,
               silent = TRUE)

# pheatmap 返回包裹对象，实际 gtable 在 $gtable 中
grob_tbl <- if (inherits(gt, "pheatmap")) gt$gtable else gt
grob_tbl$grobs <- lapply(grob_tbl$grobs, apply_font)

# 从 gtable 推算实际画布尺寸，另加 0.1 英寸边距防止图例文字字形下缘贴边
fig_w <- grid::convertWidth(gtable::gtable_width(grob_tbl), "inch", valueOnly = TRUE) + 0.1
fig_h <- grid::convertHeight(gtable::gtable_height(grob_tbl), "inch", valueOnly = TRUE) + 0.1

cairo_pdf(paste0(args[2], ".pdf"), width = fig_w, height = fig_h)
grid.draw(grob_tbl); dev.off()

showtext_opts(dpi = 600)
agg_png(paste0(args[2], ".png"), width = fig_w, height = fig_h,
        units = "in", res = 600, background = "white")
grid.draw(grob_tbl); dev.off()
grDevices::dev.off()  # 关闭开头打开的空设备

cat("完成:", args[2], ".pdf/.png |", nrow(mat), "基因 x", ncol(mat), "样本 |",
    "分组:", paste(groups, collapse = " vs "), "| 尺寸",
    round(fig_w, 2), "x", round(fig_h, 2), "英寸\n")
