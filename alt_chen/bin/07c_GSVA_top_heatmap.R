# ============================================================================
# GSVA 两个方向 top 通路热图 (KEGG / GOBP)
# 方向定义: Model vs HY0094_Model 组均值差 (eff = mean(Model) - mean(HY0094_Model))
#   Up   = eff 最大的 top15 (模型激活&给药抑制); Down = eff 最小的 top15 (模型抑制&给药激活)
# 出图: 一张热图, 上半上调 top / 下半下调 top (row_split 分隔), 行 z-score, 三组 18 样本
#   单元格小且正方形, 通路名字体加大(8pt)
# 输出: 4.Diff_Enrichment/{region}/GSVA_GSEA/GSVA_top/{DB}_directional_GSVA_heatmap.{pdf,png}
#       {DB}_directional_GSVA_scores.xls
# ============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize); library(grid); library(ragg)
})

args <- commandArgs(TRUE)
if (length(args) != 1)
  stop("用法: Rscript 07b_GSVA_top_heatmap.R <region>")
rg <- args[1]

OUT_GSVA <- Sys.getenv("PP_OUT_GSVA")
IDMAP    <- Sys.getenv("PP_IDMAP")
TOPN     <- as.numeric(Sys.getenv("PP_GSVA_TOPN"))

groups <- strsplit(Sys.getenv("PP_GROUPS"), " ")[[1]]
# 按组名动态读取 PP_COLOR_<GROUP>; 未配置则回退到固定调色板
group_colors <- vapply(groups, function(g) Sys.getenv(paste0("PP_COLOR_", g)), character(1),
                       USE.NAMES = TRUE)
empty_col <- nchar(group_colors) == 0
if (any(empty_col)) {
  pal <- c("#008B45", "#EE0000", "#3B4992", "#1F77B4", "#FF7F0E", "#2CA02C")
  group_colors[empty_col] <- pal[which(empty_col)]
}
group_levels <- groups
dbs <- c("KEGG", "GOBP")

outdir <- file.path(OUT_GSVA, rg, "GSVA_top")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 通路名 -> 标准短ID 映射
id_map_tab <- read.delim(IDMAP, check.names = FALSE, stringsAsFactors = FALSE)
id_map <- setNames(id_map_tab$gs_exact_source, id_map_tab$gs_name)

fmt_name <- function(p, db) {
  if (db == "GOBP") {
    s <- tolower(gsub("_", " ", sub("^GOBP_", "", p)))
    paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
  } else p
}

shorten <- function(x, id = NULL, nchar_max = 38) {
  y <- x
  cut <- nchar(x) > nchar_max
  if (any(cut)) {
    y[cut] <- paste0(substr(x[cut], 1, nchar_max - 3), "...")
    if (!is.null(id)) {
      y[cut] <- paste0(id[cut], " / ", y[cut])
    }
  }
  y
}

for (db in dbs) {
  f <- file.path(OUT_GSVA, sprintf("%s_%s_gsva_score.xls", rg, db))
  if (!file.exists(f)) next
  tab <- read.delim(f, check.names = FALSE)
  mat <- as.matrix(tab[, -1, drop = FALSE]); rownames(mat) <- tab$path
  storage.mode(mat) <- "double"
  grp_pattern <- paste0("^(", paste(groups, collapse = "|"), ")_.*$")
  grp <- sub(grp_pattern, "\\1", colnames(mat))
  grp <- factor(grp, levels = group_levels)
  idx <- unlist(lapply(group_levels, function(g) which(grp == g)))
  mat <- mat[, idx, drop = FALSE]; grp <- grp[idx]

  # 疾病效应方向: 组均值差 (原始 GSVA score) = 组2 均值 - 组3 均值
  gC <- Sys.getenv("PP_GROUP_CONTROL")  # 对照组
  gT <- Sys.getenv("PP_GROUP_TREAT")    # 实验组
  mT <- rowMeans(mat[, grp == gT, drop = FALSE])  # 实验组均值
  mC <- rowMeans(mat[, grp == gC, drop = FALSE])  # 对照组均值
  eff <- mT - mC   # >0 实验组相对对照上调, <0 下调 (实验组比对照组)

  up   <- head(names(sort(eff, decreasing = TRUE)), TOPN)
  down <- head(names(sort(eff, decreasing = FALSE)), TOPN)
  sel  <- c(up, down)
  direction <- c(rep("Up", length(up)), rep("Down", length(down)))

  sub <- mat[sel, , drop = FALSE]
  rz <- t(scale(t(sub))); rz[!is.finite(rz)] <- 0
  orig_id <- rownames(rz)
  rn <- vapply(rownames(rz), fmt_name, character(1), db = db)
  rn <- shorten(rn, unname(id_map[orig_id]))
  if (anyDuplicated(rn)) rn <- make.unique(rn)
  rownames(rz) <- rn

  eff_col <- paste0("effect_", gT, "_minus_", gC)
  out <- data.frame(path = sel, display_name = rn, direction = direction,
                    x = round(eff[sel], 4), rz, check.names = FALSE)
  names(out)[names(out) == "x"] <- eff_col
  write.table(out, file.path(outdir, sprintf("%s_directional_GSVA_scores.xls", db)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  row_dir <- factor(direction, levels = c("Up", "Down"))

  top_ann <- HeatmapAnnotation(Group = grp,
    col = list(Group = group_colors[group_levels]),
    simple_anno_size = unit(3, "mm"), show_annotation_name = FALSE, show_legend = FALSE)
  col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))

  # 单元格小且正方形 (3mm x 3mm), 列宽=实际样本数 x 3mm, 行高=通路数 x 3mm
  cell <- unit(3, "mm")
  n_col <- ncol(rz)
  name_w <- 62
  ht <- Heatmap(rz, name = "GSVA\nz-score", col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_column_dend = FALSE, show_row_dend = FALSE,
    column_split = grp, column_gap = unit(1, "mm"),
    row_split = row_dir, row_gap = unit(3, "mm"),
    row_title = NULL,
    top_annotation = top_ann,
    width = unit(n_col * 3, "mm"),
    height = unit(nrow(rz) * 3, "mm"),
    show_row_names = TRUE,
    row_names_gp = gpar(fontfamily = "Arial", fontsize = 8, lineheight = 0.9),
    row_names_max_width = unit(name_w, "mm"),
    show_column_names = TRUE,
    column_names_gp = gpar(fontfamily = "Arial", fontsize = 8),
    column_names_rot = 45,
    column_title = NULL,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
      title_gp = gpar(fontfamily = "Arial", fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontfamily = "Arial", fontsize = 8),
      legend_height = unit(20, "mm")))

  fig_w <- (n_col * 3 + name_w + 15 + 14) / 25.4
  fig_h <- (nrow(rz) * 3 + 20 + 26) / 25.4
  prefix <- file.path(outdir, sprintf("%s_directional_GSVA_heatmap", db))
  cairo_pdf(paste0(prefix, ".pdf"), width = fig_w, height = fig_h, family = "Arial")
  draw(ht, heatmap_legend_side = "right", padding = unit(c(2, 20, 2, 4), "mm"))
  dev.off()
  agg_png(paste0(prefix, ".png"), width = fig_w, height = fig_h, units = "in",
          res = 600, background = "white")
  draw(ht, heatmap_legend_side = "right", padding = unit(c(2, 20, 2, 4), "mm"))
  dev.off()
  cat(sprintf("%s %s: 上调 %d + 下调 %d GSVA 热图 (%.1fx%.1f in)\n",
              rg, db, length(up), length(down), fig_w, fig_h))
}
cat(sprintf("%s: GSVA 两方向 top 通路热图完成\n", rg))
