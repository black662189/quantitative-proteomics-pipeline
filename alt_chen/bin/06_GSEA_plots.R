# ============================================================================
# 两库 (KEGG / GOBP) 全量 GSEA 出图: 激活 top 与抑制 top 分开的两张柱状图
# 输入: GSEA/{DB}_{contrast}_fgsea_results.xls
# 出图: nature-figure R 规范 (theme_classic 5.5pt Arial 细线 紧凑);
#   激活(NES>0) top10 红色 #EE0000 降序; 抑制(NES<0) top10 蓝色 #3B4992
# 输出: {DB}_{cc}_NES_bar_up.png / _NES_bar_down.png
# ============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(ggplot2); library(ragg) })

args <- commandArgs(TRUE)
if (length(args) != 1)
  stop("用法: Rscript 06_GSEA_plots.R <region>")
rg <- args[1]

OUT_GSEA <- Sys.getenv("PP_OUT_GSEA")
IDMAP    <- Sys.getenv("PP_IDMAP")
PCUT     <- as.numeric(Sys.getenv("PP_PCUT"))
TOPN     <- as.numeric(Sys.getenv("PP_GSEA_BAR_TOPN"))
CONTRASTS <- strsplit(Sys.getenv("PP_CONTRASTS"), " ")[[1]]

gdir <- file.path(OUT_GSEA, rg)
pdir <- file.path(gdir, "plots")
dir.create(pdir, recursive = TRUE, showWarnings = FALSE)

dbs <- c("KEGG", "GOBP")

# 通路名 -> 标准短ID 映射
id_map_tab <- read.delim(IDMAP, check.names = FALSE, stringsAsFactors = FALSE)
id_map <- setNames(id_map_tab$gs_exact_source, id_map_tab$gs_name)

fmt_name <- function(p, db) {
  if (db == "GOBP") {
    s <- tolower(gsub("_", " ", sub("^GOBP_", "", p)))
    paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
  } else p
}

shorten <- function(x, id = NULL, nchar_max = 28) {
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

theme_nature <- theme_classic(base_size = 6, base_family = "Arial") %+replace%
  theme(
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.text = element_text(colour = "black", size = 6),
    axis.title = element_text(colour = "black", size = 6.5),
    plot.title = element_text(colour = "black", size = 7, face = "bold", hjust = 0),
    legend.position = "none",
    plot.margin = margin(1, 1, 1, 1, "mm"),
    panel.grid = element_blank()
  )

# 出单个横向柱状图: 只画对应方向轴(激活0~2 / 抑制-2~0), 柱子区宽度固定, 截断时短ID放前, 无标题
draw_bar <- function(dat, title, fill, fname, dir, with_id = TRUE) {
  sid <- if (with_id) unname(id_map[dat$pathway]) else NULL
  dat$name <- shorten(dat$name, sid)
  dat$name <- factor(dat$name, levels = rev(unique(dat$name)))
  maxn <- max(nchar(as.character(dat$name)))
  # x 轴范围按数据自适应, 留 8% 余量, 避免 NES 柱子被 limits 截断
  xext <- max(abs(dat$NES), na.rm = TRUE) * 1.08
  xlim <- if (dir == "up") c(0, xext) else c(-xext, 0)
  w <- maxn * 0.045 + 1.3
  h <- nrow(dat) * 0.14 + 0.4
  p <- ggplot(dat, aes(NES, name)) +
    geom_col(fill = fill, width = 0.5) +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    labs(x = "NES", y = "") +
    theme_nature
  prefix <- file.path(pdir, fname)
  ggsave(paste0(prefix, ".pdf"), p, width = w, height = h, device = cairo_pdf)
  ggsave(paste0(prefix, ".png"), p, width = w, height = h, dpi = 600, device = ragg::agg_png)
}

for (db in dbs) {
  for (cc in CONTRASTS) {
    f <- file.path(gdir, sprintf("%s_%s_fgsea_results.xls", db, cc))
    if (!file.exists(f)) next
    r <- read.delim(f, check.names = FALSE)
    if (nrow(r) == 0) next
    r$name <- vapply(r$pathway, fmt_name, character(1), db = db)
    sig <- r[r$pval < PCUT & !is.na(r$pval), , drop = FALSE]
    cat(sprintf("%s %s %s: 显著 %d\n", rg, db, cc, nrow(sig)))
    if (nrow(sig) == 0) next

    up <- head(sig[sig$NES > 0, ][order(-sig$NES[sig$NES > 0]), ], TOPN)
    dn <- head(sig[sig$NES < 0, ][order(sig$NES[sig$NES < 0]), ], TOPN)

    all_name <- c(as.character(up$name), as.character(dn$name))

    if (nrow(up) > 0)
      draw_bar(up, sprintf("%s GSEA: %s (%s) activated", db, cc, rg), "#EE0000",
               sprintf("%s_%s_NES_bar_up", db, cc), dir = "up", with_id = TRUE)
    if (nrow(dn) > 0)
      draw_bar(dn, sprintf("%s GSEA: %s (%s) suppressed", db, cc, rg), "#3B4992",
               sprintf("%s_%s_NES_bar_down", db, cc), dir = "down", with_id = TRUE)
  }
}
cat(sprintf("%s: 两库 GSEA 激活/抑制柱状图完成\n", rg))
