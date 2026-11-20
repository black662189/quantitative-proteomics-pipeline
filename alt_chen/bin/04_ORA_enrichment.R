# ============================================================================
# 三库 ORA 富集分析 (KEGG / GOBP / Reactome) + 气泡图 (nature 紧凑规范)
# 输入: DE_{contrast}_gene_0.5p.xls; universe = 全部定量基因
# 数据库: hsa.kegg_optimal / c5.go.bp.v7.5.1 / c2.cp.reactome.v7.5.1
# 过滤: pvalue<0.05 且 Count>=3
# 出图: nature-figure R 规范 (theme_classic, 6.5pt, Arial, 细线 0.35, 紧凑留白);
#       气泡图配色保留红蓝渐变 (-log10 pvalue)
# ============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(clusterProfiler); library(ggplot2); library(RColorBrewer); library(ragg)
})

args <- commandArgs(TRUE)
if (length(args) != 1)
  stop("用法: Rscript 04_ORA_enrichment.R <region>")
rg <- args[1]

OUT_DE  <- Sys.getenv("PP_OUT_DE")
OUT_ORA <- Sys.getenv("PP_OUT_ORA")
IDMAP   <- Sys.getenv("PP_IDMAP")
PCUT    <- as.numeric(Sys.getenv("PP_PCUT"))
TOPN    <- as.numeric(Sys.getenv("PP_ORA_TOPN"))
CONTRASTS <- strsplit(Sys.getenv("PP_CONTRASTS"), " ")[[1]]

indir <- file.path(OUT_DE, rg)
outdir <- file.path(OUT_ORA, rg)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 通路名 -> 标准短ID 映射 (GO:xxx / hsa04714 / R-HSA-xxx)
id_map_tab <- read.delim(IDMAP, check.names = FALSE, stringsAsFactors = FALSE)
id_map <- setNames(id_map_tab$gs_exact_source, id_map_tab$gs_name)

dbs <- list(
  KEGG     = Sys.getenv("PP_GMT_KEGG"),
  GOBP     = Sys.getenv("PP_GMT_GOBP"),
  Reactome = Sys.getenv("PP_GMT_REACTOME")
)

gene_tab <- read.delim(file.path(indir, sprintf("%s_gene_expression_log2_all.xls", rg)), check.names = FALSE)
universe <- unique(gene_tab$gene)

read_gmt <- function(f) {
  lines <- readLines(f)
  do.call(rbind, lapply(lines, function(l) {
    x <- strsplit(l, "\t", fixed = TRUE)[[1]]
    if (length(x) < 3) return(NULL)
    data.frame(term = x[1], gene = x[-(1:2)], stringsAsFactors = FALSE)
  }))
}
t2g <- lapply(dbs, read_gmt)

fmt_name <- function(id, db) {
  out <- id
  if (db == "GOBP") out <- sub("^GOBP_", "", id)
  if (db == "Reactome") out <- sub("^REACTOME_", "", id)
  if (db != "KEGG") {
    out <- tolower(gsub("_", " ", out))
    out <- paste0(toupper(substr(out, 1, 1)), substr(out, 2, nchar(out)))
  }
  out
}

# 长通路名截断 (保持图不因长名被拉宽)
shorten <- function(x, id = NULL, nchar_max = 40) {
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

# nature 紧凑主题 (6pt Arial 可读, 细线, 无网格)
theme_nature <- theme_classic(base_size = 6, base_family = "Arial") %+replace%
  theme(
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.text = element_text(colour = "black", size = 6),
    axis.title = element_text(colour = "black", size = 6.5),
    plot.title = element_text(colour = "black", size = 7, face = "bold", hjust = 0),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5.5),
    legend.key.size = unit(2.5, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(1, 1, 1, 1, "mm"),
    panel.grid = element_blank()
  )

for (cc in CONTRASTS) {
  deg_file <- file.path(indir, sprintf("DE_%s_gene_0.5p.xls", cc))
  if (!file.exists(deg_file)) next
  deg <- unique(read.delim(deg_file, check.names = FALSE)$gene)
  cat(sprintf("\n%s %s: %d DEG\n", rg, cc, length(deg)))
  if (length(deg) < 5) { cat("  DEG<5 跳过\n"); next }

  for (db in names(dbs)) {
    ek <- enricher(deg, universe = universe, TERM2GENE = t2g[[db]],
                   pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                   minGSSize = 3, maxGSSize = 500)
    ek <- ek[!is.na(ek$pvalue) & ek$pvalue < PCUT & ek$Count >= 3, ]
    df <- as.data.frame(ek)
    if (nrow(df) > 0) df$Description <- fmt_name(df$ID, db)
    write.table(df, file.path(outdir, sprintf("%s_%s_enrichment.xls", cc, db)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  %s: %d 条富集通路\n", db, nrow(df)))
    if (nrow(df) == 0) next

    df <- df[order(df$pvalue), , drop = FALSE]
    df$GeneRatio_num <- vapply(strsplit(df$GeneRatio, "/"),
                               function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
    show <- head(df, TOPN)
    show$mlp <- -log10(show$pvalue)
    show$Description <- shorten(show$Description, if (db == "KEGG") NULL else id_map[show$ID])
    show$Description <- factor(show$Description, levels = rev(show$Description))

    # 紧凑图幅: 名称截断后每行 ~0.09in, 绘图区收窄
    maxn <- max(nchar(as.character(show$Description)))
    w <- maxn * 0.060 + 1.3
    h <- nrow(show) * 0.09 + 0.4

    p <- ggplot(show, aes(GeneRatio_num, Description)) +
      geom_point(aes(size = Count, colour = mlp)) +
      scale_colour_gradientn(colours = rev(brewer.pal(9, "RdYlBu")[c(1:8)]),
                             name = "-log10(pvalue)") +
      scale_size_continuous(range = c(0.6, 2.6), name = "Count") +
      labs(x = "Gene Ratio", y = "") +
      theme_nature

    prefix <- file.path(outdir, sprintf("%s_%s_bubble", cc, db))
    ggsave(paste0(prefix, ".pdf"), p, width = w, height = h, device = cairo_pdf)
    ggsave(paste0(prefix, ".png"), p, width = w, height = h, dpi = 600, device = ragg::agg_png)
    cat(sprintf("    气泡图 top%d (%.2fx%.2f in)\n", nrow(show), w, h))
  }
}
cat(sprintf("%s: 三库 ORA 富集完成\n", rg))
