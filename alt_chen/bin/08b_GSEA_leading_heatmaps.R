# ============================================================================
# 两库 (KEGG/GOBP) GSEA 显著通路的领头亚基 (leading edge) 基因热图 + 映射
# 输入: GSEA/{DB}_{contrast}_fgsea_results.xls (leadingEdge 列 = ; 分隔领头基因)
# 选择: 各 contrast x db 的 top8 显著通路 (pval<0.05, 不足则全部)
# 出图: 复用 heatmap_light.r (行 z-score 三组 18 样本, 项目组色)
# 命名: lead_{rg}_{cc}_{db}_{i}.png; 输出 map: lead_map.tsv
# ============================================================================
options(stringsAsFactors = FALSE)

args <- commandArgs(TRUE)
if (length(args) != 1)
  stop("用法: Rscript 08b_GSEA_leading_heatmaps.R <region>")
rg <- args[1]

OUT_DE   <- Sys.getenv("PP_OUT_DE")
OUT_GSEA <- Sys.getenv("PP_OUT_GSEA")
OUT_HEAT <- Sys.getenv("PP_OUT_HEATMAP")
RSCRIPT  <- Sys.getenv("PP_RSCRIPT")
HM_SCRIPT <- Sys.getenv("PP_HEATMAP_SCRIPT")
TOPN     <- as.numeric(Sys.getenv("PP_GSEA_LEAD_TOPN"))
PCUT     <- as.numeric(Sys.getenv("PP_PCUT"))
GROUPS   <- strsplit(Sys.getenv("PP_GROUPS"), " ")[[1]]
CONTRASTS <- strsplit(Sys.getenv("PP_CONTRASTS"), " ")[[1]]

indir <- file.path(OUT_DE, rg)
gdir <- file.path(OUT_GSEA, rg)
ldir <- file.path(OUT_HEAT, rg, "GSEA")
dir.create(ldir, recursive = TRUE, showWarnings = FALSE)

dbs <- c("KEGG", "GOBP")
expr <- read.delim(file.path(indir, sprintf("%s_gene_expression_log2_all.xls", rg)), check.names = FALSE)

fmt_name <- function(p, db) {
  if (db == "GOBP") {
    s <- tolower(gsub("_", " ", sub("^GOBP_", "", p)))
    paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
  } else p
}

map_rows <- list()
made <- 0L
for (db in dbs) {
  for (cc in CONTRASTS) {
    f <- file.path(gdir, sprintf("%s_%s_fgsea_results.xls", db, cc))
    if (!file.exists(f)) next
    r <- read.delim(f, check.names = FALSE)
    r <- r[!is.na(r$pval) & r$pval < PCUT, , drop = FALSE]
    if (nrow(r) == 0) next
    r <- head(r[order(r$pval), ], TOPN)
    for (i in seq_len(nrow(r))) {
      term <- fmt_name(r$pathway[i], db)
      genes <- unique(unlist(strsplit(as.character(r$leadingEdge[i]), ";", fixed = TRUE)))
      genes <- genes[nzchar(genes) & genes %in% expr$gene]
      if (length(genes) < 2) {
        map_rows[[length(map_rows) + 1]] <- data.frame(
          region = rg, contrast = cc, db = db, pathway = term,
          heatmap = NA_character_, stringsAsFactors = FALSE)
        next
      }
      fname <- sprintf("lead_%s_%s_%s_%d.png", rg, cc, db, i)
      prefix <- file.path(ldir, sprintf("lead_%s_%s_%s_%d", rg, cc, db, i))
      sub <- expr[match(genes, expr$gene), , drop = FALSE]
      write.table(sub, paste0(prefix, ".xls"), sep = "\t", quote = FALSE, row.names = FALSE)
      gcols <- vapply(GROUPS, function(g) Sys.getenv(paste0("PP_COLOR_", g)), character(1))
      if (any(nchar(gcols) == 0)) {
        pal <- c("#008B45", "#EE0000", "#3B4992", "#1F77B4", "#FF7F0E", "#2CA02C")
        gcols[nchar(gcols) == 0] <- pal[which(nchar(gcols) == 0)]
      }
      anno_colors <- paste(gcols, collapse = ",")
      hm_args <- c(shQuote(HM_SCRIPT), shQuote(paste0(prefix, ".xls")), shQuote(prefix),
                   sprintf("%s_%s", GROUPS, rg),
                   paste0("--anno-colors=", anno_colors),
                   shQuote(paste0("--title=", db, ": ", term)))
      system2(RSCRIPT, hm_args)
      map_rows[[length(map_rows) + 1]] <- data.frame(
        region = rg, contrast = cc, db = db, pathway = term,
        heatmap = fname, stringsAsFactors = FALSE)
      made <- made + 1L
    }
  }
}
map <- do.call(rbind, map_rows)
write.table(map, file.path(ldir, "lead_map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: GSEA 领头亚基热图共 %d 张, 映射 %d 行\n", rg, made, nrow(map)))
