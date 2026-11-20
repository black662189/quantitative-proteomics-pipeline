# ============================================================================
# ORA 显著富集通路的通路蛋白热图 (每条通路独立出图) + 映射文件
# 选择: 各 contrast x db 的 top5 显著通路 (不足则全部)
# 出图: 复用 heatmap_light.r (行 z-score, 三组 18 样本, 项目组色)
# 命名: phm_{rg}_{cc}_{db}_{i}.png (i=1..top5, 短名避免长路径)
# 输出 extra: map.tsv (region, contrast, db, term, ID, heatmap_file)
# ============================================================================
options(stringsAsFactors = FALSE)

args <- commandArgs(TRUE)
if (length(args) != 1)
  stop("用法: Rscript 08a_pathway_protein_heatmaps.R <region>")
rg <- args[1]

OUT_DE   <- Sys.getenv("PP_OUT_DE")
OUT_ORA  <- Sys.getenv("PP_OUT_ORA")
OUT_HEAT <- Sys.getenv("PP_OUT_HEATMAP")
RSCRIPT  <- Sys.getenv("PP_RSCRIPT")
HM_SCRIPT <- Sys.getenv("PP_HEATMAP_SCRIPT")
TOPN     <- as.numeric(Sys.getenv("PP_ORA_HEATMAP_TOPN"))
GROUPS   <- strsplit(Sys.getenv("PP_GROUPS"), " ")[[1]]
CONTRASTS <- strsplit(Sys.getenv("PP_CONTRASTS"), " ")[[1]]

indir <- file.path(OUT_DE, rg)
oradir <- file.path(OUT_ORA, rg)
hdir <- file.path(OUT_HEAT, rg, "ORA")
dir.create(hdir, recursive = TRUE, showWarnings = FALSE)

dbs <- c("KEGG", "GOBP", "Reactome")
expr <- read.delim(file.path(indir, sprintf("%s_gene_expression_log2_all.xls", rg)), check.names = FALSE)

map_rows <- list()
made <- 0L
for (cc in CONTRASTS) {
  for (db in dbs) {
    f <- file.path(oradir, sprintf("%s_%s_enrichment.xls", cc, db))
    if (!file.exists(f)) next
    tab <- read.delim(f, check.names = FALSE)
    if (nrow(tab) == 0) next
    tab <- head(tab[order(tab$pvalue), ], TOPN)
    for (i in seq_len(nrow(tab))) {
      genes <- unlist(strsplit(tab$geneID[i], "/"))
      genes <- unique(trimws(genes[nzchar(genes)]))
      genes <- genes[genes %in% expr$gene]
      term <- tab$Description[i]
      if (is.na(term) || !nzchar(term)) term <- tab$ID[i]
      if (length(genes) < 2) {
        map_rows[[length(map_rows) + 1]] <- data.frame(
          region = rg, contrast = cc, db = db, term = term, ID = tab$ID[i],
          heatmap = NA_character_, stringsAsFactors = FALSE)
        next
      }
      fname <- sprintf("phm_%s_%s_%s_%d.png", rg, cc, db, i)
      prefix <- file.path(hdir, sprintf("phm_%s_%s_%s_%d", rg, cc, db, i))
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
        region = rg, contrast = cc, db = db, term = term, ID = tab$ID[i],
        heatmap = fname, stringsAsFactors = FALSE)
      made <- made + 1L
    }
  }
}
map <- do.call(rbind, map_rows)
write.table(map, file.path(hdir, "map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: ORA 通路蛋白热图共 %d 张, 映射 %d 行\n", rg, made, nrow(map)))
