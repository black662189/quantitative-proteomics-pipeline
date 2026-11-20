# ============================================================================
# 全量 KEGG + GOBP 基因集 GSEA (不踢样本版)
# 输入: 01 的 DE 全表 (18 样本 limma, 全部基因 log2FoldChange) 作为排序统计量
# 基因集: 完整 GOBP (c5.go.bp.v7.5.1) + 完整 KEGG (hsa.kegg_optimal)
# 输出 (4.Diff_Enrichment/{region}/GSEA/):
#   {DB}_{contrast}_fgsea_results.xls   全量 fgsea 结果 (pathway/pval/padj/ES/NES/size/leadingEdge)
#   {DB}_full_GSEA_summary.xls          各对比显著计数汇总 (总/激活/抑制)
# 注: 全量 GSVA 结果沿用 3.GSVA (18 样本), 见 3.GSVA/{rg}_{GOBP,KEGG}_gsva_score.xls (软连接)
# 对比: Model_vs_Normal, HY0094_Model_vs_Model, HY0094_Model_vs_Normal
# ============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(fgsea))
set.seed(42)  # fgseaMultilevel 含随机采样, 固定种子保证可重复

args <- commandArgs(TRUE)
if (length(args) != 1) {
  stop("用法: Rscript 05_full_GSEA.R <region>")
}
rg <- args[1]

OUT_DE  <- Sys.getenv("PP_OUT_DE")
OUT_GSEA <- Sys.getenv("PP_OUT_GSEA")
CONTRASTS <- strsplit(Sys.getenv("PP_CONTRASTS"), " ")[[1]]

indir <- file.path(OUT_DE, rg)
gdir <- file.path(OUT_GSEA, rg)
dir.create(gdir, recursive = TRUE, showWarnings = FALSE)

dbs <- list(GOBP = Sys.getenv("PP_GMT_GOBP"),
            KEGG = Sys.getenv("PP_GMT_KEGG"))

for (db in names(dbs)) {
  gmt <- fgsea::gmtPathways(dbs[[db]])
  summary_rows <- list()
  for (cc in CONTRASTS) {
    de <- read.delim(file.path(indir, sprintf("DE_%s.xls", cc)), check.names = FALSE)
    stats <- setNames(de$log2FoldChange, de$gene)
    stats <- sort(stats[is.finite(stats) & !duplicated(names(stats))], decreasing = TRUE)

    fg <- as.data.frame(fgsea(pathways = gmt, stats = stats, minSize = 3, maxSize = Inf))
    fg <- fg[order(fg$pval), ]
    fg$leadingEdge <- vapply(fg$leadingEdge, paste, character(1), collapse = ";")
    write.table(fg, file.path(gdir, sprintf("%s_%s_fgsea_results.xls", db, cc)),
                sep = "\t", quote = FALSE, row.names = FALSE)

    summary_rows[[cc]] <- data.frame(
      db = db, contrast = cc, n_tested = nrow(fg),
      n_NES_up = sum(fg$NES > 0), n_NES_down = sum(fg$NES < 0),
      n_sig_up = sum(fg$NES > 0 & fg$pval < 0.05),
      n_sig_down = sum(fg$NES < 0 & fg$pval < 0.05),
      mean_NES = round(mean(fg$NES), 4), median_pval = round(median(fg$pval), 4)
    )
    cat(sprintf("%s %s: %d 通路, p<0.05: 激活 %d / 抑制 %d\n", db, cc, nrow(fg),
                summary_rows[[cc]]$n_sig_up, summary_rows[[cc]]$n_sig_down))
  }
  sm <- do.call(rbind, summary_rows)
  write.table(sm, file.path(gdir, sprintf("%s_full_GSEA_summary.xls", db)),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
cat(sprintf("%s: 全量 KEGG + GOBP GSEA 完成\n", rg))
