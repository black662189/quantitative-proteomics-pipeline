# ============================================================================
# limma 差异表达分析 (通用 pipeline 版)
# 输入: $PP_OUT_STD/{region}_protein_expression_log2_imputed.xls
# 输出 ($PP_OUT_DE/{region}/):
#   {region}_gene_expression_log2_all.xls  基因水平 log2 矩阵
#   DE_{contrast}.xls                      全部基因 + log2FoldChange/pvalue/qvalue
#   DE_{contrast}_gene_0.5p.xls            差异基因 (pvalue<pcut & |log2FC|>=fccut)
# 参数: 区域名 (命令行第1参数); 组别/对比/阈值从环境变量读取
# ============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(limma))

args <- commandArgs(TRUE)
if (length(args) != 1) stop("用法: Rscript 03_DE_limma.R <region>")
rg <- args[1]

OUT_STD <- Sys.getenv("PP_OUT_STD")
OUT_DE  <- Sys.getenv("PP_OUT_DE")
GROUPS  <- strsplit(Sys.getenv("PP_GROUPS"), " ")[[1]]
CONTRASTS <- strsplit(Sys.getenv("PP_CONTRASTS"), " ")[[1]]
PCUT  <- as.numeric(Sys.getenv("PP_PCUT"))
FCCUT <- as.numeric(Sys.getenv("PP_FCCUT"))

infile <- file.path(OUT_STD, sprintf("%s_protein_expression_log2_imputed.xls", rg))
outdir <- file.path(OUT_DE, rg)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

expr <- read.delim(infile, check.names = FALSE)

# 基因水平矩阵: Genes 列第一个符号, 空符号删除, 同基因多蛋白组取均值
genes <- trimws(sub(";.*", "", expr$Genes))
keep <- nzchar(genes)
mat <- as.matrix(expr[keep, -(1:4), drop = FALSE])
storage.mode(mat) <- "double"
genes <- genes[keep]
sm <- rowsum(mat, group = genes)
gene_mat <- sm / as.numeric(table(genes))
if (any(!is.finite(gene_mat))) stop("基因水平矩阵存在非有限值")

grp <- sub(sprintf("_%s-[0-9]+$", rg), "", colnames(gene_mat))
group <- factor(grp, levels = GROUPS)
if (anyNA(group)) stop("存在无法识别分组的样本列: ", paste(colnames(gene_mat)[is.na(group)], collapse = ", "))
cat(sprintf("%s: %d 基因 x %d 样本\n", rg, nrow(gene_mat), ncol(gene_mat)))

# limma 拟合与对比 (组别/对比均从配置解析, 不限组名字符)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(gene_mat, design)

parsed <- lapply(CONTRASTS, function(cc) strsplit(cc, "_vs_")[[1]])
contrast_matrix <- sapply(parsed, function(ab) {
  v <- setNames(rep(0, length(GROUPS)), GROUPS)
  v[ab[1]] <- v[ab[1]] + 1
  v[ab[2]] <- v[ab[2]] - 1
  v
})
colnames(contrast_matrix) <- CONTRASTS
fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))

gene_tab <- data.frame(gene = rownames(gene_mat), gene_mat, check.names = FALSE)
write.table(gene_tab, file.path(outdir, sprintf("%s_gene_expression_log2_all.xls", rg)),
            sep = "\t", quote = FALSE, row.names = FALSE)

for (cc in CONTRASTS) {
  tt <- topTable(fit2, coef = cc, number = Inf, sort.by = "P")
  out <- data.frame(gene = rownames(tt),
                    gene_mat[rownames(tt), , drop = FALSE],
                    log2FoldChange = tt$logFC,
                    pvalue = tt$P.Value,
                    qvalue = tt$adj.P.Val,
                    check.names = FALSE)
  write.table(out, file.path(outdir, sprintf("DE_%s.xls", cc)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  sig <- out[!is.na(out$pvalue) & out$pvalue < PCUT & abs(out$log2FoldChange) >= FCCUT, , drop = FALSE]
  write.table(sig, file.path(outdir, sprintf("DE_%s_gene_0.5p.xls", cc)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  %s: %d/%d 基因 p<%g & |log2FC|>=%g (up: %d, down: %d)\n",
              cc, nrow(sig), nrow(out),
              PCUT, FCCUT,
              sum(sig$log2FoldChange > 0), sum(sig$log2FoldChange < 0)))
}
