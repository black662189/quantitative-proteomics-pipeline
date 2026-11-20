# 构建各部位 GSVA 输入矩阵 (基因水平, log2, 已填充, 全部 18 样本)
# 输入: 2.standard_data/{region}_protein_expression_log2_imputed.xls (蛋白组水平)
# 输出: 3.GSVA/{region}_log2_gene_GSVA_input.xls (gene/description/样本列)
# 规则: 蛋白组取 Genes 列首个基因符号 (多基因取首个, 无基因符号的行剔除);
#       同一基因对应多个蛋白组时按样本取均值合并 (避免任意丢弃蛋白组)
options(stringsAsFactors = FALSE)

args <- commandArgs(TRUE)
if (length(args) != 1) stop("用法: Rscript 07a_prepare_GSVA_input.R <region>")
rg <- args[1]

OUT_STD <- Sys.getenv("PP_OUT_STD")
OUT_GSVA <- Sys.getenv("PP_OUT_GSVA")

expr <- read.delim(file.path(OUT_STD, sprintf("%s_protein_expression_log2_imputed.xls", rg)),
                   check.names = FALSE)
samples <- colnames(expr)[-(1:4)]

# 基因符号: 取 Genes 列首个; 空值剔除
gene <- sub(";.*$", "", expr$Genes)
gene <- trimws(gene)
keep <- gene != "" & !is.na(gene)
mat <- cbind(data.frame(gene = gene[keep], stringsAsFactors = FALSE),
             expr[keep, -(1:4), drop = FALSE])

# 同基因多蛋白组: 按样本取均值
if (anyDuplicated(mat$gene)) {
  mat <- aggregate(. ~ gene, data = mat, FUN = mean)
}
rownames(mat) <- mat$gene
mat <- mat[, -1, drop = FALSE]

# 校验: 数值完整 (无 NA)
stopifnot(identical(colnames(mat), samples),
          all(vapply(as.matrix(mat), function(x) is.finite(x), logical(1))))

out <- data.frame(gene = rownames(mat), description = "na",
                  as.matrix(mat), check.names = FALSE)
outf <- file.path(OUT_GSVA, sprintf("%s_log2_gene_GSVA_input.xls", rg))
write.table(out, outf, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("%s: %d 蛋白组 -> %d 基因 x %d 样本 (剔除无基因符号 %d, 多基因字段取首个)\n",
            rg, nrow(expr), nrow(mat), ncol(mat), sum(!keep)))
