# ============================================================================
# 蛋白质组标准数据预处理 (通用 pipeline 版)
# 输入:  $PP_MATRIX (DIA-NN report.pg_matrix.tsv, 含 Protein.Group/Genes 等注释列 + 样本定量列)
# 输出:  $PP_OUT_STD/ 下各部位的 log2 与 log2+填充 矩阵 + sample_metadata.xls + process_QA.txt
# 处理:  蛋白过滤(任一 部位x组别 定量>=floor(replicates*2/3)) -> log2 -> 线性最小值/2 填充
#        基因水平按 Genes 字段(下游脚本另行处理)
# ============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(readr))

ROOT     <- Sys.getenv("PP_ROOT")
MATRIX   <- Sys.getenv("PP_MATRIX")
REG_CSV  <- Sys.getenv("PP_REG_CSV")
OUT_STD  <- Sys.getenv("PP_OUT_STD")
REGIONS  <- strsplit(Sys.getenv("PP_REGIONS"), " ")[[1]]
GROUPS   <- strsplit(Sys.getenv("PP_GROUPS"), " ")[[1]]
REP      <- as.integer(Sys.getenv("PP_REPLICATES"))
dir.create(OUT_STD, showWarnings = FALSE, recursive = TRUE)

qa_lines <- character()
qa <- function(...) qa_lines <<- c(qa_lines, sprintf(...))

# ---------------- 样本映射构建 (REGIONS x GROUPS x REP) ----------------
map <- data.frame()
k <- 0
for (tis in REGIONS) for (cond in GROUPS) for (i in 1:REP) {
  k <- k + 1
  map <- rbind(map, data.frame(region = tis, group = cond, replicate = i,
                               sample = sprintf("%s_%s-%d", cond, tis, i),
                               stringsAsFactors = FALSE))
}
n_samp <- nrow(map)
qa("样本映射: %d 样本 = %d 部位 x %d 组 x %d 重复", n_samp, length(REGIONS), length(GROUPS), REP)

# ---------------- 读取 DIA-NN 蛋白组矩阵 ----------------
mat_raw <- read.delim(MATRIX, check.names = FALSE, stringsAsFactors = FALSE)
anno_cols <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description")
stopifnot(all(anno_cols %in% colnames(mat_raw)))
prot_info <- mat_raw[, anno_cols, drop = FALSE]

# 定量列: 前 6 列是注释(Protein.Group/Names/Genes/First.Protein.Description/N.Sequences/N.Proteotypic.Sequences)
quant_cols <- mat_raw[, -(1:6), drop = FALSE]
stopifnot(ncol(quant_cols) == n_samp)

colnames(quant_cols) <- map$sample
expr_linear <- as.matrix(quant_cols)
storage.mode(expr_linear) <- "double"
rownames(expr_linear) <- prot_info$Protein.Group

# ---------------- 蛋白过滤: 任一 部位x组别 组合中定量 >= floor(REP*2/3) ----------------
min_valid <- max(1L, floor(REP * 2 / 3))
group_key <- paste(map$group, map$region, sep = "_")
valid_counts <- sapply(unique(group_key), function(g) {
  rowSums(!is.na(expr_linear[, group_key == g, drop = FALSE]))
})
keep_proteins <- rowSums(valid_counts >= min_valid) > 0
qa("蛋白过滤: %d/%d 蛋白在至少一个 部位x组别 中定量 >=%d/%d 样本, 保留",
   sum(keep_proteins), nrow(expr_linear), min_valid, REP)
qa("  过滤前缺失比例: %.2f%%", 100 * mean(is.na(expr_linear)))
expr_linear <- expr_linear[keep_proteins, , drop = FALSE]
prot_info <- prot_info[keep_proteins, , drop = FALSE]
qa("  过滤后缺失比例: %.2f%%", 100 * mean(is.na(expr_linear)))
qa("  基因注释缺失: %d/%d 蛋白无 Genes 符号", sum(is.na(prot_info$Genes) | prot_info$Genes == ""), nrow(prot_info))

# ---------------- log2 转换与缺失值填充 ----------------
expr_log2 <- log2(expr_linear)
min_nonzero <- min(expr_linear[!is.na(expr_linear)])
fill_linear <- min_nonzero / 2
fill_log2 <- log2(fill_linear)
qa("缺失值填充: 线性最小值 %.4g, 填充值 %.4g (log2 = %.3f), 填充缺失 %d 个",
   min_nonzero, fill_linear, fill_log2, sum(is.na(expr_log2)))
expr_imputed <- expr_log2
expr_imputed[is.na(expr_imputed)] <- fill_log2
stopifnot(!any(is.na(expr_imputed)))

# ---------------- 输出 ----------------
write_meta <- data.frame(
  sample = map$sample,
  region = map$region,
  group = map$group,
  replicate = map$replicate,
  stringsAsFactors = FALSE)
write.table(write_meta, file.path(OUT_STD, "sample_metadata.xls"),
            sep = "\t", quote = FALSE, row.names = FALSE)

for (tis in REGIONS) {
  sidx <- map$region == tis
  sub_samples <- unlist(lapply(GROUPS, function(g) map$sample[sidx & map$group == g]))
  out_na <- data.frame(prot_info, round(expr_log2[, sub_samples, drop = FALSE], 4),
                       check.names = FALSE)
  write.table(out_na, file.path(OUT_STD, sprintf("%s_protein_expression_log2.xls", tis)),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  out_imp <- data.frame(prot_info, round(expr_imputed[, sub_samples, drop = FALSE], 4),
                        check.names = FALSE)
  write.table(out_imp, file.path(OUT_STD, sprintf("%s_protein_expression_log2_imputed.xls", tis)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  gm <- sapply(GROUPS, function(g) {
    cols <- map$sample[sidx & map$group == g]
    sprintf("%s %.1f%%", g, 100 * mean(is.na(expr_log2[, cols, drop = FALSE])))
  })
  qa("%s: %d 蛋白 x %d 样本 | 缺失率: %s", tis, nrow(out_na), length(sub_samples), paste(gm, collapse = ", "))
  # 逐样本缺失率 (过滤后, 填充前)
  sm <- sprintf("%s %.1f%%", sub_samples, 100 * colMeans(is.na(expr_log2[, sub_samples, drop = FALSE])))
  qa("  样本缺失率: %s", paste(sm, collapse = " | "))
}

qa("log2 数值范围: %.2f ~ %.2f (填充后)", min(expr_imputed), max(expr_imputed))
qa("部位蛋白集一致 (全局过滤): %d 蛋白", nrow(expr_imputed))
writeLines(qa_lines, file.path(OUT_STD, "process_QA.txt"))
cat(paste(qa_lines, collapse = "\n"), "\n")
cat("标准数据已输出至:", OUT_STD, "\n")
