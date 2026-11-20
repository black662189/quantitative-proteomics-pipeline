# GSVA 打分 (GSVA2 引擎, maxDiff)
# 用法: Rscript 07b_GSVA_calc.R <geneset.gmt> <gene_expr.xls> <out_score.xls> [workers]
options(stringsAsFactors = FALSE)
args <- commandArgs(TRUE)
if (length(args) < 3) stop("用法: Rscript 07b_GSVA_calc.R <gmt> <expr> <out> [workers]")

genesetgmt <- args[1]
geneexp    <- args[2]
outfile    <- args[3]
workers    <- if (length(args) >= 4) as.integer(args[4]) else 4

library(GSEABase)
library(GSVAdata)
library(GSVA)
library(limma)
library(BiocParallel)

geneSets <- getGmt(genesetgmt)
mydata <- read.table(file = geneexp, header = TRUE, check.names = FALSE, sep = "\t")

index <- duplicated(mydata[, 1])
fildup <- mydata[!index, ]
exp <- fildup[, -1]
exp <- subset(exp, select = -c(description))
name <- fildup[, 1]
rownames(exp) <- name

allmydata <- as.matrix(exp)

params <- gsvaParam(
  exprData = allmydata,
  geneSets = geneSets,
  minSize = 10,
  maxSize = 500,
  maxDiff = TRUE
)

es.diff <- gsva(params, BPPARAM = MulticoreParam(workers = workers))

all_GSVA_score <- cbind(path = rownames(es.diff), es.diff)
write.table(all_GSVA_score, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)

cat("GSVA score written:", outfile, "\n")
