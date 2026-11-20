# 改编自参考流程 cluster_analysis_v2.r
# (原始: /media/nfs/nfs02/chenyuanbao/project/2025.07/HTP_ER_PRF_2503_0069_CKO_merge/IMA94_merge/phenotype_filter_CKO_5/filter_CKO_3/filter_WT5610/cluster_analysis_v2.r)
# 相对原始代码的改动:
#   1) 输入直接使用 2.standard_data/{region}_protein_expression_log2_imputed.xls
#      (log2 尺度, 前 4 列为注释列), 不再另建 all.genes.fpkm.xls
#   2) 蛋白以 Protein.Group 自身标识 (重复 accession 按均值合并), 不做基因符号转换
#   3) 分组匹配改为前缀 "^{Group}_" (原 \b 正则不适配 "Normal_LVPW-1" 命名:
#      下划线为词字符, \bNormal\b 无法匹配)
#   4) 输入已是 log2 尺度, 去掉 log2(x+1) 变换
#   5) 输出平铺到当前目录 (即部位目录), 去掉 01_cluster/1.QC 嵌套
#   6) PCA/t-SNE 主图面板正方形 (aspect.ratio=1); 标签 max.overlaps=Inf + 更强斥力防重叠
#   7) 树状图: 画布宽度公式与参考脚本一致 (n*40/100 in), 标签字号恢复参考值 0.7
#      (dendextend 竖排标签自动适应槽位, 18 样本 19 字符名验证无遮挡不裁切)
#   8) PCA/t-SNE 带标签版画布放大至 8x6 in (样本名长, 4x3 面板太小导致标签遮挡);
#      无标签版保持参考尺寸 4x3 in
library(gmodels)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(ggrepel)
library(Rtsne)
library(dendextend)
library(ape)

# --------------------------
# 动态识别分组数量：从args第7个参数开始，所有后续参数均为目标分组
# --------------------------
args <- commandArgs(TRUE)
input_file <- args[1]  # 第1个参数: 输入 log2+填充矩阵
target_groups <- args[2:length(args)]  # 第2个及以后参数: 目标分组（数量动态）

if (length(target_groups) < 2) {
  stop("请至少输入2个目标分组（从第7个参数开始）")
}

# 导入蛋白表达矩阵 (log2 尺度, 前 4 列注释: Protein.Group/Genes/Protein.Names/First.Protein.Description)
expr <- read.delim(input_file, header = TRUE, sep = "\t", check.names = FALSE)
id <- expr$Protein.Group
expr <- expr[, -(1:4), drop = FALSE]
if (any(duplicated(id))) {
  # 蛋白组 accession 重复时按均值合并 (保留参考流程的去重行为)
  da <- cbind(id = id, expr)
  da <- aggregate(. ~ id, data = da, FUN = mean)
  rownames(da) <- da$id
  expr <- da[, -1, drop = FALSE]
} else {
  rownames(expr) <- id
}

# --------------------------
# 动态筛选样本：仅保留目标分组的样本 (前缀匹配 "{Group}_{Region}-{rep}")
# --------------------------
keep_cols <- c()
for (g in target_groups) {
  cols <- grep(paste0("^", g, "_"), colnames(expr))
  keep_cols <- c(keep_cols, cols)
}
keep_cols <- unique(keep_cols)

if (length(keep_cols) == 0) {
  stop("未找到匹配目标分组的样本，请检查分组名是否正确（区分大小写）")
}
expr <- expr[, keep_cols, drop = FALSE]

# 生成分组标签
group <- character(ncol(expr))
for (g in target_groups) {
  group[grep(paste0("^", g, "_"), colnames(expr))] <- g
}
sorted <- target_groups

# 低表达过滤 (输入为 log2 尺度: 均值>0.5 且无负值; 填充矩阵实际全部通过, 保留参考流程逻辑)
keep_mean <- apply(expr, 1, mean) > 0.5
keep_no_zero <- apply(expr, 1, function(row) !any(row < 0))
expr <- expr[keep_mean & keep_no_zero, , drop = FALSE]
# 输入已是 log2 尺度, 不再做 log2(x+1)
data <- t(as.matrix(expr))

###########################################################
# PCA 分析及绘图（动态分组; 正方形主图 + 标签防重叠）
###########################################################
data.pca <- prcomp(data, retx = TRUE, scale = FALSE, center = TRUE)

a <- summary(data.pca)
tmp <- a[["importance"]]
pro1 <- as.numeric(sprintf("%.3f", tmp[2, 1])) * 100
pro2 <- as.numeric(sprintf("%.3f", tmp[2, 2])) * 100

pc <- as.data.frame(a[["x"]])
pc$group <- group
pc$names <- rownames(pc)

xlab <- paste("PC1(", pro1, "%)", sep = "")
ylab <- paste("PC2(", pro2, "%)", sep = "")

# 正方形主图面板: aspect.ratio=1 强制面板框为正方形 (coord_fixed 只锁定数据单位比例, 面板不一定方)
theme_font <- theme(
  aspect.ratio = 1,
  axis.text = element_text(size = 10, color = "black", face = "bold"),
  axis.title = element_text(size = 12, color = "black", face = "bold"),
  legend.title = element_text(size = 12, color = "black", face = "bold"),
  legend.text = element_text(size = 10, color = "black", face = "bold")
)

# ggrepel 标签间距 (单位: 行): 默认 0.5 适配 18 样本全量布局;
# 样本数少的子集 (如每组 3 个的 n3) 点更稀疏、个别近重合点 (IVS t-SNE 的 HY-2/3) 需更大间距,
# 通过环境变量 REPEL_BP 覆盖 (如 02.cluster_n3.sh 中 REPEL_BP=1.5)
repel_bp <- as.numeric(Sys.getenv("REPEL_BP", unset = "0.5"))

# 带标签的PCA图 (aspect.ratio=1: 面板正方形; 标签不丢弃 + 强斥力防重叠)
set.seed(1234)
pca_with_label <- ggplot(pc, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  stat_ellipse(geom = "polygon", alpha = 1/10, aes(fill = group, colour = group)) +
  geom_point(size = 2, pch = 21, colour = "black", aes(fill = group)) +
  labs(x = xlab, y = ylab) +
  theme_bw() +
  theme_font +
  labs(color = "Group", fill = "Group") +
  geom_rug(aes(colour = group)) +
  scale_color_aaas(limits = sorted) +
  scale_fill_aaas(limits = sorted) +
  geom_text_repel(aes(label = names), size = 2.2, force = 40, max.overlaps = Inf,
                  box.padding = repel_bp, point.padding = 0.2, segment.size = 0.2,
                  min.segment.length = 0, seed = 7)

# 无标签的PCA图
pca_no_label <- ggplot(pc, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  stat_ellipse(geom = "polygon", alpha = 1/10, aes(fill = group, colour = group)) +
  geom_point(size = 2, pch = 21, colour = "black", aes(fill = group)) +
  labs(x = xlab, y = ylab) +
  theme_bw() +
  theme_font +
  labs(color = "Group", fill = "Group") +
  geom_rug(aes(colour = group)) +
  scale_color_aaas(limits = sorted) +
  scale_fill_aaas(limits = sorted)

out1 <- paste(sorted, collapse = "-vs-")
out2 <- paste("PCA.", out1, sep = "")
out2_no_label <- paste("PCA.", out1, "_no_label", sep = "")

# 保存PCA图 (带标签版 6x4.5 in: 长样本名需要更大面板避免遮挡; 无标签版保持参考尺寸 4x3)
ggsave(paste0(out2, ".pdf"), pca_with_label, width = 6, height = 4.5)
ggsave(paste0(out2_no_label, ".pdf"), pca_no_label, width = 4, height = 3)

###########################################################
# t-SNE 分析及绘图（aspect.ratio=1 面板正方形 + 标签防重叠）
###########################################################
set.seed(12345)
data_unique <- unique(data)

tsne_out <- Rtsne(data_unique, perplexity = 1)

dm <- as.data.frame(tsne_out$Y)
dm$group <- group
dm$names <- rownames(data)
colnames(dm)[1:2] <- c("Dimension 1", "Dimension 2")

tsne_with_label <- ggplot(dm, aes(`Dimension 1`, `Dimension 2`)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_point(size = 2, aes(colour = group, fill = group)) +
  labs(x = "Dimension 1", y = "Dimension 2") +
  theme_bw() +
  theme_font +
  labs(color = "Group", fill = "Group") +
  stat_ellipse(geom = "polygon", alpha = 1/10, linetype = "dashed", aes(fill = group, colour = group)) +
  geom_rug(aes(color = group)) +
  scale_color_nejm(limits = sorted) +
  scale_fill_nejm(limits = sorted) +
  geom_text_repel(aes(label = names), size = 2.2, force = 40, max.overlaps = Inf,
                  box.padding = repel_bp, point.padding = 0.2, segment.size = 0.2,
                  min.segment.length = 0, seed = 7)

tsne_no_label <- ggplot(dm, aes(`Dimension 1`, `Dimension 2`)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_point(size = 2, aes(colour = group, fill = group)) +
  labs(x = "Dimension 1", y = "Dimension 2") +
  theme_bw() +
  theme_font +
  labs(color = "Group", fill = "Group") +
  stat_ellipse(geom = "polygon", alpha = 1/10, linetype = "dashed", aes(fill = group, colour = group)) +
  geom_rug(aes(color = group)) +
  scale_color_nejm(limits = sorted) +
  scale_fill_nejm(limits = sorted)

out3 <- paste("t-sne.", out1, sep = "")
out3_no_label <- paste("t-sne.", out1, "_no_label", sep = "")

# 保存t-SNE图 (带标签版 6x4.5 in, 同 PCA; 无标签版保持参考尺寸 4x3)
ggsave(paste0(out3, ".pdf"), tsne_with_label, width = 6, height = 4.5)
ggsave(paste0(out3_no_label, ".pdf"), tsne_no_label, width = 4, height = 3)

###########################################################
# 层次聚类分析及绘图（动态分组; 尺寸与字号严格同参考脚本, 标签竖排自适应）
###########################################################
method <- c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")
for(i in 1:length(method)){
  hc <- hclust(dist(data),method=method[i])
  labels <- names(cutree(hc,1))
  color <- rep(0,length(labels))
  for(m in seq_along(target_groups)) {
    # 前缀匹配, 避免 "Model" 误匹配 "HY0094_Model_..."
    colnum <- grep(paste0("^", target_groups[m], "_"), labels)
    for(n in 1:length(colnum)) {
      color[colnum[n]] <- rainbow(length(target_groups)+1)[m]
    }
  }
  out4 <- paste("hclust_",method[i],".",out1,sep="")
  color <- rep(0,length(labels(hc)))
  for(m in seq_along(target_groups)) {
    colnum <- grep(paste0("^", target_groups[m], "_"), labels(hc))
    for(n in 1:length(colnum)) {
      color[colnum[n]] <- pal_aaas("default")(length(target_groups)+1)[m]
    }
  }
  dend <- as.dendrogram(hc) %>%
    set("labels_col",color) %>%
    set("leaves_pch",17) %>%
    set("labels_cex",0.7) %>%
    set("leaves_cex",0.7) %>%
    set("leaves_col",color)
  out4 <- paste("hclust_",method[i],".",out1,sep="")
  # 画布宽度公式与参考脚本严格一致: n*40/100 in
  # 高度 3.3 in (参考为 3): 19 字符竖排标签在 cex 0.7 下底端达 ~206.5pt,
  # 高度 3.0 in (216pt) 时最后一字符被绘图区域边界裁掉 (LVPW 丢失首字母 H), 3.3 in 验证完整
  rw <- length(labels(hc))*40
  pdf(file=paste0(out4,".pdf"),width=rw/100,height=3.3)
  par(mar=c(6,5,2,1),oma=c(1,1,1,1))
  plot(dend,ylab="Height",cex.axis=1.2,cex.lab=1.2,font.axis=2,font.lab = 2)
  dev.off()
}

cat("分析完成！分组：", paste(target_groups, collapse = ", "), "\n")
cat("结果保存至：", getwd(), "\n", sep = "")
