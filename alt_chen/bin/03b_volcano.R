# ============================================================================
# 差异蛋白火山图 (nature 紧凑规范)
# 输入: DE_{contrast}.xls; 阈值 pvalue<0.05; 上红 #EE0000 下蓝 #3B4992 灰 #BBBBBB
# 出图: theme_classic 6.5pt Arial 细线 紧凑; 标注 top10 |log2FC| 显著基因
# ============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(ggplot2); library(ragg) })

args <- commandArgs(TRUE)
if (length(args) != 1)
  stop("用法: Rscript 03b_volcano.R <region>")
rg <- args[1]

OUT_DE  <- Sys.getenv("PP_OUT_DE")
FCCUT   <- as.numeric(Sys.getenv("PP_FCCUT"))
PCUT    <- as.numeric(Sys.getenv("PP_PCUT"))
CONTRASTS <- strsplit(Sys.getenv("PP_CONTRASTS"), " ")[[1]]

indir <- file.path(OUT_DE, rg)
outdir <- file.path(OUT_DE, rg, "Volcano")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

theme_vol <- theme_classic(base_size = 6.5, base_family = "Arial") %+replace%
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text = element_text(colour = "black", size = 6),
    axis.title = element_text(colour = "black", size = 6.5),
    plot.title = element_text(colour = "black", size = 7, face = "bold", hjust = 0),
    plot.margin = margin(2, 2, 2, 2, "mm"),
    panel.grid = element_blank(),
    legend.position = "none"
  )

for (cc in CONTRASTS) {
  de <- read.delim(file.path(indir, sprintf("DE_%s.xls", cc)), check.names = FALSE)
  de <- de[is.finite(de$log2FoldChange) & is.finite(de$pvalue), , drop = FALSE]
  fc_cut <- FCCUT
  de$sig <- ifelse(de$pvalue >= PCUT | abs(de$log2FoldChange) < fc_cut, "n.s.",
                   ifelse(de$log2FoldChange > 0, "up", "down"))
  de$mlp <- -log10(de$pvalue)
  n_up <- sum(de$sig == "up"); n_down <- sum(de$sig == "down"); n_ns <- sum(de$sig == "n.s.")

  sig <- de[de$sig != "n.s.", , drop = FALSE]
  sig_up <- sig[sig$log2FoldChange > 0, , drop = FALSE]
  sig_dn <- sig[sig$log2FoldChange < 0, , drop = FALSE]
  top_up <- head(sig_up[order(-sig_up$log2FoldChange), ]$gene, 10)
  top_dn <- head(sig_dn[order(sig_dn$log2FoldChange), ]$gene, 10)
  de$lab <- ifelse(de$gene %in% c(top_up, top_dn), de$gene, NA)

  xlim <- max(abs(de$log2FoldChange)) * 1.08
  p <- ggplot(de, aes(log2FoldChange, mlp, colour = sig)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_colour_manual(values = c("up" = "#EE0000", "down" = "#3B4992", "n.s." = "#BBBBBB")) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70", linetype = "dashed") +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linewidth = 0.25, colour = "grey50", linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), linewidth = 0.25, colour = "grey70", linetype = "dashed") +
    ggrepel::geom_text_repel(aes(label = lab), size = 1.9, colour = "black", na.rm = TRUE,
              family = "Arial", max.overlaps = 30, min.segment.length = 0,
              segment.size = 0.2, segment.colour = "grey50",
              box.padding = 0.25, point.padding = 0.2, seed = 1) +
    labs(x = "log2 fold change", y = "-log10(pvalue)") +
    xlim(-xlim, xlim) +
    theme_vol

  lab <- c("n.s." = sprintf("n.s. (%d)", n_ns),
           "up" = sprintf("up (%d)", n_up),
           "down" = sprintf("down (%d)", n_down))
  # 计数图例置于左上角(负 log2FC 侧, 该区显著点与标签较少), 避免与右上"up"显著点/标签重叠
  ymax_plot <- max(de$mlp)
  p <- p +
    annotate("text", x = -xlim * 0.02, y = ymax_plot, hjust = 0, vjust = 1.0,
             label = lab[["n.s."]], colour = "#888888", size = 2.1, family = "Arial") +
    annotate("text", x = -xlim * 0.02, y = ymax_plot * 0.94, hjust = 0, vjust = 1.0,
             label = lab[["down"]], colour = "#3B4992", size = 2.1, family = "Arial") +
    annotate("text", x = -xlim * 0.02, y = ymax_plot * 0.88, hjust = 0, vjust = 1.0,
             label = lab[["up"]], colour = "#EE0000", size = 2.1, family = "Arial")

  prefix <- file.path(outdir, sprintf("%s_volcano", cc))
  ggsave(paste0(prefix, ".pdf"), p, width = 3.3, height = 2.8, device = cairo_pdf)
  ggsave(paste0(prefix, ".png"), p, width = 3.3, height = 2.8, dpi = 600, device = ragg::agg_png)
  cat(sprintf("%s %s: up=%d down=%d ns=%d\n", rg, cc, n_up, n_down, n_ns))
}
cat(sprintf("%s: 火山图完成\n", rg))
