#!/usr/bin/env Rscript
# 竖直棒棒糖图：GSEA 表型通路 NES 可视化（x 轴=通路，y 轴=-NES，杆竖直、顶端圆点）
# 用法: Rscript lollipop_plot.r <input.xls> <output_prefix> [top_n_per_class]
# 输入格式（无表头 TSV，同 GSEA barplot_input）：
#   V1=表型分类  V2=GOBP 通路名  V6=NES

suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
})
library(showtext)

# 注册 Arial（ggplot theme 中显式指定 family="Arial"）
font_add("Arial",
         regular    = path.expand("~/fonts/Arial/ARIAL.TTF"),
         bold       = path.expand("~/fonts/Arial/ARIALBD.TTF"),
         italic     = path.expand("~/fonts/Arial/ARIALI.TTF"),
         bolditalic = path.expand("~/fonts/Arial/ARIALBI.TTF"))
showtext_auto()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("用法: lollipop_plot.r <input.xls> <output_prefix> [top_n_per_class]")
}
infile <- args[1]
outprefix <- args[2]
top_n <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_

data <- read.table(infile, sep = "\t", header = FALSE, check.names = FALSE,
                   quote = "\"", fill = TRUE)
if (nrow(data) == 0) stop("输入文件为空: ", infile)
if (!all(c("V1", "V2", "V6") %in% names(data))) stop("输入列不足，需要 V1/V2/V6")

# 可选：每个表型分类仅保留 |NES| 最大的前 N 条通路
if (!is.na(top_n)) {
  data <- do.call(rbind, lapply(split(data, data$V1), function(d) {
    d <- d[order(-abs(d$V6)), , drop = FALSE]
    d[seq_len(min(top_n, nrow(d))), , drop = FALSE]
  }))
  rownames(data) <- NULL
}

# 通路名格式化：去 GOBP_ 前缀、下划线转空格、每词首字母大写
format_term <- function(x) {
  x <- sub("^GOBP_", "", x)
  x <- gsub("_", " ", x)
  vapply(strsplit(x, " ", fixed = TRUE), function(w) {
    paste(vapply(w, function(s) {
      paste0(toupper(substr(s, 1, 1)), tolower(substr(s, 2, nchar(s))))
    }, character(1)), collapse = " ")
  }, character(1))
}
data$term <- format_term(as.character(data$V2))

# 分类显示名与配色（呼应 RdBu 浅色热图：纤维化砖红 / 蛋白加工橙 / 心肌肥厚蓝）
class_levels <- unique(as.character(data$V1))
class_palette <- c(
  "fibrosis"            = "#B2182B",
  "protein_processing"  = "#E6550D",
  "cardiac_hypertrophy" = "#2166AC"
)
unknown <- setdiff(class_levels, names(class_palette))
if (length(unknown) > 0) stop("出现未定义配色的分类: ", paste(unknown, collapse = ", "))
class_labels <- c(
  "fibrosis"            = "Fibrosis",
  "protein_processing"  = "Protein Processing",
  "cardiac_hypertrophy" = "Cardiac Hypertrophy"
)
data$class <- factor(data$V1, levels = class_levels)

# x 轴排序：按分类分组（保持文件出现顺序），组内按 -NES 从高到低
data <- data[order(data$class, data$V6), ]
data$term <- factor(data$term, levels = data$term)

p <- ggplot(data, aes(x = term, y = -V6, colour = class)) +
  geom_segment(aes(xend = term, yend = 0), linewidth = 0.6, show.legend = FALSE) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = class_palette, labels = class_labels, name = NULL) +
  # y 轴：NES=0 线直接贴住上边框（顶部余量为 0）；底部留 0.25 单位防止底端圆点被裁剪
  scale_y_continuous(expand = expansion(add = c(0.25, 0))) +
  # 顶部圆点骑在上边框上绘制（不裁剪越界部分）
  coord_cartesian(clip = "off") +
  labs(y = "GSEA -NES", x = NULL, colour = NULL) +
  theme_bw() +
  theme(
    text             = element_text(family = "Arial", colour = "black"),
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                    size = 7, colour = "black"),
    axis.text.y      = element_text(size = 7, colour = "black"),
    axis.title.y     = element_text(size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "black", linewidth = 0.3, fill = NA),
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 7),
    legend.key       = element_blank(),
    legend.spacing.x = unit(2, "pt"),
    legend.key.width = unit(8, "pt"),
    # 图例贴近 panel 顶部边框，仅留最小间距
    legend.box.margin = margin(0, 0, 2, 0, "pt"),
    legend.box.spacing = unit(2, "pt"),
    plot.margin      = margin(4, 4, 4, 4, "pt")
  )

# 图形尺寸（紧凑版，通路名单行 90 度竖排）：
# 间距取 0.13 in/个（点径 7.4pt + 最小间隙）保证点不粘连；
# 高度 = 顶部图例/边距 0.35 + 图体 ~0.8 + 单行标签区（7pt 字符宽约 0.05 in）
n_term <- nrow(data)
max_nchar <- max(nchar(as.character(data$term)))
fig_w <- max(3.6, 0.13 * n_term + 1.0)
fig_h <- 1.15 + 0.05 * max_nchar

cairo_pdf(paste0(outprefix, ".pdf"), width = fig_w, height = fig_h)
print(p); dev.off()

showtext_opts(dpi = 600)
agg_png(paste0(outprefix, ".png"), width = fig_w, height = fig_h,
        units = "in", res = 600, background = "white")
print(p); dev.off()

cat("完成:", outprefix, ".pdf/.png |", n_term, "个通路 |",
    "尺寸", round(fig_w, 2), "x", round(fig_h, 2), "英寸\n")
