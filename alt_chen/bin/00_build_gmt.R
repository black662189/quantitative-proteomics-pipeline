# ============================================================================
# 通用基因集 GMT 生成器 (human / mouse / rat)
# KEGG:   原生 (KEGGREST) -> org.<Sp>.eg.db 转 SYMBOL
# GOBP/Reactome:
#   human: msigdbr 默认 <- 人类原生 (category C5/C2)
#   mouse: msigdbr db_species="MM" <- 原生小鼠 (category M5/M2)
#   rat:   msigdbr 默认 species="Rattus norvegicus" <- 人类库同源映射到 RGNC 符号 (category C5/C2)
# 用法: PP_SPECIES=human|mouse|rat  PP_GMT_OUTPUT=<dir>  Rscript 00_build_gmt.R
# ============================================================================
options(stringsAsFactors = FALSE)

SPECIES <- tolower(Sys.getenv("PP_SPECIES", "human"))
OUT     <- Sys.getenv("PP_GMT_OUTPUT")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# db: 是否用 db_species(原生小鼠); catp: msigdbr category 前缀(C 或 M)
sp <- list(
  human = list(kegg = "hsa", org = "org.Hs.eg.db",
               ms = "Homo sapiens", db = "HS", catp = "C",
               tag = "hsa", label = "Homo sapiens"),
  mouse = list(kegg = "mmu", org = "org.Mm.eg.db",
               ms = "Mus musculus", db = "MM", catp = "M",
               tag = "mmu", label = "Mus musculus"),
  rat   = list(kegg = "rno", org = "org.Rn.eg.db",
               ms = "Rattus norvegicus", db = "HS", catp = "C",
               tag = "rno", label = "Rattus norvegicus")
)
if (! SPECIES %in% names(sp)) stop("PP_SPECIES 必须是 human/mouse/rat, 当前=", SPECIES)
S <- sp[[SPECIES]]

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(KEGGREST)
  library(msigdbr)
  library(S$org, character.only = TRUE)
})

f_kegg  <- file.path(OUT, paste0(S$tag, "_kegg.symbol.gmt"))
f_gobp  <- file.path(OUT, paste0(S$tag, "_go_bp.symbol.gmt"))
f_react <- file.path(OUT, paste0(S$tag, "_reactome.symbol.gmt"))

## 获取 msigdbr 结果 (每条含 gs_name / gene_symbol)
get_msig <- function(cat, sub) {
  if (S$catp == "M") {
    msigdbr(species = S$ms, db_species = "MM", category = cat, subcategory = sub)
  } else {
    # 人类原生(human) 或 人类库同源映射(rat): 都不传 db_species
    msigdbr(species = S$ms, category = cat, subcategory = sub)
  }
}
write_gmt <- function(fname, df) {
  df <- df[, c("gs_name", "gene_symbol")]
  sets <- split(df$gene_symbol, df$gs_name)
  con <- file(fname, "w", encoding = "UTF-8")
  n <- 0
  for (nm in names(sets)) {
    g <- unique(sets[[nm]]); g <- g[nzchar(g)]
    if (length(g) < 3) next
    writeLines(paste(c(nm, "na", g), collapse = "\t"), con, useBytes = TRUE)
    n <- n + 1
  }
  close(con)
  n
}

# ---------- 1) KEGG (原生) ----------
cat("[", SPECIES, "] 构建 KEGG ...\n")
pw_list <- keggList("pathway", S$kegg)
pw_ann  <- keggLink("pathway", S$kegg)
gene2pw <- sub("path:", "", pw_ann)
names(gene2pw) <- sub(paste0(S$kegg, ":"), "", names(pw_ann))
entrez <- unique(names(gene2pw))
sym_map <- suppressMessages(AnnotationDbi::mapIds(get(S$org),
        keys = entrez, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
names(gene2pw) <- sym_map[names(gene2pw)]
gene2pw <- gene2pw[!is.na(names(gene2pw))]

con <- file(f_kegg, "w", encoding = "UTF-8")
n <- 0
for (pw in names(pw_list)) {
  genes <- unique(names(gene2pw)[gene2pw == pw])
  if (length(genes) < 3) next
  name <- sub(paste0(" - ", S$label, ".*$"), "", pw_list[[pw]])
  writeLines(paste(c(name, "na", genes), collapse = "\t"), con, useBytes = TRUE)
  n <- n + 1
}
close(con)
cat("  KEGG 通路数:", n, "\n")

# ---------- 2) GOBP ----------
cat("[", SPECIES, "] 构建 GOBP ...\n")
gobp <- get_msig(paste0(S$catp, "5"), "GO:BP")
cat("  GOBP 通路数:", write_gmt(f_gobp, gobp), "\n")

# ---------- 3) Reactome ----------
cat("[", SPECIES, "] 构建 Reactome ...\n")
rea <- get_msig(paste0(S$catp, "2"), "CP:REACTOME")
cat("  Reactome 通路数:", write_gmt(f_react, rea), "\n")

cat("完成 ->", OUT, "\n")
