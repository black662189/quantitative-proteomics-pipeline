args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(BiocParallel))
cfg <- fromJSON(args[1],simplifyVector=FALSE)
out <- args[2]
dest <- file.path(out,"5.Pathways")
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
write_tsv <- function(x,p) write.table(x,p,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
read_gmt <- function(p) {
  a <- strsplit(readLines(p,warn=FALSE),"\t",fixed=TRUE)
  ans <- list()
  for (v in a) if(length(v)>=3) ans[[v[1]]] <- unique(c(ans[[v[1]]],v[-c(1,2)]))
  ans
}
meta <- read.delim(file.path(out,"0.InputQC/samples.tsv"),check.names=FALSE)
meta <- meta[meta$role=="sample",,drop=FALSE]
x <- as.matrix(read.delim(file.path(out,"1.Preprocess/normalized_log2.tsv"),row.names=1,check.names=FALSE))
ann <- read.delim(file.path(out,"0.InputQC/annotation.tsv"),check.names=FALSE,na.strings="")
ann <- ann[match(rownames(x),ann$protein_id),]
valid <- !is.na(ann$gene) & nzchar(ann$gene)
x <- x[valid,meta$sample,drop=FALSE]; genes <- ann$gene[valid]
# Collapse only for pathway analyses; protein-level statistics remain unchanged.
groups <- split(seq_along(genes),genes)
gx <- t(vapply(groups,function(i) apply(x[i,,drop=FALSE],2,function(v) {
  if(all(is.na(v))) NA_real_ else median(v,na.rm=TRUE)
}),numeric(ncol(x))))
colnames(gx) <- colnames(x)
complete <- rowSums(!is.finite(gx))==0
variable <- apply(gx,1,function(v) if(all(is.finite(v))) sd(v)>0 else FALSE)
write_tsv(data.frame(gene=rownames(gx),complete=complete,variable=variable),
          file.path(dest,"GSVA_gene_eligibility.tsv"))
gx <- gx[complete & variable,,drop=FALSE]
write_tsv(data.frame(gene=rownames(gx),gx,check.names=FALSE),file.path(dest,"GSVA_gene_matrix.tsv"))
statuses <- list()
for (gs in cfg$gene_sets) {
  sets <- read_gmt(gs$path)
  minsize <- if(is.null(gs$min_size)) 5L else as.integer(gs$min_size)
  maxsize <- if(is.null(gs$max_size)) 500L else as.integer(gs$max_size)
  if(isTRUE(cfg$gsea)) for(co in cfg$comparisons) {
    r <- read.delim(file.path(out,"2.Comparison",paste0(co$name,".all_proteins.tsv")),check.names=FALSE)
    r <- r[tolower(as.character(r$tested))=="true" & !is.na(r$gene) & nzchar(r$gene) & is.finite(r$t),,drop=FALSE]
    ranks <- tapply(r$t,r$gene,median)
    ranks <- sort(ranks,decreasing=TRUE)
    set.seed(20261020)
    eligible <- sets[lengths(lapply(sets,intersect,names(ranks)))>=minsize &
                     lengths(lapply(sets,intersect,names(ranks)))<=maxsize]
    if(length(eligible)>0 && length(ranks)>1 && any(ranks>0) && any(ranks<0)) {
      res <- as.data.frame(fgseaMultilevel(pathways=eligible,stats=ranks,minSize=minsize,
                          maxSize=maxsize,eps=0,nPermSimple=1000,BPPARAM=SerialParam()))
      res$leadingEdge <- vapply(res$leadingEdge,paste,collapse=";",FUN.VALUE="")
      res <- res[order(res$padj),]
      write_tsv(res,file.path(dest,paste0(co$name,".",gs$name,".GSEA.tsv")))
      status <- "DONE"
    } else {res <- data.frame();status <- "NO_TEST: insufficient eligible sets or two-sided ranks"}
    write_tsv(data.frame(gene=names(ranks),statistic=as.numeric(ranks)),
              file.path(dest,paste0(co$name,".",gs$name,".GSEA_ranking.tsv")))
    statuses[[length(statuses)+1]] <- list(module="GSEA",database=gs$name,comparison=co$name,
       genes=length(ranks),tested_terms=nrow(res),status=status)
  }
  if(isTRUE(cfg$gsva)) {
    eligible <- lapply(sets,intersect,rownames(gx))
    eligible <- eligible[lengths(eligible)>=minsize & lengths(eligible)<=maxsize]
    if(length(eligible)>0 && nrow(gx)>1) {
      param <- gsvaParam(gx,eligible,kcdf="Gaussian",minSize=minsize,maxSize=maxsize)
      scores <- gsva(param,verbose=FALSE,BPPARAM=SerialParam())
      write_tsv(data.frame(pathway=rownames(scores),scores,check.names=FALSE),
                file.path(dest,paste0(gs$name,".GSVA_scores.tsv")))
      for(co in cfg$comparisons) {
        m <- meta[meta$group %in% c(co$control,co$treatment),,drop=FALSE]
        y <- scores[,m$sample,drop=FALSE]
        if(isTRUE(co$paired)) {
          subjects <- sort(unique(m$subject))
          y <- do.call(cbind,lapply(subjects,function(s)
            y[,m$subject==s & m$group==co$treatment] - y[,m$subject==s & m$group==co$control]))
          rownames(y) <- rownames(scores)
          design <- matrix(1,nrow=length(subjects),ncol=1)
          coef <- 1
        } else {
          m$condition <- factor(m$group,levels=c(co$control,co$treatment))
          if(length(unique(m$plex))>1) {
            m$batch <- factor(m$plex);design <- model.matrix(~condition+batch,m)
          } else design <- model.matrix(~condition,m)
          coef <- 2
        }
        if(qr(design)$rank<ncol(design)) stop("GSVA design confounded")
        if(nrow(y)>=2) {
          fit <- eBayes(lmFit(y,design),trend=FALSE)
          tab <- data.frame(pathway=rownames(y),score_difference=fit$coefficients[,coef],
                  t=fit$t[,coef],P.Value=fit$p.value[,coef],
                  adj.P.Val=p.adjust(fit$p.value[,coef],"BH"))
          write_tsv(tab[order(tab$adj.P.Val),],file.path(dest,paste0(co$name,".",gs$name,".GSVA_difference.tsv")))
        }
      }
      status <- "DONE"
    } else status <- "NO_TEST: insufficient complete variable genes or eligible pathways"
    statuses[[length(statuses)+1]] <- list(module="GSVA",database=gs$name,genes=nrow(gx),
      tested_terms=length(eligible),status=status)
  }
}
write_json(statuses,file.path(dest,"pathway_status.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(dest,"R_sessionInfo.txt"))
