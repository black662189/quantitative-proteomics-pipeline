args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(jsonlite))
cfg <- fromJSON(args[1], simplifyVector=FALSE)
x <- as.matrix(read.delim(args[2],row.names=1,check.names=FALSE))
meta <- read.delim(args[3],check.names=FALSE,stringsAsFactors=FALSE)
dir.create(args[4],recursive=TRUE,showWarnings=FALSE)
logs <- list()
for (co in cfg$comparisons) {
  m <- meta[meta$group %in% c(co$control,co$treatment) & meta$role=="sample",,drop=FALSE]
  y <- x[,m$sample,drop=FALSE]
  a <- which(m$group==co$control); b <- which(m$group==co$treatment)
  na <- rowSums(is.finite(y[,a,drop=FALSE])); nb <- rowSums(is.finite(y[,b,drop=FALSE]))
  keep <- na>=max(cfg$min_replicates,ceiling(length(a)*cfg$observed_fraction)) &
          nb>=max(cfg$min_replicates,ceiling(length(b)*cfg$observed_fraction))
  complete <- rep(NA_integer_,nrow(y))
  if (isTRUE(co$paired)) {
    subjects <- sort(unique(m$subject))
    delta <- sapply(subjects,function(s) y[,which(m$subject==s & m$group==co$treatment)] -
                                          y[,which(m$subject==s & m$group==co$control)])
    if (is.null(dim(delta))) delta <- matrix(delta,ncol=length(subjects))
    rownames(delta)<-rownames(y)
    complete <- rowSums(is.finite(delta))
    keep <- keep & complete>=max(cfg$min_replicates,ceiling(length(subjects)*cfg$observed_fraction))
    design <- matrix(1,nrow=length(subjects),ncol=1,dimnames=list(subjects,"effect"))
    yfit <- delta[keep,,drop=FALSE]
    coefficient <- 1
  } else {
    m$condition <- factor(ifelse(m$group==co$treatment,"treated","control"),levels=c("control","treated"))
    if(length(unique(m$plex))>1) {
      m$batch<-factor(m$plex)
      design<-model.matrix(~condition+batch,m)
    } else design<-model.matrix(~condition,m)
    if(qr(design)$rank<ncol(design)) stop("Confounded/rank-deficient design: ",co$name)
    if(nrow(design)-ncol(design)<1)stop("No residual degrees of freedom: ",co$name)
    coefficient <- which(colnames(design)=="conditiontreated")
    yfit <- y[keep,,drop=FALSE]
  }
  eligibility <- data.frame(protein_id=rownames(y),n_control=na,n_treatment=nb,
                            n_complete_pairs=complete,tested=keep)
  write.table(eligibility,file.path(args[4],paste0(co$name,".eligibility.tsv")),sep="\t",row.names=FALSE,quote=FALSE)
  if(sum(keep)<2)stop("Fewer than two testable proteins: ",co$name)
  fit<-lmFit(yfit,design)
  fit<-eBayes(fit,trend=TRUE,robust=FALSE)
  tab<-data.frame(protein_id=rownames(yfit),log2FC=fit$coefficients[,coefficient],
                  AveExpr=rowMeans(y[keep,,drop=FALSE],na.rm=TRUE),
                  t=fit$t[,coefficient],P.Value=fit$p.value[,coefficient],
                  adj.P.Val=p.adjust(fit$p.value[,coefficient],method="BH"))
  if(any(!is.finite(tab$P.Value)))stop("Non-finite test results: ",co$name)
  write.table(tab,file.path(args[4],paste0(co$name,".statistics.tsv")),sep="\t",row.names=FALSE,quote=FALSE)
  logs[[co$name]]<-list(tested=nrow(tab),paired=isTRUE(co$paired),design=colnames(design),
                        limma=as.character(packageVersion("limma")))
}
write_json(logs,file.path(args[4],"models.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(args[4],"R_sessionInfo.txt"))
