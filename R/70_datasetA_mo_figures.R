#!/usr/bin/env Rscript
# Clean tables + figures for the Dataset A 3-omics integration.
## ---- locate repo root + load machine config (R/config_local.R, gitignored) ----
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT <- if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..")) else
  if (basename(getwd()) == "R") dirname(getwd()) else getwd()
if (!file.exists(file.path(ROOT, "R", "config_local.R")))
  stop("R/config_local.R not found: copy R/config_example.R to R/config_local.R and edit it.")
source(file.path(ROOT, "R", "config_local.R"))
suppressWarnings(suppressMessages({ library(MOFA2) }))
LD<-A_MO_DATA_DIR; RES<-RES_DIR; FIGS<-FIGS_DIR
x <- readRDS(file.path(RES,paste0(OUT_PREFIX_A,"_3view.rds")))

## maps
rmap<-read.delim(file.path(LD,"rna_gene_map.tsv")); rsym<-setNames(rmap$symbol,rmap$gene_id)
pmap<-read.delim(file.path(LD,"prot_gene_map.tsv")); psym<-setNames(trimws(sub("\\|.*","",pmap$symbol)),pmap$FeatureID)
an<-read.delim(file.path(LD,"metab_annot.tsv"),check.names=FALSE)
# find which annot col matches expr rownames (RT.._MZ..)
exr<-names(x$rec$Metabolomics)
idcol<-names(which.max(sapply(an, function(c) mean(exr %in% as.character(c)))))
msym<-setNames(ifelse(is.na(an$Name)|an$Name=="", as.character(an[[idcol]]), an$Name), as.character(an[[idcol]]))

lab<-function(v,map,k=15){ o<-order(-abs(v))[1:k]; nm<-if(is.null(map)) names(v)[o] else {n<-map[names(v)[o]]; ifelse(is.na(n)|n=="",names(v)[o],n)}; data.frame(name=nm, loading=v[o]) }
TR<-lab(x$rec$RNA,rsym); TP<-lab(x$rec$Proteomics,psym); TM<-lab(x$rec$Metabolomics,msym)
write.csv(TR,file.path(RES,paste0(OUT_PREFIX_A,"_MO_top_RNA.csv")),row.names=FALSE)
write.csv(TP,file.path(RES,paste0(OUT_PREFIX_A,"_MO_top_protein.csv")),row.names=FALSE)
write.csv(TM,file.path(RES,paste0(OUT_PREFIX_A,"_MO_top_metabolite.csv")),row.names=FALSE)

## ---- Figure: top loadings per omic (diverging bars) ----
NAVY<-"#243b53"; TEAL<-"#1FB8C4"; CORAL<-"#E8735A"; GOLD<-"#f0a202"
bar<-function(d,title,accent){
  d<-d[order(d$loading),]; cols<-ifelse(d$loading>0,TEAL,CORAL)
  par(mar=c(4,8.5,3,1))
  bp<-barplot(d$loading,horiz=TRUE,col=cols,border=NA,names.arg=d$name,las=1,cex.names=.8,
              xlab="shared-factor loading",main=title,col.main=NAVY,font.main=2)
  abline(v=0,col="grey60")
}
png(file.path(FIGS,paste0(OUT_PREFIX_A,"_mo_loadings.png")),width=1650,height=650,res=150)
par(mfrow=c(1,3))
bar(lab(x$rec$RNA,rsym,12),"RNA (top genes)",TEAL)
bar(lab(x$rec$Proteomics,psym,12),"Proteins (top)",TEAL)
bar(lab(x$rec$Metabolomics,msym,12),"Metabolites (top)",TEAL)
dev.off()

## ---- Figure: real matched prot+metab MOFA variance explained ----
png(file.path(FIGS,paste0(OUT_PREFIX_A,"_mo_consistency.png")),width=900,height=520,res=150)
par(mar=c(4.2,4.2,3,1))
plot(x$real$Proteomics[names(x$rec$Proteomics)], x$rec$Proteomics, pch=16,col="#243b5399",cex=.5,
     xlab="REAL matched prot+metab MOFA (protein loading)", ylab="3-view (RNA folded in)",
     main=sprintf("Folding RNA in preserves the matched signal\nprotein |r|=%.3f",x$consist[1]))
abline(0,coef(lm(x$rec$Proteomics~x$real$Proteomics[names(x$rec$Proteomics)]))[2],col=GOLD,lwd=2)
dev.off()

cat("figures + tables written\n")
cat("\nTop proteins on shared group factor:\n"); print(TP[order(-abs(TP$loading)),][1:12,],row.names=FALSE)
cat("\nTop metabolites (named where available):\n"); print(TM[order(-abs(TM$loading)),][1:10,],row.names=FALSE)
cat(sprintf("\nmetab id column matched: %s | named top-15: %d/15\n", idcol, sum(!grepl("^RT",TM$name[order(-abs(TM$loading))][1:15]))))
