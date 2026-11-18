#!/usr/bin/env Rscript
# Dataset A 3-omics integration.
# Centerpiece: REAL matched MOFA on Proteomics + Metabolomics (n=10, 5 vs 5).
# Added on top: RNA (3v3, separate prep) folded in via permuted pairing (prot&metab locked).
# No matched ground truth for Dataset A -> we report the integrated factor + cross-omics concordance,
# and check the permuted prot/metab loadings against the REAL matched MOFA (B-real).
## ---- locate repo root + load machine config (R/config_local.R, gitignored) ----
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT <- if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..")) else
  if (basename(getwd()) == "R") dirname(getwd()) else getwd()
if (!file.exists(file.path(ROOT, "R", "config_local.R")))
  stop("R/config_local.R not found: copy R/config_example.R to R/config_local.R and edit it.")
source(file.path(ROOT, "R", "config_local.R"))
suppressWarnings(suppressMessages(source(file.path(ROOT, "R", "lib_mofa.R"))))
set.seed(1)
LD <- A_MO_DATA_DIR; RES <- RES_DIR; FIGS <- FIGS_DIR
rd <- function(f) as.matrix(read.delim(file.path(LD, f), check.names = FALSE, row.names = 1))

## ---- load & feature-select ----
RNA <- rd("rna_logCPM_3v3.tsv")                       # 21992 x 6 (3 vs 3)
PROT<- rd("prot_expr.tsv")                            # 6314 x 10 (5 vs 5), has NAs
METAB<- rd("metab_expr.tsv")                          # 310 x 10
hv <- function(m, n) { v <- apply(m, 1, var, na.rm = TRUE); m[order(-v)[seq_len(min(n, sum(v > 0, na.rm = TRUE)))], , drop = FALSE] }
# Balanced feature selection: cap each view at ~500 top-variance features so the
# large RNA space cannot dominate the shared factor (metab has only ~310, its ceiling).
# Combined with scale_views=TRUE in fit_mofa_multi(), all three views contribute comparably.
NFEAT <- 500
RNA <- hv(RNA, NFEAT)
keepP <- rowSums(!is.na(PROT)) >= 7 & apply(PROT, 1, function(x) var(x, na.rm = TRUE) > 0)
PROT <- hv(PROT[keepP, , drop = FALSE], NFEAT)
METAB <- METAB[apply(METAB, 1, function(x) var(x, na.rm = TRUE) > 0), , drop = FALSE]

## symbol maps
rmap <- read.delim(file.path(LD, "rna_gene_map.tsv")); rsym <- setNames(rmap$symbol, rmap$gene_id)
pmap <- read.delim(file.path(LD, "prot_gene_map.tsv")); psym <- setNames(trimws(sub("\\|.*","",pmap$symbol)), pmap$FeatureID)

G1 <- A_MO_GROUP1; G2 <- A_MO_GROUP2
grp <- function(s) ifelse(grepl(A_MO_G2_REGEX, s), G2, G1)
rna_s <- colnames(RNA); pm_s <- colnames(PROT)                       # prot & metab share sample names

## ---- CENTERPIECE: real matched MOFA on prot + metab (n=10) ----
cat("Fitting real matched Proteomics+Metabolomics MOFA (n=10)...\n")
ref_prot <- { cm <- sapply(split(pm_s, grp(pm_s)), function(s) rowMeans(PROT[, s, drop=FALSE], na.rm=TRUE)); cm[,G1]-cm[,G2] }
fitR <- fit_mofa_multi(list(Proteomics = PROT, Metabolomics = METAB), K = 4, seed = 1L)
sfR  <- shared_factor_multi(fitR, c("Proteomics","Metabolomics"), anchor = "Proteomics", ref = ref_prot)
real_prot <- sfR$Proteomics; real_metab <- sfR$Metabolomics

## ---- RNA folded in: permuted 3-view (prot & metab locked) ----
# inject 3 RNA -> 3 of 5 prot/metab per group; sample NP injections
inj <- function() do.call(rbind, lapply(c(G2,G1), function(g) {
  r <- rna_s[grp(rna_s) == g]; pm <- sample(pm_s[grp(pm_s) == g], 3)
  data.frame(rna = r, pm = pm, stringsAsFactors = FALSE)
}))
NP <- 150
Wr <- matrix(NA_real_, nrow(RNA), NP, dimnames=list(rownames(RNA), NULL))
Wp <- matrix(NA_real_, nrow(PROT), NP, dimnames=list(rownames(PROT), NULL))
Wm <- matrix(NA_real_, nrow(METAB), NP, dimnames=list(rownames(METAB), NULL))
for (i in seq_len(NP)) {
  pr <- inj(); pid <- paste0("s", seq_len(nrow(pr)))
  Rv <- RNA[, pr$rna];  colnames(Rv) <- pid
  Pv <- PROT[, pr$pm];  colnames(Pv) <- pid
  Mv <- METAB[, pr$pm]; colnames(Mv) <- pid
  fit <- fit_mofa_multi(list(RNA=Rv, Proteomics=Pv, Metabolomics=Mv), K=3, seed=1L)
  sf  <- shared_factor_multi(fit, c("RNA","Proteomics","Metabolomics"), anchor="Proteomics", ref=ref_prot)
  Wr[,i]<-sf$RNA; Wp[,i]<-sf$Proteomics; Wm[,i]<-sf$Metabolomics
  cat(sprintf("\r  RNA-folded fit %d/%d", i, NP)); flush.console()
}
cat("\n")
rec_rna<-apply(Wr,1,median); rec_prot<-apply(Wp,1,median); rec_metab<-apply(Wm,1,median)

## consistency: 3-view prot/metab vs real matched MOFA
consist <- c(prot = abs(cor(rec_prot, real_prot[names(rec_prot)], use="complete.obs")),
             metab= abs(cor(rec_metab, real_metab[names(rec_metab)], use="complete.obs")))

## ---- integrated feature tables (top loaders per view) ----
top_tab <- function(load, symmap=NULL, k=25) {
  o <- order(-abs(load)); d <- data.frame(feature=names(load)[o][1:k], loading=round(load[o][1:k],3))
  if (!is.null(symmap)) d$name <- symmap[d$feature]
  d
}
saveRDS(list(rec=list(RNA=rec_rna,Proteomics=rec_prot,Metabolomics=rec_metab),
             real=list(Proteomics=real_prot,Metabolomics=real_metab),
             consist=consist, Wr=Wr,Wp=Wp,Wm=Wm), file.path(RES,paste0(OUT_PREFIX_A,"_3view.rds")))
write.csv(top_tab(rec_rna,rsym), file.path(RES,paste0(OUT_PREFIX_A,"_MO_top_RNA.csv")), row.names=FALSE)
write.csv(top_tab(rec_prot,psym), file.path(RES,paste0(OUT_PREFIX_A,"_MO_top_protein.csv")), row.names=FALSE)
write.csv(top_tab(rec_metab), file.path(RES,paste0(OUT_PREFIX_A,"_MO_top_metabolite.csv")), row.names=FALSE)

sink(file.path(RES,paste0(OUT_PREFIX_A,"_3view_summary.txt")))
cat("=== Dataset A 3-omics integration ===\n\n")
cat(sprintf("RNA %d(3v3) | Prot %d | Metab %d feats. Centerpiece=matched prot+metab MOFA; RNA folded via %d permuted injections.\n\n", nrow(RNA),nrow(PROT),nrow(METAB),NP))
cat(sprintf("3-view permuted prot/metab loadings vs REAL matched prot+metab MOFA: prot |r|=%.3f, metab |r|=%.3f\n",consist[1],consist[2]))
cat(sprintf("  (high = RNA-folding did not distort the real matched signal)\n\n"))
cat("Top RNA genes on the shared group factor:\n"); print(head(top_tab(rec_rna,rsym),12),row.names=FALSE)
cat("\nTop proteins:\n"); print(head(top_tab(rec_prot,psym),12),row.names=FALSE)
cat("\nTop metabolites:\n"); print(head(top_tab(rec_metab),12),row.names=FALSE)
sink()
cat(readLines(file.path(RES,paste0(OUT_PREFIX_A,"_3view_summary.txt"))), sep="\n")
