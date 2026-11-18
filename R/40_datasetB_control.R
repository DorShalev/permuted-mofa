#!/usr/bin/env Rscript
# POSITIVE CONTROL + SELECTION variant, on Dataset B (matched -> deliberately unmatched).
# Ground truth = matched-MOFA Factor1 loadings (RNA + Proteomics).
# Q1: does unmatched-MOFA (AVERAGE over pairings) recover the matched Factor1 features?
# Q2: does SELECT-best-pairing (Wasserstein-Procrustes style) do better/differently?
## ---- locate repo root + load machine config (R/config_local.R, gitignored) ----
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT <- if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..")) else
  if (basename(getwd()) == "R") dirname(getwd()) else getwd()
if (!file.exists(file.path(ROOT, "R", "config_local.R")))
  stop("R/config_local.R not found: copy R/config_example.R to R/config_local.R and edit it.")
source(file.path(ROOT, "R", "config_local.R"))
suppressWarnings(suppressMessages({
  source(file.path(ROOT, "R", "lib_mofa.R"))
}))
set.seed(1)
ED <- B_DATA_DIR; RES <- RES_DIR; FIGS <- FIGS_DIR

rd <- function(f) as.matrix(read.delim(file.path(ED, f), check.names = FALSE, row.names = 1))
RNA <- rd("RNA.tsv"); PROT <- rd("Proteomics.tsv")
meta <- read.delim(file.path(ED, "meta.tsv"), stringsAsFactors = FALSE)
# RNA & Proteomics measured for 2 of the 3 conditions (the third is all-NA in these views) -> clean 4v4
keep <- colnames(RNA)[colSums(is.na(RNA)) == 0 & colSums(is.na(PROT)) == 0]
RNA <- RNA[, keep, drop = FALSE]; PROT <- PROT[, keep, drop = FALSE]
meta <- meta[meta$sample %in% keep, ]
cat("Dataset B RNA/Prot usable samples:", length(keep), "\n")
gtR <- read.delim(file.path(ED, "gt_weights_RNA.tsv"), check.names = FALSE, row.names = 1)
gtP <- read.delim(file.path(ED, "gt_weights_Proteomics.tsv"), check.names = FALSE, row.names = 1)
gt_rna <- setNames(gtR$Factor1, rownames(gtR)); gt_prot <- setNames(gtP$Factor1, rownames(gtP))
cond <- setNames(meta$condition, meta$sample)

## condition-centroid direction per view (GT-independent reference & baseline)
centroid_dir <- function(mat) {
  cm <- sapply(split(colnames(mat), cond[colnames(mat)]), function(s) rowMeans(mat[, s, drop = FALSE]))
  prcomp(t(cm), scale. = FALSE)$rotation[, 1]     # PC1 of condition means = dominant condition axis
}
ref_rna <- centroid_dir(RNA); ref_prot <- centroid_dir(PROT)

## enumerate condition-respecting pairings (RNA sample -> Prot sample within condition); sample N
rand_pairing <- function() {
  do.call(rbind, lapply(unique(meta$condition), function(g) {
    a <- meta$sample[meta$condition == g]; b <- sample(meta$sample[meta$condition == g])
    data.frame(rna = a, prot = b, stringsAsFactors = FALSE)
  }))
}
NP <- 120
pairings <- replicate(NP, rand_pairing(), simplify = FALSE)

## fit each pairing, collect shared-factor loadings + balanced-R2 (selection score)
Wr <- matrix(NA_real_, nrow(RNA), NP, dimnames = list(rownames(RNA), NULL))
Wp <- matrix(NA_real_, nrow(PROT), NP, dimnames = list(rownames(PROT), NULL))
score <- numeric(NP)
for (i in seq_len(NP)) {
  pr <- pairings[[i]]; pid <- paste0("s", seq_len(nrow(pr)))
  R <- RNA[, pr$rna, drop = FALSE];  colnames(R) <- pid
  P <- PROT[, pr$prot, drop = FALSE]; colnames(P) <- pid
  fit <- fit_mofa_generic(R, P, "RNA", "Proteomics", K = 3, seed = 1L)
  sf  <- shared_factor(fit, "RNA", "Proteomics", ref = ref_rna)
  Wr[, i] <- sf$v1; Wp[, i] <- sf$v2; score[i] <- sf$balanced_r2
  cat(sprintf("\r  fit %d/%d", i, NP)); flush.console()
}
cat("\n")

## AVERAGE (consensus) and SELECT (best balanced-R2 pairing)
avg_rna <- apply(Wr, 1, median); avg_prot <- apply(Wp, 1, median)
best <- which.max(score)
sel_rna <- Wr[, best]; sel_prot <- Wp[, best]

## metrics vs ground-truth matched Factor1
topK <- function(v, k = 30, dec = TRUE) names(sort(abs(v), decreasing = TRUE))[1:k]
jac  <- function(a, b) length(intersect(a, b)) / length(union(a, b))
gt_top_rna <- topK(gt_rna); gt_top_prot <- topK(gt_prot)

metric <- function(rec_rna, rec_prot, lab) {
  data.frame(
    approach = lab,
    rna_cor_vs_GT  = cor(rec_rna[names(gt_rna)],  gt_rna,  use = "complete.obs"),
    prot_cor_vs_GT = cor(rec_prot[names(gt_prot)], gt_prot, use = "complete.obs"),
    rna_top30_overlap  = jac(topK(rec_rna),  gt_top_rna),
    prot_top30_overlap = jac(topK(rec_prot), gt_top_prot))
}
cen_rna <- ref_rna[names(gt_rna)]; cen_prot <- ref_prot[names(gt_prot)]
# align centroid sign to GT for fair correlation reporting
cen_rna <- cen_rna * sign(cor(cen_rna, gt_rna, use = "complete.obs"))
cen_prot<- cen_prot* sign(cor(cen_prot, gt_prot, use = "complete.obs"))

tab <- rbind(
  metric(avg_rna, avg_prot, "AVERAGE (consensus)"),
  metric(sel_rna, sel_prot, "SELECT (best pairing)"),
  metric(setNames(cen_rna, names(gt_rna)), setNames(cen_prot, names(gt_prot)), "condition-centroid"))
tab[,-1] <- round(tab[,-1], 3)

## also: does AVERAGE == centroid here too?
avg_vs_cen <- c(rna = cor(avg_rna[names(ref_rna)], ref_rna, use="complete.obs"),
                prot= cor(avg_prot[names(ref_prot)], ref_prot, use="complete.obs"))

saveRDS(list(Wr=Wr, Wp=Wp, score=score, tab=tab, avg_vs_cen=avg_vs_cen,
             gt_top_rna=gt_top_rna, gt_top_prot=gt_top_prot), file.path(RES, paste0(OUT_PREFIX_B, "_control.rds")))

sink(file.path(RES, paste0(OUT_PREFIX_B, "_summary.txt")))
cat("=== Dataset B POSITIVE CONTROL: unmatched-MOFA vs matched ground truth ===\n\n")
cat(sprintf("Data: RNA 2000x12, Proteomics 2000x12, 3 conditions x4 (2 measured in both views). %d sampled pairings, K=3.\n", NP))
cat("Ground truth = matched-MOFA Factor1 loadings (dominant shared factor).\n\n")
cat("Recovery of the matched Factor1 signature (higher = better):\n"); print(tab, row.names = FALSE)
cat(sprintf("\nAVERAGE consensus vs condition-centroid: RNA r=%.3f prot r=%.3f", avg_vs_cen[1], avg_vs_cen[2]))
cat(sprintf("\nSelection score (balanced R2) range across %d pairings: %.1f%% .. %.1f%% (best pairing #%d)\n",
            NP, 100*min(score), 100*max(score), best))
cat("\nTop-30 RNA features recovered by AVERAGE that are in matched GT top-30:\n")
cat(paste(intersect(topK(avg_rna), gt_top_rna), collapse=", "), "\n")
sink()
cat(readLines(file.path(RES, paste0(OUT_PREFIX_B, "_summary.txt"))), sep = "\n")

## figure
png(file.path(FIGS, paste0(OUT_PREFIX_B, "_control.png")), width = 1500, height = 520, res = 150)
par(mfrow = c(1, 3), mar = c(4.2,4.2,3,1))
plot(gt_rna[names(avg_rna)], avg_rna, pch=16, col="#243b5399", cex=.4,
     xlab="matched GT Factor1 (RNA)", ylab="unmatched AVERAGE (RNA)",
     main=sprintf("RNA recovery (AVERAGE)\nr=%.2f, top30 Jaccard=%.2f", tab$rna_cor_vs_GT[1], tab$rna_top30_overlap[1]))
plot(gt_prot[names(avg_prot)], avg_prot, pch=16, col="#1FB8C499", cex=.4,
     xlab="matched GT Factor1 (prot)", ylab="unmatched AVERAGE (prot)",
     main=sprintf("Protein recovery (AVERAGE)\nr=%.2f, top30 Jaccard=%.2f", tab$prot_cor_vs_GT[1], tab$prot_top30_overlap[1]))
plot(sel_rna[names(gt_rna)], avg_rna[names(gt_rna)], pch=16, col="#E8735A99", cex=.4,
     xlab="SELECT best pairing (RNA)", ylab="AVERAGE (RNA)",
     main=sprintf("select vs average (RNA)\nr=%.2f", cor(sel_rna[names(gt_rna)], avg_rna[names(gt_rna)])))
dev.off()
cat("\nfigure ->", file.path(FIGS, paste0(OUT_PREFIX_B, "_control.png")), "\n")
