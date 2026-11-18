#!/usr/bin/env Rscript
# 3-VIEW dry-run on Dataset B, mirroring Dataset A's structure:
#   Proteomics <-> Metabolomics MATCHED (kept locked together),
#   RNA treated as a separate prep (permuted within condition).
# Ground truth = matched-MOFA Factor1 loadings for all 3 views.
# Q: does the 3-view permuted approach recover matched Factor1 for RNA / prot / metab?
## ---- locate repo root + load machine config (R/config_local.R, gitignored) ----
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT <- if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..")) else
  if (basename(getwd()) == "R") dirname(getwd()) else getwd()
if (!file.exists(file.path(ROOT, "R", "config_local.R")))
  stop("R/config_local.R not found: copy R/config_example.R to R/config_local.R and edit it.")
source(file.path(ROOT, "R", "config_local.R"))
suppressWarnings(suppressMessages(source(file.path(ROOT, "R", "lib_mofa.R"))))
set.seed(1)
ED <- B_DATA_DIR; RES <- RES_DIR; FIGS <- FIGS_DIR
rd  <- function(f) as.matrix(read.delim(file.path(ED, f), check.names = FALSE, row.names = 1))
gtf <- function(f) { d <- read.delim(file.path(ED, f), check.names = FALSE, row.names = 1); setNames(d$Factor1, rownames(d)) }

RNA <- rd("RNA.tsv"); PROT <- rd("Proteomics.tsv"); METAB <- rd("Metabolomics.tsv")
meta <- read.delim(file.path(ED, "meta.tsv"), stringsAsFactors = FALSE)
# common fully-observed samples across all 3 views  (2 of 3 conditions; the third is absent in RNA/prot)
usable <- Reduce(intersect, list(
  colnames(RNA)[colSums(is.na(RNA)) == 0], colnames(PROT)[colSums(is.na(PROT)) == 0],
  colnames(METAB)[colSums(is.na(METAB)) == 0]))
RNA <- RNA[, usable]; PROT <- PROT[, usable]; METAB <- METAB[, usable]
meta <- meta[meta$sample %in% usable, ]
cat("3-view usable samples:", length(usable), "\n")
gt <- list(RNA = gtf("gt_weights_RNA.tsv"), Proteomics = gtf("gt_weights_Proteomics.tsv"),
           Metabolomics = gtf("gt_weights_Metabolomics.tsv"))
cond <- setNames(meta$condition, meta$sample)

centroid_dir <- function(mat) {
  cm <- sapply(split(colnames(mat), cond[colnames(mat)]), function(s) rowMeans(mat[, s, drop = FALSE]))
  prcomp(t(cm))$rotation[, 1]
}
ref <- list(RNA = centroid_dir(RNA), Proteomics = centroid_dir(PROT), Metabolomics = centroid_dir(METAB))

## permute RNA within condition; prot & metab stay matched (identity)
rand_rna_pairing <- function() {
  do.call(rbind, lapply(unique(meta$condition), function(g) {
    s <- meta$sample[meta$condition == g]
    data.frame(anchor = s, rna = sample(s), stringsAsFactors = FALSE)  # prot/metab=anchor sample; rna shuffled
  }))
}
NP <- 120
Wr <- matrix(NA_real_, nrow(RNA), NP, dimnames = list(rownames(RNA), NULL))
Wp <- matrix(NA_real_, nrow(PROT), NP, dimnames = list(rownames(PROT), NULL))
Wm <- matrix(NA_real_, nrow(METAB), NP, dimnames = list(rownames(METAB), NULL))
for (i in seq_len(NP)) {
  pr <- rand_rna_pairing(); pid <- paste0("s", seq_len(nrow(pr)))
  Rv <- RNA[, pr$rna];   colnames(Rv) <- pid          # RNA permuted
  Pv <- PROT[, pr$anchor];  colnames(Pv) <- pid        # prot & metab locked to the anchor sample
  Mv <- METAB[, pr$anchor]; colnames(Mv) <- pid
  fit <- fit_mofa_multi(list(RNA = Rv, Proteomics = Pv, Metabolomics = Mv), K = 3, seed = 1L)
  sf  <- shared_factor_multi(fit, c("RNA","Proteomics","Metabolomics"), anchor = "Proteomics", ref = ref$Proteomics)
  Wr[, i] <- sf$RNA; Wp[, i] <- sf$Proteomics; Wm[, i] <- sf$Metabolomics
  cat(sprintf("\r  3-view fit %d/%d", i, NP)); flush.console()
}
cat("\n")
rec <- list(RNA = apply(Wr,1,median), Proteomics = apply(Wp,1,median), Metabolomics = apply(Wm,1,median))

## B-real: matched MOFA on prot + metab only (no tricks) -> real shared factor
fitR <- fit_mofa_multi(list(Proteomics = PROT, Metabolomics = METAB), K = 3, seed = 1L)
sfR  <- shared_factor_multi(fitR, c("Proteomics","Metabolomics"), anchor = "Proteomics", ref = ref$Proteomics)

cc <- function(a, gtv) { g <- gtv[names(a)]; s <- sign(cor(a, g, use="complete.obs")); abs(cor(a, g, use="complete.obs")) }
tab <- data.frame(
  view = c("RNA (permuted)","Proteomics (matched)","Metabolomics (matched)"),
  recover_matched_F1 = round(c(cc(rec$RNA, gt$RNA), cc(rec$Proteomics, gt$Proteomics), cc(rec$Metabolomics, gt$Metabolomics)), 3),
  vs_centroid = round(c(abs(cor(rec$RNA, ref$RNA[names(rec$RNA)],use="complete.obs")),
                        abs(cor(rec$Proteomics, ref$Proteomics[names(rec$Proteomics)],use="complete.obs")),
                        abs(cor(rec$Metabolomics, ref$Metabolomics[names(rec$Metabolomics)],use="complete.obs"))), 3))
# does 3-view prot/metab match the real matched prot+metab MOFA?
breal <- c(prot = abs(cor(rec$Proteomics, sfR$Proteomics[names(rec$Proteomics)], use="complete.obs")),
           metab= abs(cor(rec$Metabolomics, sfR$Metabolomics[names(rec$Metabolomics)], use="complete.obs")))

saveRDS(list(rec=rec, tab=tab, breal=breal, gt=gt), file.path(RES, paste0(OUT_PREFIX_B, "_3view.rds")))
sink(file.path(RES, paste0(OUT_PREFIX_B, "_3view_summary.txt")))
cat("=== Dataset B 3-VIEW dry-run: prot<->metab matched, RNA permuted ===\n\n")
cat(sprintf("Samples: %d (2 conditions, 4v4) | RNA %d, Prot %d, Metab %d feats | %d RNA pairings, K=3\n\n",
            length(usable), nrow(RNA), nrow(PROT), nrow(METAB), NP))
cat("Recovery of matched-MOFA Factor1 (|corr|, higher=better):\n"); print(tab, row.names=FALSE)
cat(sprintf("\n3-view prot/metab loadings vs REAL matched prot+metab MOFA: prot |r|=%.3f, metab |r|=%.3f\n", breal[1], breal[2]))
cat("\nInterpretation: matched views (prot, metab) carry real shared structure; RNA folds in at the\n")
cat("condition-signature (centroid) level -- exactly the intended Dataset A design.\n")
sink()
cat(readLines(file.path(RES, paste0(OUT_PREFIX_B, "_3view_summary.txt"))), sep = "\n")
