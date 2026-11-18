#!/usr/bin/env Rscript
# permuted-MOFA decisive experiment.
# Q: does the 36-pairing MOFA consensus contain cross-omics info BEYOND the group centroid?
# Outputs: results/consensus.rds, results/summary.txt, figs/decisive.png
## ---- locate repo root + load machine config (R/config_local.R, gitignored) ----
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT <- if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..")) else
  if (basename(getwd()) == "R") dirname(getwd()) else getwd()
if (!file.exists(file.path(ROOT, "R", "config_local.R")))
  stop("R/config_local.R not found: copy R/config_example.R to R/config_local.R and edit it.")
source(file.path(ROOT, "R", "config_local.R"))
suppressWarnings(suppressMessages({
  source(file.path(ROOT, "R", "lib_mofa.R"))
  source(file.path(ROOT, "R", "10_pairings.R"))
}))
set.seed(1)
RES <- RES_DIR; FIGS <- FIGS_DIR
dir.create(RES, showWarnings = FALSE); dir.create(FIGS, showWarnings = FALSE)

V <- load_views()
rmap <- read.delim(file.path(DATA, "rna_gene_map.tsv"),  stringsAsFactors = FALSE)
pmap <- read.delim(file.path(DATA, "prot_gene_map.tsv"), stringsAsFactors = FALSE)
rsym <- setNames(rmap$symbol, rmap$gene_id)
psym <- setNames(pmap$symbol, pmap$FeatureID)

## ---- enumerate pairings ----
pair_gr <- enumerate_pairings(V$rmeta, V$pmeta, respect_group = TRUE)     # 36
pair_gb_all <- enumerate_pairings(V$rmeta, V$pmeta, respect_group = FALSE)# 720
pair_gb <- pair_gb_all[sample(length(pair_gb_all), 36)]                   # null: 36 group-breaking
cat(sprintf("group-respecting pairings: %d | group-breaking sampled: %d\n", length(pair_gr), length(pair_gb)))

## ---- fit consensus over a set of pairings ----
run_consensus <- function(pairs, label) {
  W_rna <- matrix(NA_real_, nrow(V$rna), length(pairs), dimnames = list(rownames(V$rna), NULL))
  W_prot<- matrix(NA_real_, nrow(V$prot), length(pairs), dimnames = list(rownames(V$prot), NULL))
  for (i in seq_along(pairs)) {
    f <- fit_pairing(V$rna, V$prot, pairs[[i]], V$rmeta, K = 2, seed = 1L)
    W_rna[, i] <- f$rna; W_prot[, i] <- f$prot
    cat(sprintf("\r  [%s] fit %d/%d", label, i, length(pairs))); flush.console()
  }
  cat("\n")
  list(
    rna_med  = apply(W_rna, 1, median), rna_iqr  = apply(W_rna, 1, IQR),
    prot_med = apply(W_prot,1, median), prot_iqr = apply(W_prot,1, IQR),
    W_rna = W_rna, W_prot = W_prot)
}

t0 <- Sys.time()
con_gr <- run_consensus(pair_gr, "group-respecting")
con_gb <- run_consensus(pair_gb, "group-breaking null")
cat(sprintf("all fits done in %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))

## ---- B4 centroid baseline (group-mean difference per view) ----
cen <- centroid_signature(V$rna, V$prot, V$rmeta, V$pmeta)

## ---- comparisons ----
cc <- function(a, b) c(pearson = cor(a, b, use = "complete.obs"),
                       spearman = cor(a, b, method = "spearman", use = "complete.obs"))
cmp <- list(
  rna_consensus_vs_centroid  = cc(con_gr$rna_med,  cen$rna),
  prot_consensus_vs_centroid = cc(con_gr$prot_med, cen$prot),
  rna_gr_vs_gb   = cc(con_gr$rna_med,  con_gb$rna_med),
  prot_gr_vs_gb  = cc(con_gr$prot_med, con_gb$prot_med),
  rna_centroid_vs_gbnull  = cc(cen$rna,  con_gb$rna_med),
  prot_centroid_vs_gbnull = cc(cen$prot, con_gb$prot_med)
)

## ---- biology recovery: where do the known marker genes sit in the consensus? ----
rank_of <- function(vec, sym_map, gene, decreasing = TRUE) {
  s <- sym_map[names(vec)]
  ord <- order(vec, decreasing = decreasing)
  hit <- which(s[ord] == gene)
  if (!length(hit)) return(NA_integer_)
  hit[1]
}
bio <- data.frame(
  gene = c(A_MARKER_UP, A_MARKER_DOWN),
  rna_consensus_rank_up = c(rank_of(con_gr$rna_med, rsym, A_MARKER_UP),
                            rank_of(con_gr$rna_med, rsym, A_MARKER_DOWN, decreasing = FALSE)),
  prot_consensus_rank   = c(rank_of(con_gr$prot_med, psym, A_MARKER_UP),
                            rank_of(con_gr$prot_med, psym, A_MARKER_DOWN, decreasing = FALSE))
)

## ---- stability ----
stab <- c(rna_median_IQR  = median(con_gr$rna_iqr),
          prot_median_IQR = median(con_gr$prot_iqr),
          rna_med_absload = median(abs(con_gr$rna_med)),
          prot_med_absload= median(abs(con_gr$prot_med)))

saveRDS(list(con_gr = con_gr, con_gb = con_gb, cen = cen, cmp = cmp, bio = bio, stab = stab),
        file.path(RES, "consensus.rds"))

## ---- verdict ----
r1 <- cmp$rna_consensus_vs_centroid["pearson"]; r2 <- cmp$prot_consensus_vs_centroid["pearson"]
collapses <- (r1 > 0.95 && r2 > 0.95)
gr_gb_diff <- cmp$rna_gr_vs_gb["pearson"]

sink(file.path(RES, "summary.txt"))
cat("=== permuted-MOFA decisive experiment ===\n\n")
cat(sprintf("Views: RNA %d feats x 6 | protein %d feats x 6 (balanced 3v3)\n", nrow(V$rna), nrow(V$prot)))
cat(sprintf("Pairings: 36 group-respecting, 36 group-breaking (null), K=2 factors/fit\n\n"))
cat("B4 -- consensus vs group-centroid (THE decisive test):\n")
cat(sprintf("   RNA  loadings  Pearson=%.3f Spearman=%.3f\n", cmp$rna_consensus_vs_centroid[1], cmp$rna_consensus_vs_centroid[2]))
cat(sprintf("   prot loadings  Pearson=%.3f Spearman=%.3f\n\n", cmp$prot_consensus_vs_centroid[1], cmp$prot_consensus_vs_centroid[2]))
cat("N1 -- group-respecting vs group-breaking consensus (should differ if real):\n")
cat(sprintf("   RNA  Pearson=%.3f | prot Pearson=%.3f\n\n", cmp$rna_gr_vs_gb[1], cmp$prot_gr_vs_gb[1]))
cat("Stability (lower IQR = more stable across the 36 fits):\n")
cat(sprintf("   RNA median IQR=%.3f (median |load|=%.3f) | prot median IQR=%.3f (median |load|=%.3f)\n\n",
            stab[1], stab[3], stab[2], stab[4]))
cat("Biology recovery (rank in consensus, 1=strongest in expected direction):\n"); print(bio, row.names = FALSE)
cat("\n---------------------------------------------\n")
cat(sprintf("VERDICT: consensus %s the group centroid (RNA r=%.2f, prot r=%.2f).\n",
            ifelse(collapses, "COLLAPSES TO", "DIVERGES FROM"), r1, r2))
cat(sprintf("  -> %s\n", ifelse(collapses,
   "NO-GO: permutation machinery adds no info beyond a group-centroid analysis.",
   "PROMISING: consensus carries structure beyond the centroid -- investigate further.")))
sink()
cat(readLines(file.path(RES, "summary.txt")), sep = "\n")

## ---- figure ----
png(file.path(FIGS, "decisive.png"), width = 1500, height = 620, res = 150)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1))
plot(cen$rna, con_gr$rna_med, pch = 16, col = "#243b5399", cex = .5,
     xlab = "group-centroid loading (RNA)", ylab = "36-pairing consensus (RNA)",
     main = sprintf("RNA: consensus vs centroid\nr = %.2f", cmp$rna_consensus_vs_centroid[1]))
abline(0, coef(lm(con_gr$rna_med ~ cen$rna))[2], col = "#f0a202", lwd = 2)
plot(cen$prot, con_gr$prot_med, pch = 16, col = "#1FB8C499", cex = .5,
     xlab = "group-centroid loading (protein)", ylab = "36-pairing consensus (protein)",
     main = sprintf("protein: consensus vs centroid\nr = %.2f", cmp$prot_consensus_vs_centroid[1]))
abline(0, coef(lm(con_gr$prot_med ~ cen$prot))[2], col = "#f0a202", lwd = 2)
plot(con_gb$rna_med, con_gr$rna_med, pch = 16, col = "#E8735A99", cex = .5,
     xlab = "group-BREAKING consensus (RNA)", ylab = "group-respecting consensus (RNA)",
     main = sprintf("null control (RNA)\nr = %.2f", cmp$rna_gr_vs_gb[1]))
dev.off()
cat("\nfigure -> figs/decisive.png\n")
