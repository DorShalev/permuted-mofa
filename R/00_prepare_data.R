#!/usr/bin/env Rscript
# Extract a self-contained snapshot of the Dataset A RNA (3v3) + proteomics
# matrices for the permuted-MOFA experiment. Base R only. Run once.
#
# Writes into DATA_DIR. The project reads only from there afterwards.
## ---- locate repo root + load machine config (R/config_local.R, gitignored) ----
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT <- if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..")) else
  if (basename(getwd()) == "R") dirname(getwd()) else getwd()
if (!file.exists(file.path(ROOT, "R", "config_local.R")))
  stop("R/config_local.R not found: copy R/config_example.R to R/config_local.R and edit it.")
source(file.path(ROOT, "R", "config_local.R"))

RNAD <- A_SRC_RNA_DIR
PROTD<- A_SRC_PROT_DIR
OUT  <- DATA_DIR
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- RNA (3v3) ----
fr <- read.delim(file.path(RNAD, "Datasets/final_results.tsv"), check.names = FALSE)
rna_samples <- A_RNA_SAMPLES
counts <- fr[, c("Gene", rna_samples)]
write.table(counts, file.path(OUT, "rna_counts_3v3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# normalized logCPM matrix (feature x sample)
norm <- read.delim(file.path(RNAD, "Datasets/rna_norm_TMMlogCPM.tsv"), check.names = FALSE)
write.table(norm, file.path(OUT, "rna_logCPM_3v3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# annotation (gene_id -> symbol) so views can be labeled
ann <- read.csv(file.path(RNAD, "Enrichment/gene_annotation.csv"), stringsAsFactors = FALSE)
write.table(ann[, c("gene_id","symbol")], file.path(OUT, "rna_gene_map.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

rna_meta <- data.frame(
  sample = rna_samples,
  group  = ifelse(grepl(A_RNA_G2_REGEX, rna_samples), A_GROUP2, A_GROUP1),
  layer  = "RNA", stringsAsFactors = FALSE)
write.table(rna_meta, file.path(OUT, "rna_meta_3v3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- Proteomics (5 vs 5) ----
fp <- read.delim(file.path(PROTD, "Datasets/final_results.tsv"), check.names = FALSE)
prot_samples <- A_PROT_SAMPLES
prot_samples <- intersect(prot_samples, colnames(fp))
prot <- fp[, c("FeatureID", prot_samples)]
write.table(prot, file.path(OUT, "prot_intensity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# FeatureID -> gene symbol (split "SYMBOL | description")
gmap <- data.frame(FeatureID = fp$FeatureID,
                   symbol = trimws(sub("\\|.*", "", fp$Genes)),
                   stringsAsFactors = FALSE)
write.table(gmap, file.path(OUT, "prot_gene_map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

prot_meta <- data.frame(
  sample = prot_samples,
  group  = ifelse(grepl(A_PROT_G2_REGEX, prot_samples), A_GROUP2, A_GROUP1),
  layer  = "protein", stringsAsFactors = FALSE)
write.table(prot_meta, file.path(OUT, "prot_meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Data snapshot written to", OUT, "\n")
cat(sprintf("  RNA logCPM : %d features x %d samples\n", nrow(norm), length(rna_samples)))
cat(sprintf("  Protein    : %d features x %d samples\n", nrow(prot), length(prot_samples)))
cat(sprintf("  shared symbols (RNA∩prot): %d\n",
    length(intersect(na.omit(ann$symbol), na.omit(gmap$symbol)))))
