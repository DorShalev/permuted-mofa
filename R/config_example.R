# Machine-specific configuration for permuted-MOFA.
#
# Copy this file to R/config_local.R (gitignored) and fill in your own values.
# Every script sources R/config_local.R via a path relative to the repo root,
# so nothing machine- or dataset-specific lives in the committed code.
#
# Two datasets are referred to throughout:
#   Dataset A — RNA + proteomics (3 vs 3 after balancing), plus a matched
#               metabolomics layer used by scripts 60-80. No sample matching
#               between RNA and the other layers.
#   Dataset B — RNA + proteomics + metabolomics, matched samples, 4 vs 4 in the
#               two conditions measured in every layer (positive control).

## ---- Python used by MOFA2 via reticulate ----
PYTHON_BIN <- "/usr/bin/python3"

## ---- project directories (default: inside the repo, all gitignored) ----
DATA_DIR <- file.path(ROOT, "data")      # Dataset A RNA/protein snapshot (written by 00)
RES_DIR  <- file.path(ROOT, "results")
FIGS_DIR <- file.path(ROOT, "figs")
A_MO_DATA_DIR <- file.path(ROOT, "data", "datasetA_mo")  # Dataset A 3-omics matrices (60-80)
B_DATA_DIR    <- file.path(ROOT, "data", "datasetB")     # Dataset B matrices + matched-MOFA GT weights

## ---- output file prefixes ----
OUT_PREFIX_A <- "datasetA"
OUT_PREFIX_B <- "datasetB"

## ---- Dataset A: source pipeline outputs read once by 00_prepare_data.R ----
A_SRC_RNA_DIR  <- "/path/to/datasetA/RNA/results"         # contains Datasets/, Enrichment/
A_SRC_PROT_DIR <- "/path/to/datasetA/Proteomics/results"  # contains Datasets/

## ---- Dataset A: sample names and group labels ----
# Loadings are signed so that positive = higher in A_GROUP1.
A_GROUP1 <- "COND1"
A_GROUP2 <- "COND2"
A_RNA_SAMPLES  <- c("rna_c2_1", "rna_c2_2", "rna_c2_3", "rna_c1_1", "rna_c1_2", "rna_c1_3")
A_RNA_G2_REGEX <- "_c2_"                 # matches the A_GROUP2 RNA samples
A_PROT_SAMPLES <- c("C2_1", "C2_2", "C2_3", "C2_4", "C2_5", "C1_1", "C1_2", "C1_3", "C1_4", "C1_5")
A_PROT_G2_REGEX <- "^C2"                 # matches the A_GROUP2 protein samples

# 3-omics matrices (60-80): group is derived from the sample name.
A_MO_GROUP1 <- "COND1"
A_MO_GROUP2 <- "COND2"
A_MO_G2_REGEX <- "C2"                    # matches A_MO_GROUP2 sample names

# Marker genes checked for "biological recovery" in 20_run_experiment.R
# (expected up / down in A_GROUP1, taken from the per-layer DE).
A_MARKER_UP   <- "GENE_UP"
A_MARKER_DOWN <- "GENE_DOWN"
