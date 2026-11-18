#!/usr/bin/env Rscript
# Shared helpers for permuted-MOFA. No side effects on source().
# Expects R/config_local.R to be loaded first (entry scripts do this).
if (!exists("DATA_DIR") || !exists("PYTHON_BIN"))
  stop("load R/config_local.R before sourcing lib_mofa.R")
suppressWarnings(suppressMessages({
  library(reticulate)
  Sys.setenv(RETICULATE_PYTHON = PYTHON_BIN)
  reticulate::use_python(PYTHON_BIN, required = FALSE)
  library(MOFA2)
}))

DATA <- DATA_DIR

#' Load views, select features. RNA = top-var logCPM; protein = well-covered top-var.
load_views <- function(n_rna = 2000, n_prot = 2000, min_prot_obs = 4) {
  r <- as.matrix(read.delim(file.path(DATA, "rna_logCPM_3v3.tsv"), check.names = FALSE, row.names = 1))
  p <- as.matrix(read.delim(file.path(DATA, "prot_intensity.tsv"),  check.names = FALSE, row.names = 1))
  rmeta <- read.delim(file.path(DATA, "rna_meta_3v3.tsv"), stringsAsFactors = FALSE)
  pmeta <- read.delim(file.path(DATA, "prot_meta.tsv"),   stringsAsFactors = FALSE)
  # balanced 3v3 protein subset (first 3 per group) -> clean bijections
  psub <- unlist(lapply(split(pmeta$sample, pmeta$group), function(s) s[1:3]))
  p <- p[, psub, drop = FALSE]; pmeta <- pmeta[match(psub, pmeta$sample), ]
  # feature selection
  rv <- apply(r, 1, var)
  r <- r[order(-rv)[seq_len(min(n_rna, sum(rv > 0)))], , drop = FALSE]
  keep_p <- rowSums(!is.na(p)) >= min_prot_obs & apply(p, 1, function(x) var(x, na.rm = TRUE) > 0)
  p <- p[keep_p, , drop = FALSE]
  pv <- apply(p, 1, var, na.rm = TRUE)
  p <- p[order(-pv)[seq_len(min(n_prot, nrow(p)))], , drop = FALSE]
  list(rna = r, prot = p, rmeta = rmeta, pmeta = pmeta)
}

#' Fit MOFA on one imputed pairing; return the group-factor loadings per view.
#' pairing: data.frame(rna=<rna sample>, prot=<prot sample>) in pseudo-sample order.
fit_pairing <- function(rna, prot, pairing, rmeta, K = 2, seed = 1L) {
  pid <- paste0("s", seq_len(nrow(pairing)))
  Rv <- rna[, pairing$rna, drop = FALSE];  colnames(Rv) <- pid
  Pv <- prot[, pairing$prot, drop = FALSE]; colnames(Pv) <- pid
  grp <- rmeta$group[match(pairing$rna, rmeta$sample)]  # pseudo-sample group (well-defined for group-respecting)

  obj <- create_mofa(list(RNA = Rv, protein = Pv))
  do  <- get_default_data_options(obj);  do$scale_views <- TRUE
  mo  <- get_default_model_options(obj); mo$num_factors <- K; mo$likelihoods[] <- "gaussian"
  to  <- get_default_training_options(obj)
  to$seed <- seed; to$convergence_mode <- "fast"; to$verbose <- FALSE; to$maxiter <- 1000
  obj <- suppressWarnings(prepare_mofa(obj, data_options = do, model_options = mo, training_options = to))
  tf  <- tempfile(fileext = ".hdf5")
  invisible(reticulate::py_capture_output(
    capture.output(mod <- run_mofa(obj, outfile = tf, use_basilisk = FALSE, save_data = FALSE))))
  on.exit(unlink(tf), add = TRUE)

  Z <- get_factors(mod)[[1]]                 # samples x K
  W <- get_weights(mod)                      # list: view -> features x K
  # identify the group factor: largest |mean(group1) - mean(group2)| in Z, scaled by sd
  is_g1 <- grp == A_GROUP1
  score <- apply(Z, 2, function(z) {
    m <- mean(z[is_g1]) - mean(z[!is_g1]); s <- sd(z); if (is.na(s) || s == 0) 0 else abs(m) / s
  })
  gf <- which.max(score)
  # sign so that positive loading = up in group 1
  sgn <- sign(mean(Z[is_g1, gf]) - mean(Z[!is_g1, gf])); if (sgn == 0) sgn <- 1
  list(
    rna  = setNames(W$RNA[, gf]     * sgn, rownames(rna)),
    prot = setNames(W$protein[, gf] * sgn, rownames(prot)),
    gf_score = score[gf], gf_index = gf
  )
}

#' Generic MOFA fit on an imputed pairing of two views (features x pseudo-samples,
#' same colnames). Returns loadings, per-view/per-factor variance, and factors.
fit_mofa_generic <- function(V1, V2, v1name = "V1", v2name = "V2", K = 3, seed = 1L) {
  vl <- list(V1, V2); names(vl) <- c(v1name, v2name)
  obj <- create_mofa(vl)
  do  <- get_default_data_options(obj);  do$scale_views <- TRUE
  mo  <- get_default_model_options(obj); mo$num_factors <- K; mo$likelihoods[] <- "gaussian"
  to  <- get_default_training_options(obj)
  to$seed <- seed; to$convergence_mode <- "fast"; to$verbose <- FALSE; to$maxiter <- 1000
  obj <- suppressWarnings(prepare_mofa(obj, data_options = do, model_options = mo, training_options = to))
  tf  <- tempfile(fileext = ".hdf5")
  invisible(reticulate::py_capture_output(
    capture.output(mod <- run_mofa(obj, outfile = tf, use_basilisk = FALSE, save_data = FALSE))))
  on.exit(unlink(tf), add = TRUE)
  r2 <- get_variance_explained(mod)$r2_per_factor[[1]]   # factors x views
  list(W = get_weights(mod), Z = get_factors(mod)[[1]], R2 = r2)
}

#' Pick the most-shared factor (max of the smaller per-view R2) and sign-align its
#' loadings so RNA loadings correlate positively with `ref` (a GT-independent direction).
shared_factor <- function(fit, v1 = "RNA", v2 = "protein", ref) {
  bal <- apply(fit$R2[, c(v1, v2), drop = FALSE], 1, min)   # balanced sharing per factor
  gf  <- which.max(bal)
  w1  <- fit$W[[v1]][, gf]; w2 <- fit$W[[v2]][, gf]
  common <- intersect(names(ref), names(w1))
  sgn <- sign(cor(w1[common], ref[common])); if (is.na(sgn) || sgn == 0) sgn <- 1
  list(v1 = setNames(w1 * sgn, rownames(fit$W[[v1]])),
       v2 = setNames(w2 * sgn, rownames(fit$W[[v2]])),
       balanced_r2 = bal[gf], gf = gf)
}

#' Fit MOFA on a list of >=2 views (each features x pseudo-samples, same colnames).
fit_mofa_multi <- function(views, K = 3, seed = 1L) {
  obj <- create_mofa(views)
  do  <- get_default_data_options(obj);  do$scale_views <- TRUE
  mo  <- get_default_model_options(obj); mo$num_factors <- K; mo$likelihoods[] <- "gaussian"
  to  <- get_default_training_options(obj)
  to$seed <- seed; to$convergence_mode <- "fast"; to$verbose <- FALSE; to$maxiter <- 1000
  obj <- suppressWarnings(prepare_mofa(obj, data_options = do, model_options = mo, training_options = to))
  tf  <- tempfile(fileext = ".hdf5")
  invisible(reticulate::py_capture_output(
    capture.output(mod <- run_mofa(obj, outfile = tf, use_basilisk = FALSE, save_data = FALSE))))
  on.exit(unlink(tf), add = TRUE)
  list(W = get_weights(mod), Z = get_factors(mod)[[1]],
       R2 = get_variance_explained(mod)$r2_per_factor[[1]])
}

#' Pick the most-shared factor across `vnames` (max of the smallest per-view R2),
#' sign-align loadings so view `anchor` correlates positively with `ref`.
shared_factor_multi <- function(fit, vnames, anchor, ref) {
  bal <- apply(fit$R2[, vnames, drop = FALSE], 1, min)
  gf  <- which.max(bal)
  wa  <- fit$W[[anchor]][, gf]
  common <- intersect(names(ref), rownames(fit$W[[anchor]]))
  sgn <- sign(cor(setNames(wa, rownames(fit$W[[anchor]]))[common], ref[common]))
  if (is.na(sgn) || sgn == 0) sgn <- 1
  loads <- lapply(vnames, function(v) setNames(fit$W[[v]][, gf] * sgn, rownames(fit$W[[v]])))
  names(loads) <- vnames
  c(loads, list(balanced_r2 = bal[gf], gf = gf))
}

#' Group-mean-difference signature per view (the B4 centroid baseline reduces to this).
centroid_signature <- function(rna, prot, rmeta, pmeta) {
  md <- function(mat, meta) {
    g <- meta$group[match(colnames(mat), meta$sample)]
    rowMeans(mat[, g == A_GROUP1, drop = FALSE], na.rm = TRUE) -
    rowMeans(mat[, g == A_GROUP2, drop = FALSE], na.rm = TRUE)
  }
  list(rna = md(rna, rmeta), prot = md(prot, pmeta))
}
