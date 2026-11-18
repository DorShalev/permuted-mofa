#!/usr/bin/env Rscript
# Core combinatorics for permuted-MOFA: enumerate group-respecting (and
# group-breaking, for the null control) RNA<->protein sample pairings.
#
# For the experiment we down-select proteomics to a balanced 3v3 so both layers
# have equal n and every pairing is a clean bijection. With 3 + 3 samples per
# layer: 3! * 3! = 36 group-respecting pairings (fully enumerable).

suppressWarnings(library(utils))

#' All bijections (permutations) mapping vector `a` onto `b`.
#' @return list of data.frames with columns rna, prot
bijections <- function(a, b) {
  stopifnot(length(a) == length(b))
  perms <- gtools_permutations(length(b), length(b))  # rows = orderings of b
  lapply(seq_len(nrow(perms)), function(i)
    data.frame(rna = a, prot = b[perms[i, ]], stringsAsFactors = FALSE))
}

# tiny self-contained permutation generator (avoids a gtools dependency)
gtools_permutations <- function(n, r) {
  if (r == 1) return(matrix(seq_len(n), ncol = 1))
  sub <- gtools_permutations(n, r - 1)
  do.call(rbind, lapply(seq_len(n), function(i) {
    rest <- sub[apply(sub, 1, function(row) !(i %in% row)), , drop = FALSE]
    if (nrow(rest) == 0) return(NULL)
    cbind(i, rest)
  }))
}

#' Enumerate group-respecting pairings: group-g RNA -> group-g prot, for each group g.
enumerate_pairings <- function(rna_meta, prot_meta, respect_group = TRUE) {
  grps <- if (respect_group) unique(rna_meta$group) else "ALL"
  per_group <- lapply(grps, function(g) {
    if (respect_group) {
      a <- rna_meta$sample[rna_meta$group == g]
      b <- prot_meta$sample[prot_meta$group == g]
    } else {
      a <- rna_meta$sample; b <- prot_meta$sample
    }
    bijections(a, b)
  })
  # cartesian product of within-group bijections
  grid <- expand.grid(lapply(per_group, seq_along))
  lapply(seq_len(nrow(grid)), function(i) {
    do.call(rbind, Map(function(gl, idx) gl[[idx]], per_group, as.integer(grid[i, ])))
  })
}

if (sys.nframe() == 0) {
  ## ---- locate repo root + load machine config (R/config_local.R, gitignored) ----
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  ROOT <- if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..")) else
    if (basename(getwd()) == "R") dirname(getwd()) else getwd()
  if (!file.exists(file.path(ROOT, "R", "config_local.R")))
    stop("R/config_local.R not found: copy R/config_example.R to R/config_local.R and edit it.")
  source(file.path(ROOT, "R", "config_local.R"))
  rna  <- read.delim(file.path(DATA_DIR, "rna_meta_3v3.tsv"), stringsAsFactors = FALSE)
  prot <- read.delim(file.path(DATA_DIR, "prot_meta.tsv"), stringsAsFactors = FALSE)

  # balanced 3v3 protein subset (first 3 per group) so pairings are bijections
  set.seed(1)
  prot3 <- do.call(rbind, lapply(split(prot, prot$group),
                                 function(d) d[seq_len(3), ]))

  gr <- enumerate_pairings(rna, prot3, respect_group = TRUE)
  gb <- enumerate_pairings(rna, prot3, respect_group = FALSE)
  cat("group-respecting pairings:", length(gr), "(expect 36)\n")
  cat("group-breaking pairings  :", length(gb), "(expect 720 = 6!)\n")
  cat("\nExample group-respecting pairing #1:\n")
  print(gr[[1]], row.names = FALSE)
  # sanity: every RNA sample used exactly once in each pairing
  ok <- all(vapply(gr, function(p) !anyDuplicated(p$rna) && !anyDuplicated(p$prot), logical(1)))
  cat("\nall pairings are valid bijections:", ok, "\n")
  # sanity: group respected in every group-respecting pairing
  rg <- setNames(rna$group, rna$sample); pg <- setNames(prot3$group, prot3$sample)
  resp <- all(vapply(gr, function(p) all(rg[p$rna] == pg[p$prot]), logical(1)))
  cat("group respected in all 36:", resp, "\n")
}
