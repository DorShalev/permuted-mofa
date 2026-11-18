#!/usr/bin/env Rscript
# 10-sample MOFA factor plot for Dataset A: matched Proteomics+Metabolomics core (all 10),
# balanced 500-feature/view + scaled — same setup as the RNA-folded run (RNA contributes
# via loadings; the 10 sample scores come from the matched prot+metab core).
## ---- locate repo root + load machine config (R/config_local.R, gitignored) ----
.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT <- if (length(.f)) normalizePath(file.path(dirname(.f[1]), "..")) else
  if (basename(getwd()) == "R") dirname(getwd()) else getwd()
if (!file.exists(file.path(ROOT, "R", "config_local.R")))
  stop("R/config_local.R not found: copy R/config_example.R to R/config_local.R and edit it.")
source(file.path(ROOT, "R", "config_local.R"))
suppressWarnings(suppressMessages(source(file.path(ROOT, "R", "lib_mofa.R"))))
suppressWarnings(suppressMessages(library(ggplot2)))
set.seed(1)
LD <- A_MO_DATA_DIR; FIGS <- FIGS_DIR
rd <- function(f) as.matrix(read.delim(file.path(LD, f), check.names = FALSE, row.names = 1))
hv <- function(m, n) { v <- apply(m, 1, var, na.rm = TRUE); m[order(-v)[seq_len(min(n, sum(v > 0, na.rm = TRUE)))], , drop = FALSE] }

PROT <- rd("prot_expr.tsv"); METAB <- rd("metab_expr.tsv")
NFEAT <- 500
keepP <- rowSums(!is.na(PROT)) >= 7 & apply(PROT, 1, function(x) var(x, na.rm = TRUE) > 0)
PROT  <- hv(PROT[keepP, , drop = FALSE], NFEAT)
METAB <- METAB[apply(METAB, 1, function(x) var(x, na.rm = TRUE) > 0), , drop = FALSE]

G1 <- A_MO_GROUP1; G2 <- A_MO_GROUP2
grp <- function(s) ifelse(grepl(A_MO_G2_REGEX, s), G2, G1)
ref_prot <- { cm <- sapply(split(colnames(PROT), grp(colnames(PROT))),
                           function(s) rowMeans(PROT[, s, drop = FALSE], na.rm = TRUE)); cm[,G1] - cm[,G2] }
fit <- fit_mofa_multi(list(Proteomics = PROT, Metabolomics = METAB), K = 4, seed = 1L)
sf  <- shared_factor_multi(fit, c("Proteomics","Metabolomics"), anchor = "Proteomics", ref = ref_prot)
gf  <- sf$gf                                            # group-separating factor index
Z   <- fit$Z                                            # samples x factors
tot <- colSums(fit$R2)                                  # total variance per factor
f2  <- setdiff(order(-tot), gf)[1]                      # 2nd factor = next most variance
cat(sprintf("group factor = F%d (bal R2 %.2f); 2nd = F%d\n", gf, sf$balanced_r2, f2))
cat("per-factor total R2:", paste(round(tot,1), collapse=" "), "\n")

samp <- rownames(Z)
df <- data.frame(x = Z[, gf], y = Z[, f2], Group = grp(samp), label = samp)
# orient so group 2 is on the positive side (fixed plotting convention)
if (mean(df$x[df$Group==G2]) < mean(df$x[df$Group==G1])) df$x <- -df$x

p <- ggplot(df, aes(x, y, color = Group)) +
  geom_vline(xintercept = 0, color = "grey85") + geom_hline(yintercept = 0, color = "grey85") +
  geom_point(size = 5) +
  ggrepel::geom_text_repel(aes(label = label), size = 4.3, seed = 1, show.legend = FALSE) +
  scale_color_manual(values = setNames(c("#E8735A", "#1FB8C4"), c(G2, G1))) +
  labs(title = "MOFA2: Factor 1 vs Factor 2", subtitle = "all 10 matched protein + metabolite samples (RNA folded in)",
       x = "Factor 1", y = "Factor 2") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 12, color = "grey40"))

ggsave(file.path(FIGS, "mofa_factors_10.png"), p, width = 8, height = 6, dpi = 150, bg = "white")
cat("saved", file.path(FIGS, "mofa_factors_10.png"), "\n")
