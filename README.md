# permuted-MOFA — consensus MOFA for *unmatched* multi-omics (small-n)

> **Status: early-stage research code.** Tested only on two small unpublished datasets — **Dataset A** (RNA + proteomics, 3 vs 3) and **Dataset B** (RNA + proteomics, 4 vs 4) — i.e. n=3–4 per group. The decisive test on Dataset A (n=3) was **NO-GO**: the consensus collapsed to the group centroid, and the Dataset B positive control (n=4) recovered only the condition-centroid signal (see `FINDINGS.md`). This is not a validated integration method; verification and further testing on larger cohorts are planned. The input data, result objects and figures come from unpublished collaborator datasets and are deliberately **not** included in this repository.

**Status:** ✅ DECISIVE TEST DONE → **NO-GO at n=3** (consensus ≡ group centroid, r=1.00/0.98). See `FINDINGS.md`.
Novelty scan done (`notes/literature_scan.md`): idea is novel but predicted to collapse — and it did.
**Owner:** Dor Shalev · started 2026-07-05
**Isolated on purpose:** standalone project with its own R environment — no ties to any production pipeline. This is a method experiment, not production code.

---

## The idea (being tested)

We have two omics layers measured on **different, unmatched biological samples** that share the same design:
RNA-seq (3 vs 3) and proteomics (5 vs 5, separate preps) from **Dataset A** (unpublished; condition 1 vs condition 2). Standard joint-factor
models (MOFA, DIABLO, MEFISTO) need **matched samples** — same biological unit in every layer — which we don't have.

**Proposal:** impute a sample correspondence, but don't trust any single one:

1. Enumerate every **group-respecting** RNA↔protein pairing (condition 1→condition 1, condition 2→condition 2). For 3v3 that's `3! × 3! = 36` pairings — fully enumerable, no sampling.
2. For each pairing, treat the paired columns as if matched and fit **MOFA2** (seeded).
3. Align factors across the 36 fits (sign/order ambiguity) and **aggregate loadings (median)** → a *consensus* cross-omics factor, plus a per-feature **stability** score (IQR across pairings).

Intuition: MOFA learns from **feature covariance**, not sample identity, and it does **not** require
matched *features* — so it can use the full RNA (~22k genes) and protein (~6k) spaces, not just the ~5.9k shared symbols.
Averaging over pairings marginalises out the arbitrary within-group correspondence.

## The honest open question

At n=3/group the only *reliable* shared structure is the **group axis** (condition 1 vs condition 2). Within-group covariance from 3
points may overfit. So the experiment must decide: does the consensus factor recover **real cross-omics signal beyond
the group contrast**, or does it just re-derive what simple feature/pathway-level integration already gives?

The design below is built to give a clean yes/no — with the baselines and null controls a methods paper would demand.

---

## Experimental design

### Method under test
`permMOFA`: 36-pairing enumeration → per-pairing MOFA2 → factor alignment → median-aggregated loadings + stability.

### Baselines it must beat (else it's not a new method)
- **B4 — group-centroid MOFA `★ THE DECISIVE TEST` :** collapse each layer to 1 condition-1-mean + 1 condition-2-mean column, run MOFA once on that 2-sample pseudo-pairing. The literature scan independently predicts the 36-pairing consensus will **equal this centroid** (under random group-respecting pairing, within-group cross-omics coupling averages to the group-mean coupling). If consensus ≈ centroid → the permutation machinery adds cost, not information → method fails its reason to exist. **Run this first — it's cheap and decisive.**
- **B1 — feature-level:** per-gene RNA logFC vs protein logFC correlation (already done: 5,874 genes, r≈0.15 global, 85% concordant among 52 co-significant).
- **B2 — group-mean PLS:** concatenate group-mean profiles, PLS-DA / correlate loadings across layers.
- **B3 — DE-rank overlap:** overlap of single-view DE rankings.

### Implementation hazard to handle (flagged by scan)
MOFA factors are identified only up to **sign / scale / order** and reorder between the 36 fits. Averaging loadings requires a **factor-matching + sign-alignment step first** (Procrustes or Hungarian, anchored to the group-contrast direction) — otherwise the "average" is meaningless. This is handled in `R/lib_mofa.R` (group-factor selection + sign anchoring) before aggregation.

### Null controls (credibility)
- **N1 — group-breaking pairings:** pair condition-1 RNA↔condition-2 protein. Consensus factor should collapse. If group-respecting ≈ group-breaking → no real cross-omics content.
- **N2 — label shuffle:** permute group labels, rerun; biology should vanish.
- **N3 — within-view feature permutation** null.

### Evaluation
- **Stability:** are consensus loadings stable across the 36 pairings (low IQR), or pairing-dependent noise?
- **Biological recovery:** do top consensus loadings recover the marker genes and pathway module already identified by per-layer DE (names withheld; unpublished data)? Do RNA-view and protein-view loadings of the shared factor agree in sign for known-concordant genes?
- **Value-add:** does it surface concordant features that B1 misses (esp. single-layer features, or pairs B1 ranks low that validate)?

### Go / no-go gate
Stable **AND** group-respecting ≫ group-breaking **AND** recovers biology **AND** adds value over B1 →
promising → develop (simulations + 2nd dataset + write-up). Otherwise → documented negative result; feature/pathway-level
integration remains the route for the dataset itself. Both outcomes are useful evidence; only the first is a new *method*.

---

## Data (not in this repository)
`data/`, `results/` and `figs/` are gitignored. `R/00_prepare_data.R` extracts a self-contained snapshot of
Dataset A from upstream processing outputs into `data/` (the project does not read the upstream outputs at runtime):
- `rna_logCPM_3v3.tsv`, `rna_counts_3v3.tsv`, `rna_meta_3v3.tsv`, `rna_gene_map.tsv` — RNA 3 vs 3
- `prot_intensity.tsv`, `prot_meta.tsv`, `prot_gene_map.tsv` — proteomics 5 vs 5

> Both datasets are unpublished. No raw values, sample names or results are committed.

## Configuration
All machine paths, sample names and group labels live in `R/config_local.R` (gitignored).
Copy `R/config_example.R` to `R/config_local.R` and fill it in; every script sources it via a path
relative to the repo root, so run them from anywhere as `Rscript R/<script>.R`.

## Scripts (`R/`)
- `config_example.R` — template for `config_local.R`.
- `lib_mofa.R` — shared MOFA2 helpers (fit, group-factor selection, sign alignment, centroid signature).
- `00_prepare_data.R` — Dataset A: extract & save the matrices above.
- `10_pairings.R` — enumerate group-respecting / group-breaking pairings (run alone = sanity check).
- `20_run_experiment.R` — Dataset A decisive test (B4 centroid, N1 null, stability, marker recovery).
- `40_datasetB_control.R` — Dataset B positive control: AVERAGE vs SELECT vs matched ground truth.
- `50_datasetB_3view.R` — Dataset B 3-view dry run (prot↔metab matched, RNA permuted).
- `60_datasetA_3view.R`, `70_datasetA_mo_figures.R`, `80_datasetA_factorplot_10.R` — Dataset A 3-omics run + figures.

## Dependencies
Needs `MOFA2` (Bioconductor) + `mofapy2` (python via reticulate; set `PYTHON_BIN` in the config),
plus `ggplot2` and `ggrepel` for script 80. Use a project-specific R library.

## Literature / novelty
See `notes/literature_scan.md` (step 0 — running). Do not invest in development until the novelty verdict is in.
