# Findings — permuted-MOFA decisive test (2026-07-05)

Datasets: **Dataset A** (RNA + proteomics, 3 vs 3, unpublished) and **Dataset B** (RNA + proteomics,
4 vs 4, unpublished). Condition labels and sample names are withheld.

## Result: NO-GO (at n=3). The consensus collapses to the group centroid.

We fit MOFA on all **36 group-respecting** RNA↔protein pairings (+ 36 group-breaking as a null),
extracted the group factor from each, sign-aligned, and median-aggregated the loadings. Then we compared
the consensus to the **group-centroid** signature (per-view condition-1 mean − condition-2 mean) on Dataset A.

| Test | RNA | protein | reading |
|---|---|---|---|
| **B4 — consensus vs group-centroid** | **r = 1.00** | **r = 0.98** | consensus ≡ centroid |
| N1 — group-respecting vs group-breaking (RNA) | r = 1.00 | — | RNA loadings driven by RNA's own group axis, pairing-invariant |
| N1 — group-respecting vs group-breaking (protein) | — | r = 0.14 | breaking groups destroys protein signal (null behaves) |
| Stability (median IQR across 36 fits) | 0.048 | 0.160 | stable — but stably equal to the centroid |
| Biology: known up-marker (from per-layer DE) rank in consensus | 72 (RNA) | 46 (prot) | NOT surfaced at top |
| Biology: known down-marker (from per-layer DE) rank in consensus | 25 (RNA) | — | not surfaced |

See `figs/decisive.png` — consensus points lie on the identity line with the centroid.

## Why it collapses (the math, confirmed empirically)

Within-group pairings are arbitrary. Averaging MOFA over all of them **marginalizes out** any
pairing-specific cross-view coupling; in expectation only the **group-mean difference** in each view
survives. So the 36-pairing consensus = each view's independent group centroid. The "integration" does
nothing a per-layer DE (group means) doesn't already give. The RNA null nails it: even *breaking* the
group correspondence leaves the RNA loadings unchanged (r=1.00) — the RNA view finds its own group axis
regardless of what protein it's paired with. There is no emergent cross-omics structure to find at n=3.

This confirms both the a-priori prediction and the literature-scan warning (risk #2).

## Honest verdict on the idea
The instinct — "MOFA learns from cross-feature covariance, not sample identity" — is correct in general,
but at **n=3/group there is no within-group covariance to learn**, so the only shared structure is the
group contrast, and averaging recovers exactly that. The method is not wrong; it's **degenerate at this
sample size**. It would need many samples per group (the regime where the single-cell OT methods already
operate) for within-group cross-omics coupling to be estimable — and there, selecting/optimizing the
correspondence (Wasserstein–Procrustes / OT) beats uniform averaging.

## What this is good for
- A clean **cautionary/negative result**: *"permutation-marginalized MOFA on unmatched small-n data
  provably reduces to the group centroid"* — publishable as a short methods note / caution, with the
  36-pairing proof + the null control. Not a new working integration method.
- For the dataset itself, feature-level + pathway-level integration (already done) DID surface the
  cross-omics anchors (marker genes and a pathway module) that this consensus buried.

---

# Follow-up 1 — SELECT vs AVERAGE (Wasserstein–Procrustes style)

Instead of averaging over pairings, **select** the single pairing whose top shared factor has the highest
balanced cross-view variance (the "optimize the correspondence" idea). On both Dataset A (3v3) and Dataset B (4v4):
selection lands in the same place as averaging — **no advantage, no divergence from the centroid**. At small
n there is no signal to select *toward*; the best-of-N pairing is best by chance, not by biology.

# Follow-up 2 — Dataset B POSITIVE CONTROL (the important one)

Took Dataset B's **matched** multi-omics (RNA + Proteomics, 4 condition-1 vs 4 condition-2 samples — a third
condition is unmeasured in these two views), where matched-MOFA **Factor1** is known ground truth (the dominant shared factor). Deliberately
**unmatched** it (permute within-condition sample correspondence), ran the same approach, and checked whether
we recover the matched Factor1 features.

| Approach | RNA \|cor\| vs matched F1 | prot \|cor\| vs matched F1 | GT top-30 found in recovered top-100 |
|---|---|---|---|
| **AVERAGE (unmatched)** | **0.999** | **0.98** | RNA 17/30 · prot 15/30 |
| SELECT (best pairing) | 0.99 | 0.99 | similar |
| condition-centroid | 0.996 | 0.987 | similar |
| _AVERAGE vs centroid_ | _0.998_ | _0.997_ | — |

**So: YES — treating Dataset B as unmatched recovers the matched Factor1 signature almost perfectly** (loading
direction r≈0.99; ~most important features captured within a wider top-N net, though the exact top-30 ranking
reshuffles). **But it works for the same reason it "failed" on Dataset A:** unmatched-average = condition-centroid
(r=0.998), and Dataset B's Factor1 *is* the condition axis. It recovers the treatment-associated features
because those equal the centroid — not by learning hidden cross-omics coupling.

## Unified conclusion
Unmatched-MOFA (average **or** select) ≈ **per-view group/condition centroid**. When the biologically
interesting factor is the treatment contrast (true for both Dataset A and Dataset B), the method
**does recover those features** — validated against Dataset B's matched ground truth. It will **not** find structure
orthogonal to the group contrast, and selection does not help at small n.

## Practical upshot (for the stated goal: "get interesting features from unmatched data")
Usable and validated for extracting the **condition-associated multi-omics feature list** — but it is
mathematically equivalent to *per-layer DE → combine* / a group-centroid analysis. Present it as that, not as
MOFA discovering latent cross-omics coupling. For small unmatched datasets the simpler feature/pathway
integration is the honest, equivalent, and more interpretable route.

## Reproduce
Set up `R/config_local.R` (see `R/config_example.R`), then
`Rscript R/20_run_experiment.R` (Dataset A) · `Rscript R/40_datasetB_control.R` (Dataset B control) →
`results/*.txt`, `results/*.rds`, `figs/decisive.png`, `figs/<OUT_PREFIX_B>_control.png`.
