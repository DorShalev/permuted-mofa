# Prior-Art Scan: Permutation-Averaged / Consensus MOFA for Unpaired Small-n Multi-Omics
_Step 0 novelty check · 2026-07-05 · background research agent (12 web sources)_

## VERDICT

The **specific recipe** — enumerate all group-respecting sample-to-sample pairings between two unpaired
omics, fit MOFA on each imputed pairing, then average the factor loadings into a consensus factor — does
**not appear to exist as a named, published method.** It is a genuinely novel *combination*. But each
ingredient is well-trodden: (a) the underlying "diagonal"/unpaired-integration problem is a large active
field; (b) MOFA and relatives (DIABLO, JIVE, iCluster, MEFISTO) **require row-matched samples**, so bolting
permutation-imputation onto MOFA is a workaround for a limitation the community usually solves with optimal
transport or adversarial autoencoders; and (c) "solve over unknown correspondence" exists in the
Wasserstein–Procrustes literature, but as an *optimization* over permutations, not a *uniform average* over
all valid ones. Biggest concern is not novelty but **soundness at n=3**: unpaired integration is treated as
data-hungry (built for thousands of single cells); averaging over an arbitrary correspondence at n=3/group is
likely statistically ill-posed. Novel — but publishing into a headwind of "why not just correlate at the
group/feature level."

## Closest prior art
1. **MOFA / MOFA2** (Argelaguet 2018/2020) — the tool being extended. Tolerates missing values but has **no**
   mechanism for unknown row correspondence. A Biostars thread (#9619965) shows users asking exactly how to
   feed partially-unmatched RNA vs proteomics into MOFA2 with no clean answer → the gap this targets.
2. **MODIS** (arXiv:2503.18856, 2025) — closest *problem framing*: small + unpaired + **bulk** (TCGA). But
   different engine: multi-VAE + adversarial modality classifier + relativistic GAN, diagonal integration via
   shared class labels + large reference DB. Nearest competitor to differentiate against.
3. **Wasserstein–Procrustes / OT over unknown permutations** (arXiv:2103.07458) — principled version: jointly
   optimize transform + permutation P. Our novelty seam = **average over all valid P uniformly** (marginalize)
   instead of optimizing one. Invites reviewer Q: "why average instead of learn the good pairings?"
4. **OT diagonal-integration family (SCOT, Pamona, UnionCom, MMD-MA, GLUE)** — already solve "no correspondence"
   via soft correspondence, but all built/benchmarked for **single-cell scale (1000s)**, not n=3.

## Landscape (matched-sample requirement)
- **Require matched samples (bulk/small-n friendly):** MOFA/MOFA2/MOFA+, MEFISTO, DIABLO, JIVE/AJIVE,
  iCluster(+), SNF.
- **Unpaired but single-cell scale:** MMD-MA, UnionCom, SCOT, Pamona, GLUE/scGLUE, MultiMAP, bindSC,
  Seurat-bridge, scTopoGAN.
- **Unpaired + bulk + small-n:** essentially only MODIS (different method). → this is the empty cell we sit in.

## What's genuinely novel here
- **Exhaustive enumeration is tractable *because* n is tiny** (3×3 → 36 pairings). OT/manifold methods never
  enumerate. "Small enough to enumerate the whole correspondence group" is absent from the literature.
- **Averaging the MODEL over all valid permutations** (marginalization), vs selecting/optimizing one, was not
  found for MOFA or multi-omics FA. Framing = crude Bayesian marginalization over a uniform prior on
  group-respecting P. No hits for "consensus/ensemble/permutation MOFA" as named methods.

## Risks the literature flags (must pre-empt)
1. **Data-hungry premise** — at n=3 there's little information to identify correspondence; may be ill-posed.
2. **★ Averaging may converge to the group-centroid** — under random group-respecting pairing, within-group
   cross-omics coupling averages (in expectation) to the group-mean coupling. **The consensus factor may equal
   just pairing group-mean RNA with group-mean protein.** → MUST test consensus vs n=1-per-group centroid MOFA.
3. **Factor identifiability** — MOFA factors are up to sign/scale/rotation and reorder between runs. Averaging
   36 independent fits requires factor-matching + sign alignment (Procrustes/Hungarian) FIRST, or the average
   is meaningless.
4. **Small-n replicability** — omics below a few hundred samples: low sensitivity, poor replicability; unpaired
   worsens it.
5. **Diagonal integration can manufacture false correspondences** (Bioinformatics btaf345, 2025).

## Bottom line
Method is unclaimed → citable as novel. Pitch as **"small-n permutation-marginalized MOFA."** Pre-empt two
objections: (1) show it beats a **group-centroid / DIABLO** baseline (risk #2), and (2) handle **factor
sign/order alignment** before averaging (risk #3). Differentiate vs MODIS (small unpaired bulk, different
engine) and the OT family (unpaired but single-cell, and they *optimize* not *average* correspondence).

## Key citations
- MOFA: https://www.embopress.org/doi/10.15252/msb.20178124 · PMC6010767 · MOFA+ PMC7212577 · docs https://biofam.github.io/MOFA2/
- Biostars unmatched RNA/prot in MOFA2: https://www.biostars.org/p/9619965/
- MODIS: https://arxiv.org/abs/2503.18856
- Pamona: PMC8696097 · SCOT: RG 358006341 · GLUE: https://www.nature.com/articles/s41587-022-01284-4 (PMC9546775)
- scTopoGAN: https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad171/7450146
- Diagonal false-mapping caution: https://academic.oup.com/bioinformatics/article/41/6/btaf345/8162456
- Joint multi-omics DR benchmark: https://www.nature.com/articles/s41467-020-20430-7
- OT over unknown permutations: https://arxiv.org/pdf/2103.07458 · Procrustes–Wasserstein: https://arxiv.org/html/2510.01706
- Small-n omics replicability: https://arxiv.org/pdf/2504.20431 · MO integration review: https://arxiv.org/pdf/2501.17729
