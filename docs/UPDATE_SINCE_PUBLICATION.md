# Updates since publication

**Last updated:** 2026-08-08

## Read this first

This repository is cited in the Data Availability Statement of:

> Ruhle, M., et al. (2025). *Leakage-aware machine learning reveals structured clinical and vaginal
> microbiome patterns associated with preterm birth in a Mexican cohort.*
> Front. Glob. Women's Health 7:1799518. doi: 10.3389/fgwh.2026.1799518

**After the article was published, the analysis code was reviewed and several methodological updates
were made. Some of these updates change the reported numbers.** The published article is not being
amended; this repository is the authoritative record of the updated analysis, and this document is
the changelog between the two.

If you are reading the article and want the current numbers, use the tables in this document, or the
results table in [`README.md`](../README.md#results-current) — not the figures in the published PDF,
and not the archived report in
[`results/published_version/`](../results/published_version/), which predates these updates.

Nothing here changes the study design, the cohort, the data, or the biological question. The updates
are refinements to the analysis code.

---

## 1. Headline: what changed in the results

### 1.1 The best-performing model changed identity

| | Model | Clinical variable set | Microbiome input | AUROC |
|---|---|---|---|---|
| **Published** | Random Forest | Approach 3 (data-driven) | Full microbiome | 0.813 |
| **Updated** | **Elastic net (glmnet)** | **Approach 3 (data-driven)** | **ANCOM-BC2-selected taxa** | **0.760 ± 0.270** |

*On the published value:* two different figures for the previous best model used to circulate in
this repository — **0.813**, the value the permutation-test script referenced, and **0.849**, which
appeared in the `README.md` results table and in the archived report now kept under
`results/published_version/`. The 0.849 figure comes from an older analysis run made before a
contaminant-filtering step was applied (59 taxa instead of 49) and is superseded; 0.813 is the
post-filtering value. Both predate the updates below. `README.md` now reports the updated numbers,
and the earlier render is labelled as superseded.

The previously reported best combination (Random Forest / Approach 3 / full microbiome) now scores
**0.680** and ranks 6th of 12. This reordering follows directly from update 2.1 (the CLR
transformation): once the microbiome is compositionally transformed, the relative standing of the
model combinations changes (Spearman ρ ≈ 0.47 against the earlier ranking; largest single shift
ΔAUROC = 0.24).

### 1.2 Full current results (real cohort data, all 12 combinations)

Nested cross-validation, 5 outer folds, subject-level splits stratified by outcome.
AUROC is the mean across outer folds ± SD.

| Model | Clinical variable set | Microbiome input | AUROC |
|---|---|---|---|
| glmnet | Approach 3 (data-driven) | ANCOM taxa | **0.760 ± 0.270** |
| Random Forest | Approach 3 (data-driven) | ANCOM taxa | 0.736 ± 0.147 |
| glmnet | Approach 2 (literature) | ANCOM taxa | 0.724 ± 0.187 |
| Random Forest | Approach 1 (DREAM) | ANCOM taxa | 0.718 ± 0.172 |
| glmnet | Approach 2 (literature) | Full microbiome | 0.682 ± 0.105 |
| Random Forest | Approach 3 (data-driven) | Full microbiome | 0.680 ± 0.201 |
| Random Forest | Approach 2 (literature) | ANCOM taxa | 0.680 ± 0.163 |
| Random Forest | Approach 1 (DREAM) | Full microbiome | 0.660 ± 0.149 |
| glmnet | Approach 1 (DREAM) | Full microbiome | 0.649 ± 0.191 |
| Random Forest | Approach 2 (literature) | Full microbiome | 0.649 ± 0.160 |
| glmnet | Approach 3 (data-driven) | Full microbiome | 0.636 ± 0.221 |
| glmnet | Approach 1 (DREAM) | ANCOM taxa | 0.540 ± 0.098 |

Note the wide standard deviations. With 43 subjects and 14 preterm cases, each outer test fold holds
7–9 subjects (2–3 preterm), so per-fold AUROC is inherently unstable. This was true of the published
analysis as well; it is stated here because it bears on how much weight the ranking can carry.

### 1.3 The permutation test now yields a valid p-value

A reviewer asked for a permutation test showing the model's discrimination exceeds chance. An earlier
version of that test existed but required substantial revision (see 2.4). Updated result:

| | Value |
|---|---|
| Observed AUROC | **0.760** |
| Null distribution (999 subject-level label permutations, 998 valid) | **0.502 ± 0.121** |
| Permutations reaching the observed value | 18 / 998 |
| **One-sided permutation p-value** | **0.019** |

The null is centred on 0.502 — statistically indistinguishable from the 0.5 expected under no
signal (z = 0.60), and 48.6% of null values fall below 0.5. That is the check that the null is
legitimate.

**With the earlier statistic the same data would have yielded p = 0.139** — i.e. "not significant".
That outcome reflected an inflated null distribution rather than the model's actual discrimination.

---

## 2. What changed in the analysis, in detail

### 2.1 The CLR transformation now applies as intended (changes results)

Compositional data such as relative microbial abundances must be transformed before use in a linear
model. The code intended a centred log-ratio (CLR) transform, but added a constant of **0.65** to
relative abundances before taking logs. 0.65 is a pseudocount appropriate for raw *counts*; applied
to proportions (which sum to 1) it dominates the signal, and combined with the downstream
standardisation step it reduced to a per-feature `log(x + 0.65)` rescaling. The practical effect was
that **the intended compositional transformation was not taking effect.**

**Update.** Zeros are now replaced per taxon with `zCompositions::cmultRepl` (geometric Bayesian
multiplicative method), followed by a genuine per-sample CLR (each sample's transformed values sum
to zero). Zero-replacement levels are learned once on the full microbiome and are outcome-blind, so
they cannot leak outcome information into the folds.

**Effect:** reorders the model ranking (see 1.1).

### 2.2 One genus is now retained in differential-abundance analysis (changes results)

The absolute-count matrix used for ANCOM-BC2 was read in a way that renamed `Escherichia-Shigella`
to `Escherichia.Shigella` (R's default column-name sanitisation). Because downstream code matched
genus names literally, the renamed taxon was not picked up, and it dropped out of the candidate pool
without raising any message.

**Update.** The matrix is now read preserving exact names and matched by name. The ANCOM candidate
pool goes from 48 to **49** taxa, and the ANCOM-based results were re-derived.

### 2.3 Data domains are now split by name rather than column position (robustness)

Microbiome and clinical blocks were separated positionally (`data[, 1:99]` and `data[, 98:167]`).
Besides being fragile, the ranges overlapped: columns 98–99 (`index`, `id`) fell into *both*
blocks.

**Update.** The split is now driven by an explicit list of variables in `config/data_dictionary.csv`
(`role = microbiome`). No positional indexing anywhere in the pipeline. This also makes the pipeline
usable on other datasets, which positional indexing had prevented.

### 2.4 The permutation test was rebuilt (changes the reported p-value)

Four independent aspects were revised:

**(a) The null distribution was inflated.** Per-fold AUROC was computed with pROC's
`direction = "auto"`, which infers the orientation of the score from the data (it compares the two
groups' medians rather than being fixed a priori). Under randomly permuted labels this adapts to
noise: the null distribution centred on **0.68** instead of 0.5, and only 0.4% of null values fell
below 0.5. Because the null was nearly as high as the observed value, the test had almost no power.
The AUROC is now computed with a fixed orientation (higher predicted probability = higher predicted
risk), which is decided before seeing the labels. The updated null centres on 0.502.

*Scope note:* for the winning model this changes nothing about the reported performance — we
verified fold by fold that the fixed-direction and automatic-direction AUROCs are identical
(0.760 either way). The effect is on the *null*, not on this observed value. Whether the same holds
for the other 11 combinations has not yet been examined; the pipeline engine still reports
automatic-direction AUROC, and revisiting that is tracked as pending work.

**(b) The reference value is no longer hardcoded.** The script compared the null against a
hardcoded `0.813`, which predates updates 2.1 and 2.2. The observed value is now recomputed from the
current pipeline, using exactly the same code path and statistic as the null.

**(c) It now targets the reported model.** The script was pointed at Random Forest / Approach 3 /
full microbiome, which after these updates is no longer the best-performing combination (it is 6th,
at 0.680). The test now targets the model that is actually reported (see 1.1).

**(d) It now runs end to end.** The script had drifted out of sync with the analysis engine (missing
required arguments) and would have stopped on the first call. It also overwrote the full
572-permutation result with a 10-permutation trial run while still reporting "572" in the output, so
any p-value it produced would have rested on 10 permutations. Both behaviours are resolved: there is
now a single implementation that calls the current engine directly.

Two further points, relevant to whether the null is trustworthy:

- **Label permutation must reach the feature-selection step.** ANCOM-BC2 selects taxa using the
  outcome labels. If only the model's outcome column is permuted and the labels used by ANCOM are
  not, the taxa are still selected using the *true* labels and the null is contaminated (too high).
  All three places the outcome is consumed are now permuted together.
- **Degenerate folds are handled explicitly.** Outer folds are stratified using the true labels;
  once labels are permuted, a test fold can end up containing a single outcome class, for which
  AUROC is undefined. This happened in 111 of 999 permutations (11%). Those folds are now skipped
  and the permutation is scored on its remaining folds, rather than the whole permutation being
  discarded (which would have thrown away 11% of the null).

---

## 3. Reproducibility work (does not change results)

The analysis was migrated from a single monolithic R Markdown file to a research compendium:

- **Reusable engine in `R/`.** Pure functions (CLR, feature selection, threshold optimisation,
  nested CV) with no hidden global state, no `setwd()`, no absolute paths. Verified to reproduce the
  original numbers exactly at each migration step.
- **Configuration-driven.** Paths, variable roles, model specifications, and approach definitions
  live in `config/config.yml` and `config/data_dictionary.csv` instead of being hardcoded.
- **Automated tests** (`tests/testthat/`), including checks that the leakage-sensitive steps (CLR,
  ANCOM-BC2 selection, threshold optimisation) are fitted only on training data, plus an end-to-end
  test that permuted labels collapse performance to chance.
- **Pinned dependencies** (`renv.lock`) and **continuous integration** (GitHub Actions) that
  restores the exact package set, installs the package, and runs the full test suite — including
  the no-leakage test — on every push.
- **Synthetic example data** (`data/example/`) so the pipeline can be run end to end without access
  to the restricted microbiome data.

---

## 4. How to reproduce the updated results

The full cohort data are restricted (see `docs/DATA_ACCESS.md`); `data/example/` contains synthetic
data with the same structure.

```bash
# Full pipeline (12 model combinations)
Rscript -e 'rmarkdown::render("analysis/integrated_preterm_prediction_workflow.Rmd")'

# Permutation test — observed value + per-fold supplementary table (~40 s)
PTB_PERM_MODE=verify Rscript scripts/permutation_test.R

# Permutation test — null distribution (999 permutations, ~3.6 h on 8 cores)
PTB_PERM_MODE=run PTB_PERM_N=999 PTB_PERM_CORES=8 PTB_PERM_TAG=final999 Rscript scripts/permutation_test.R

# Figure, p-value and summary tables from the saved null (seconds)
PTB_PERM_MODE=report PTB_PERM_TAG=final999 Rscript scripts/permutation_test.R
```

The permutation seeds are fixed (`10001`–`10999`), so the p-value is reproducible; this was
confirmed across independent runs.

---

## 5. Limitations of the updated analysis

- **The p-value is conditional on model selection.** The permutation test is a valid test of the
  null for the specific pipeline configuration it was run on — labels are permuted and the entire
  nested CV, including per-fold feature selection, is re-run. What it does **not** account for is
  that this configuration was itself chosen as the best of 12 combinations evaluated on the same
  43 subjects. **p = 0.019 therefore quantifies "this particular model beats chance", not "the best
  of 12 models beats chance", and it carries no adjustment for multiplicity.** We report it
  unadjusted and state the conditioning rather than substituting an adjusted number we cannot
  justify: a naive Bonferroni adjustment (0.019 × 12 ≈ 0.23) would over-correct badly, because the
  12 combinations are not independent tests — they share the same subjects, the same outer folds
  and largely the same features, so the effective number of independent comparisons is well below
  12, and no defensible estimate of it is available at this sample size. The p-value should be read
  as supporting evidence for the reported model, not as a family-wise significance claim.
- **Small sample.** 43 subjects, 14 preterm, 110 longitudinal samples. Per-fold AUROC standard
  deviations are large (see 1.2) and the ranking among mid-table combinations should not be
  over-interpreted.
- **The engine still reports automatic-direction AUROC** for its performance tables. For the winning
  model this is verified to be identical to the fixed-direction value; it has not yet been examined
  for the other combinations.

---

## 6. Where things live

| Path | What it is |
|---|---|
| `README.md` | Current results, how to run the pipeline, how to use it on other data |
| `analysis/integrated_preterm_prediction_workflow.Rmd` | Main pipeline (all 12 combinations, figures) |
| `scripts/permutation_test.R` | Permutation test — the single implementation |
| `R/` | Reusable engine (CLR, feature selection, thresholds, nested CV) |
| `config/config.yml`, `config/data_dictionary.csv` | Configuration and data contract |
| `tests/testthat/` | Test suite, including no-leakage checks |
| `data/example/` | Synthetic example data |
| `vignettes/ptbpredict.Rmd` | Narrative walkthrough of a complete run |
| `CONTRIBUTING.md` | Invariants that must be preserved; how to propose a change |
| `results/permutation_test/` | Saved permutation-test artefacts (figure, null distribution, per-fold table) |
| `results/published_version/` | Archived report and figures from before these updates |

The previous permutation-test implementation (`analysis/permutation_test_nested_cv.Rmd`) was retired
rather than revised: it duplicated the analysis engine, and that duplication is precisely how it
drifted out of sync with it. Maintaining a single implementation is what the current
`scripts/permutation_test.R` provides. The retired file remains in the Git history.
