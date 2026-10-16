# Socioeconomic indicators: education level read as a code

**Last updated:** 2026-10-16

The analysis derives three socioeconomic indicators from the clinical table (chunk
`feature_engineering` of the analysis notebook, now in
[`R/clinical_features.R`](../R/clinical_features.R)):

| Indicator | Definition ([`config/data_dictionary.csv`](../config/data_dictionary.csv)) |
|---|---|
| `low_education` | 1 if the highest educational level is secondary education or lower |
| `unmarried` | 1 if the marital status is other than married |
| `ses_risk_score` | `low_education` + `unmarried` (0 = lower risk, 2 = higher risk) |

This document reports a methodological update to `low_education`, and therefore to
`ses_risk_score`, and its effect on the 12 combinations of the main analysis (clinical variables
known at the sample visit, [`CLINICAL_AVAILABILITY_RESULTS.md`](CLINICAL_AVAILABILITY_RESULTS.md)).

All values are the mean ± SD across the 5 outer test folds. Four test folds hold 9 subjects (3
preterm) and one holds 7 (2 preterm), so one preterm–term pair ranked the other way round moves a
fold's AUROC by 0.056 (0.100 in the smaller fold).

---

## 1. What changed

In the cohort table, `nivel_academico` is stored as text with the code first:

| Value in the table | Samples | Participants |
|---|---|---|
| `1.- Primary education or less` | 2 | 2 |
| `2.- Secondary education` | 44 | 19 |
| `3.- Vocational/Technical` | 60 | 20 |
| `4.- Associates degree or higher` | 2 | 1 |
| missing | 2 | 1 |

The level is constant across the visits of each participant. The indicator was computed as
`nivel_academico <= 2`. With a text column, R compares the two values as strings, and
`"2.- Secondary education" <= "2"` is `FALSE`, so only primary education was flagged. The
synthetic example data store the bare integer code, which is why the example-based tests did not
show it.

`derive_ses_features()` now reads the leading code before comparing, accepts integer codes as
well, and stops if a value has no leading code (instead of returning a missing indicator).
`tests/testthat/test-clinical-features.R` checks the indicators on the cohort's text values; it
fails with the previous comparison.

| | Before | After |
|---|---|---|
| `low_education` = 1: samples | 2 | 46 |
| `low_education` = 1: participants (preterm / term) | 2 (2 / 0) | 21 (8 / 13) |
| `ses_risk_score` 0 / 1 / 2, samples | 59 / 48 / 1 | 37 / 48 / 23 |
| `ses_risk_score` 0 / 1 / 2, participants | 19 / 22 / 1 | 11 / 19 / 12 |
| Correlation of `ses_risk_score` with `unmarried`, across samples | 0.97 | 0.74 |

Before the update, `ses_risk_score` was therefore almost a copy of `unmarried`. `unmarried` is
unchanged. On the example data the indicators are identical before and after.

## 2. What was run

Two runs of [`scripts/run_baseline.R`](../scripts/run_baseline.R) on the cohort data, both with the
main analysis's admitted availability class (`at_visit`). "Before" is the run reported in
[`CLINICAL_AVAILABILITY_RESULTS.md`](CLINICAL_AVAILABILITY_RESULTS.md); "after" differs from it only
in the derivation of the indicators. The comparison script
([`scripts/compare_ses_runs.R`](../scripts/compare_ses_runs.R)) verifies that the two runs share:

- the outer folds (5, subject level, stratified by outcome, seed 123) and the subjects in each test
  fold;
- the taxa ANCOM-BC2 selected in each fold. ANCOM-BC2 adjusts only for maternal age and
  pre-pregnancy BMI, so the indicators do not reach it;
- the results of the four Approach 1 combinations, which use gestational age at the visit and
  maternal age only. They are identical in every fold;
- the admitted availability class: no variable of another class entered a model in either run.

The metrics are those of the main analysis: AUROC with the orientation fixed in advance (a higher
predicted probability of preterm birth is the positive call), PR-AUC for the preterm class, and the
Youden threshold chosen on the inner validation split with the same orientation.

## 3. Before and after

Which combinations use the indicators:

- **Approach 2** keeps the same list of 10 variables, `ses_risk_score` among them; only the values of
  `ses_risk_score` change.
- **Approach 3** selects its variables inside each outer fold, and the selection is the same in its
  four combinations. `ses_risk_score` is selected in 3 of the 5 folds in both runs; `low_education`
  goes from 0 to 4 folds (section 4).
- **Approach 1** does not use them.

Means across the 5 outer folds, ordered by the AUROC after the update:

| Model | Clinical variable set | Microbiome input | AUROC | ΔAUROC | PR-AUC | Sensitivity | Specificity |
|---|---|---|---|---|---|---|---|
| Elastic net (glmnet) | 3 — data-driven | ANCOM-BC2 taxa | 0.758 → 0.738 | −0.020 | 0.714 → 0.716 | 0.600 → 0.467 | 0.520 → 0.787 |
| Elastic net (glmnet) | 3 — data-driven | Full microbiome | 0.700 → 0.669 | −0.031 | 0.629 → 0.627 | 0.500 → 0.433 | 0.707 → 0.740 |
| Random Forest | 3 — data-driven | ANCOM-BC2 taxa | 0.591 → 0.580 | −0.011 | 0.586 → 0.579 | 0.600 → 0.667 | 0.707 → 0.647 |
| Random Forest | 2 — literature | Full microbiome | 0.562 → 0.562 | 0 | 0.397 → 0.404 | 0.600 → 0.600 | 0.707 → 0.673 |
| Random Forest | 1 — DREAM (minimal) | ANCOM-BC2 taxa | 0.560 → 0.560 | 0 | 0.480 → 0.480 | 0.533 → 0.533 | 0.507 → 0.507 |
| Random Forest | 2 — literature | ANCOM-BC2 taxa | 0.569 → 0.558 | −0.011 | 0.406 → 0.395 | 0.333 → 0.333 | 0.787 → 0.787 |
| Random Forest | 3 — data-driven | Full microbiome | 0.573 → 0.531 | −0.042 | 0.515 → 0.481 | 0.333 → 0.333 | 0.607 → 0.707 |
| Random Forest | 1 — DREAM (minimal) | Full microbiome | 0.496 → 0.496 | 0 | 0.423 → 0.423 | 0.267 → 0.267 | 0.647 → 0.647 |
| Elastic net (glmnet) | 1 — DREAM (minimal) | Full microbiome | 0.471 → 0.471 | 0 | 0.453 → 0.453 | 0.300 → 0.300 | 0.760 → 0.760 |
| Elastic net (glmnet) | 2 — literature | Full microbiome | 0.396 → 0.351 | −0.044 | 0.393 → 0.350 | 0.333 → 0.200 | 0.593 → 0.693 |
| Elastic net (glmnet) | 2 — literature | ANCOM-BC2 taxa | 0.389 → 0.344 | −0.044 | 0.406 → 0.350 | 0.133 → 0.067 | 0.827 → 0.760 |
| Elastic net (glmnet) | 1 — DREAM (minimal) | ANCOM-BC2 taxa | 0.262 → 0.262 | 0 | 0.225 → 0.225 | 0.067 → 0.067 | 0.853 → 0.853 |

The same table after the update, with the SD across folds:

| Model | Clinical variable set | Microbiome input | AUROC | PR-AUC | Sensitivity | Specificity |
|---|---|---|---|---|---|---|
| **Elastic net (glmnet)** | **3 — data-driven** | **ANCOM-BC2 taxa** | **0.738 ± 0.270** | 0.716 ± 0.297 | 0.467 ± 0.298 | 0.787 ± 0.307 |
| Elastic net (glmnet) | 3 — data-driven | Full microbiome | 0.669 ± 0.184 | 0.627 ± 0.235 | 0.433 ± 0.279 | 0.740 ± 0.371 |
| Random Forest | 3 — data-driven | ANCOM-BC2 taxa | 0.580 ± 0.200 | 0.579 ± 0.244 | 0.667 ± 0.333 | 0.647 ± 0.282 |
| Random Forest | 2 — literature | Full microbiome | 0.562 ± 0.236 | 0.404 ± 0.139 | 0.600 ± 0.279 | 0.673 ± 0.274 |
| Random Forest | 1 — DREAM (minimal) | ANCOM-BC2 taxa | 0.560 ± 0.268 | 0.480 ± 0.181 | 0.533 ± 0.298 | 0.507 ± 0.325 |
| Random Forest | 2 — literature | ANCOM-BC2 taxa | 0.558 ± 0.190 | 0.395 ± 0.242 | 0.333 ± 0.333 | 0.787 ± 0.157 |
| Random Forest | 3 — data-driven | Full microbiome | 0.531 ± 0.284 | 0.481 ± 0.222 | 0.333 ± 0.236 | 0.707 ± 0.283 |
| Random Forest | 1 — DREAM (minimal) | Full microbiome | 0.496 ± 0.233 | 0.423 ± 0.149 | 0.267 ± 0.365 | 0.647 ± 0.282 |
| Elastic net (glmnet) | 1 — DREAM (minimal) | Full microbiome | 0.471 ± 0.251 | 0.453 ± 0.204 | 0.300 ± 0.298 | 0.760 ± 0.434 |
| Elastic net (glmnet) | 2 — literature | Full microbiome | 0.351 ± 0.221 | 0.350 ± 0.207 | 0.200 ± 0.183 | 0.693 ± 0.317 |
| Elastic net (glmnet) | 2 — literature | ANCOM-BC2 taxa | 0.344 ± 0.250 | 0.350 ± 0.201 | 0.067 ± 0.149 | 0.760 ± 0.188 |
| Elastic net (glmnet) | 1 — DREAM (minimal) | ANCOM-BC2 taxa | 0.262 ± 0.086 | 0.225 ± 0.036 | 0.067 ± 0.149 | 0.853 ± 0.202 |

- The highest combination is the same, elastic net with the data-driven clinical set and ANCOM-BC2
  taxa: **0.758 → 0.738 ± 0.270**. By fold: 0.944, 0.667, 1.000, 0.778, 0.400 → 0.944, 0.667,
  0.944, 0.833, 0.300, that is, one pair ranked differently in each of folds 3, 4 and 5. Its
  sensitivity goes from 0.600 to 0.467 and its specificity from 0.520 to 0.787.
- The eight Approach 2 and 3 combinations change; seven lose mean AUROC (−0.011 to −0.044) and one
  keeps it (random forest, Approach 2, full microbiome, whose folds change in opposite directions).
  Mean AUROC by approach: Approach 2 0.479 → 0.454, Approach 3 0.656 → 0.629, Approach 1 0.447
  unchanged. Over the 12 combinations the mean AUROC goes from 0.527 to 0.510 and the mean PR-AUC
  from 0.469 to 0.457; seven combinations have a mean AUROC above 0.5 in both runs.
- The changes amount to one or two pairs per fold and are small next to the fold-to-fold SDs
  (0.18–0.30 in Approaches 2 and 3). They are not evidence that the indicator, as defined, lowers
  or raises discrimination; they show that the reported numbers depended on it only slightly.

## 4. Which variables changed in Approach 3

Folds (of 5) in which each variable whose selection changed was chosen; the other selected
variables are the same in both runs (`results/socioeconomic_indicators/clinical_variables_by_fold.csv`):

| Variable | Before | After |
|---|---|---|
| `low_education` | 0 | 4 (folds 2–5) |
| `unmarried` | 0 | 3 (folds 3–5) |
| `ses_risk_score` | 3 | 3 (folds 3–5) |
| `age_risk_category` | 3 | 2 (leaves fold 5) |

Variables per fold go from 9, 7, 7, 8, 9 to 9, 8, 9, 10, 10. `low_education` flagged 2
participants before the update, both preterm, and 21 after. `unmarried` was not selected in any
fold before and is selected in three after. This is consistent with the screening's collinearity
rule, which keeps only one of two candidates whose correlation across the outer-training samples
exceeds 0.95: before the update, `ses_risk_score` and `unmarried` correlated at 0.957–0.978 in the
outer-training samples of the five folds; after it, at 0.684–0.785. Which of the two the rule kept
in each fold was not traced.

## 5. Not covered here

- **The sensitivity analysis** of [`CLINICAL_AVAILABILITY_RESULTS.md`](CLINICAL_AVAILABILITY_RESULTS.md)
  (adding back the variables known only after the visit) used the previous derivation and has not
  been repeated.
- **The permutation test** will use the updated principal model: observed AUROC 0.738.
- **`nivel_academico` itself** is a candidate in Approaches 2 and 3 and was not selected in either
  run; this update does not change it. The data dictionary's description of its codes 3 and 4 is
  reviewed separately, together with other descriptions that do not match the data.

## 6. Files and reproducing these numbers

[`results/socioeconomic_indicators/`](../results/socioeconomic_indicators/) holds the numbers of
the two runs (`run` = `before`, `after`):

- `metrics_by_combination.csv` — one row per run and combination: mean and SD across the outer folds
  of AUROC, PR-AUC, sensitivity, specificity and balanced accuracy, and which indicators the
  combination used.
- `clinical_variables_by_fold.csv` — the clinical variables each combination used: the fixed list
  for Approaches 1 and 2 (`fold` = `all`) and the per-fold selection for Approach 3.

The runs need the cohort data in `data/raw/`. The "before" run is `atvisit_2026-10-07`, made with
the code of [`CLINICAL_AVAILABILITY_RESULTS.md`](CLINICAL_AVAILABILITY_RESULTS.md) §7. The "after"
run, with the current code:

```bash
PTB_RUN_TAG=sesfix_2026-10-16 Rscript scripts/run_baseline.R
```

```bash
PTB_BEFORE_TAG=atvisit_2026-10-07 PTB_AFTER_TAG=sesfix_2026-10-16 PTB_RESULTS_DIR=results/socioeconomic_indicators Rscript scripts/compare_ses_runs.R
```

The run took 3.9 minutes (Windows 11, R 4.4.2). The comparison script prints the checks of
section 2 and the tables of this document.
