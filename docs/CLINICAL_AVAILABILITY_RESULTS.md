# Results with clinical variables known at the sample visit

**Last updated:** 2026-10-07

[`CLINICAL_VARIABLE_AVAILABILITY.md`](CLINICAL_VARIABLE_AVAILABILITY.md) classifies every clinical
variable by when its value is known relative to the visit at which the vaginal sample was taken,
and shows that 8 of the 12 model combinations used variables not known at the visit. The pipeline
now reads that classification: only the classes listed under `clinical_availability$allowed` in
[`config/config.yml`](../config/config.yml) can enter a model.

This document reports the 12 combinations re-run with clinical variables known at the visit only
(**main analysis**), compares them with the run that admitted every variable, and reports a
**sensitivity analysis** that adds back the variables known only after the visit. Variables
defined by the delivery itself (`outcome_defined`) enter neither analysis.

All values are the mean ± SD across the 5 outer test folds. Four test folds hold 9 subjects (3
preterm) and one holds 7 (2 preterm), so one preterm–term pair ranked the other way round moves a
fold's AUROC by 0.056 (0.100 in the smaller fold).

---

## 1. What was run

Three runs of [`scripts/run_baseline.R`](../scripts/run_baseline.R) on the cohort data. They differ
only in the availability classes admitted as clinical predictors:

| Run | Classes admitted | Approach 2 candidates | Approach 3 pool |
|---|---|---|---|
| Before | `at_visit`, `after_visit`, `outcome_defined` | 17 | 37 |
| **Main analysis** | `at_visit` | 12 (5 excluded) | 26 (11 excluded) |
| Sensitivity analysis | `at_visit`, `after_visit` | 16 (`rpm` excluded) | 32 (5 excluded) |

How the restriction is applied:

- **Approach 3** selects its variables inside each outer fold, so its screening pool is restricted
  before the fold loop. Everything downstream (screening, fold models, final model) sees only the
  admitted variables.
- **Approach 2** has a fixed list of 10 chosen from 17 literature candidates by completeness,
  collinearity and evidence ranking (the rules are unchanged). The candidates are restricted
  *before* that selection, so the 10 slots are filled from admitted candidates.
- **Approach 1** uses gestational age at the visit and maternal age, both `at_visit`. Its four
  combinations are therefore identical in the three runs, which serves as a check.

Everything else is the same in the three runs, and the comparison script
([`scripts/compare_availability_runs.R`](../scripts/compare_availability_runs.R)) verifies it:

- the outer folds (5, subject level, stratified by outcome, seed 123) and the subjects in each test
  fold;
- the count table used by ANCOM-BC2 (`genus_abs_2026-09-20.csv`) and the taxa it selected in each
  fold. ANCOM-BC2 adjusts only for maternal age and pre-pregnancy BMI, both `at_visit`, so the
  restriction does not reach it;
- the metrics: AUROC with the orientation fixed in advance (a higher predicted probability of
  preterm birth is the positive call), PR-AUC for the preterm class, and the Youden threshold chosen
  on the inner validation split with the same orientation (invariant 9 in
  [`CONTRIBUTING.md`](../CONTRIBUTING.md)).

The "before" run uses the current code with every class admitted. It reproduces the run made
before the pipeline read the classification exactly: every fold metric, prediction and selection is
identical. The differences below therefore come only from the restriction. No variable of a class
that was not admitted entered a model in any run.

## 2. Main analysis: clinical variables known at the visit

| Model | Clinical variable set | Microbiome input | AUROC | PRAUC | Sensitivity | Specificity |
|---|---|---|---|---|---|---|
| **Elastic net (glmnet)** | **3 — data-driven** | **ANCOM-BC2 taxa** | **0.758 ± 0.240** | 0.714 ± 0.300 | 0.600 ± 0.435 | 0.520 ± 0.315 |
| Elastic net (glmnet) | 3 — data-driven | Full microbiome | 0.700 ± 0.178 | 0.629 ± 0.245 | 0.500 ± 0.289 | 0.707 ± 0.404 |
| Random Forest | 3 — data-driven | ANCOM-BC2 taxa | 0.591 ± 0.196 | 0.586 ± 0.243 | 0.600 ± 0.279 | 0.707 ± 0.283 |
| Random Forest | 3 — data-driven | Full microbiome | 0.573 ± 0.240 | 0.515 ± 0.213 | 0.333 ± 0.236 | 0.607 ± 0.266 |
| Random Forest | 2 — literature | ANCOM-BC2 taxa | 0.569 ± 0.199 | 0.406 ± 0.240 | 0.333 ± 0.333 | 0.787 ± 0.157 |
| Random Forest | 2 — literature | Full microbiome | 0.562 ± 0.236 | 0.397 ± 0.141 | 0.600 ± 0.279 | 0.707 ± 0.283 |
| Random Forest | 1 — DREAM (minimal) | ANCOM-BC2 taxa | 0.560 ± 0.268 | 0.480 ± 0.181 | 0.533 ± 0.298 | 0.507 ± 0.325 |
| Random Forest | 1 — DREAM (minimal) | Full microbiome | 0.496 ± 0.233 | 0.423 ± 0.149 | 0.267 ± 0.365 | 0.647 ± 0.282 |
| Elastic net (glmnet) | 1 — DREAM (minimal) | Full microbiome | 0.471 ± 0.251 | 0.453 ± 0.204 | 0.300 ± 0.298 | 0.760 ± 0.434 |
| Elastic net (glmnet) | 2 — literature | Full microbiome | 0.396 ± 0.254 | 0.393 ± 0.271 | 0.333 ± 0.236 | 0.593 ± 0.340 |
| Elastic net (glmnet) | 2 — literature | ANCOM-BC2 taxa | 0.389 ± 0.304 | 0.406 ± 0.281 | 0.133 ± 0.183 | 0.827 ± 0.205 |
| Elastic net (glmnet) | 1 — DREAM (minimal) | ANCOM-BC2 taxa | 0.262 ± 0.086 | 0.225 ± 0.036 | 0.067 ± 0.149 | 0.853 ± 0.202 |

The best combination is elastic net with the data-driven clinical set and ANCOM-BC2 taxa, AUROC
0.758 ± 0.240, the same combination that ranked first before the restriction. Only the two elastic
net combinations of Approach 3 have a mean AUROC above 0.62; seven of the 12 are above 0.5. Five
have a mean AUROC below 0.5: with a fixed orientation, a model that ranks preterm subjects below term
subjects is reported as such rather than reversed. The PR-AUC reference line, the share of preterm
subjects in a test fold, is 0.29–0.33.

## 3. Before and after

Means across the 5 outer folds; "before" admits every class, "after" is the main analysis.

| Model | Clinical variable set | Microbiome input | AUROC | ΔAUROC | PRAUC | Sensitivity | Specificity |
|---|---|---|---|---|---|---|---|
| Elastic net (glmnet) | 1 — DREAM (minimal) | ANCOM-BC2 taxa | 0.262 → 0.262 | +0.000 | 0.225 → 0.225 | 0.067 → 0.067 | 0.853 → 0.853 |
| Random Forest | 1 — DREAM (minimal) | ANCOM-BC2 taxa | 0.560 → 0.560 | +0.000 | 0.480 → 0.480 | 0.533 → 0.533 | 0.507 → 0.507 |
| Elastic net (glmnet) | 1 — DREAM (minimal) | Full microbiome | 0.471 → 0.471 | +0.000 | 0.453 → 0.453 | 0.300 → 0.300 | 0.760 → 0.760 |
| Random Forest | 1 — DREAM (minimal) | Full microbiome | 0.496 → 0.496 | +0.000 | 0.423 → 0.423 | 0.267 → 0.267 | 0.647 → 0.647 |
| Elastic net (glmnet) | 2 — literature | ANCOM-BC2 taxa | 0.278 → 0.389 | +0.111 | 0.235 → 0.406 | 0.067 → 0.133 | 0.633 → 0.827 |
| Random Forest | 2 — literature | ANCOM-BC2 taxa | 0.536 → 0.569 | +0.033 | 0.398 → 0.406 | 0.200 → 0.333 | 0.753 → 0.787 |
| Elastic net (glmnet) | 2 — literature | Full microbiome | 0.384 → 0.396 | +0.011 | 0.353 → 0.393 | 0.200 → 0.333 | 0.767 → 0.593 |
| Random Forest | 2 — literature | Full microbiome | 0.529 → 0.562 | +0.033 | 0.384 → 0.397 | 0.400 → 0.600 | 0.713 → 0.707 |
| Elastic net (glmnet) | 3 — data-driven | ANCOM-BC2 taxa | 0.691 → 0.758 | +0.067 | 0.649 → 0.714 | 0.767 → 0.600 | 0.413 → 0.520 |
| Random Forest | 3 — data-driven | ANCOM-BC2 taxa | 0.622 → 0.591 | −0.031 | 0.590 → 0.586 | 0.467 → 0.600 | 0.560 → 0.707 |
| Elastic net (glmnet) | 3 — data-driven | Full microbiome | 0.658 → 0.700 | +0.042 | 0.614 → 0.629 | 0.567 → 0.500 | 0.667 → 0.707 |
| Random Forest | 3 — data-driven | Full microbiome | 0.542 → 0.573 | +0.031 | 0.508 → 0.515 | 0.400 → 0.333 | 0.540 → 0.607 |

Over the 12 combinations, the mean AUROC goes from 0.502 to 0.527 and the mean PR-AUC from 0.443 to
0.469; seven combinations are above 0.5 in both runs.

- **Approach 1** does not change, as expected.
- **Approach 2**: AUROC rises in its four combinations (+0.011 to +0.111). The two elastic net
  combinations remain below 0.5.
- **Approach 3**: AUROC rises in three combinations (+0.031 to +0.067) and falls in one (Random
  Forest with ANCOM-BC2 taxa, −0.031). In the best combination the change comes from two outer
  folds: fold 1 (0.72 → 0.94), where the screening had selected `oligohidramnios`,
  `complication_count` and `any_complication`, and fold 4 (0.67 → 0.78), where it had selected
  `compvaginf`. The other three folds had also selected `oligohidramnios` or `compvaginf`, and their
  AUROC is the same as before.

Restricting the models to variables known at the visit did not lower the mean AUROC of any approach
(Approach 1 unchanged, Approach 2 +0.047, Approach 3 +0.027): in this cohort, the performance
reported before did not depend on information that is not available at the visit. The changes are
small compared with the fold-to-fold SDs (0.09–0.30), so they are not evidence that the restricted
models discriminate better.

## 4. Sensitivity analysis: adding the variables known after the visit

The sensitivity analysis admits `at_visit` and `after_visit`. Compared with the main analysis, it
adds back the undated obstetric complications (`comppreeclam`, `compvaginf`, `diabetes_gest`,
`oligohidramnios`, `rciu`) and `sex_baby`; `rpm`, `rpm_preterm` and the complication counts
(`outcome_defined`) remain excluded.

| Model | Clinical variable set | Microbiome input | AUROC | PRAUC | Sensitivity | Specificity |
|---|---|---|---|---|---|---|
| Elastic net (glmnet) | 3 — data-driven | ANCOM-BC2 taxa | 0.658 ± 0.220 | 0.584 ± 0.315 | 0.633 ± 0.415 | 0.480 ± 0.315 |
| Elastic net (glmnet) | 3 — data-driven | Full microbiome | 0.636 ± 0.207 | 0.604 ± 0.254 | 0.567 ± 0.149 | 0.667 ± 0.471 |
| Random Forest | 3 — data-driven | ANCOM-BC2 taxa | 0.589 ± 0.178 | 0.556 ± 0.217 | 0.467 ± 0.298 | 0.627 ± 0.357 |
| Random Forest | 1 — DREAM (minimal) | ANCOM-BC2 taxa | 0.560 ± 0.268 | 0.480 ± 0.181 | 0.533 ± 0.298 | 0.507 ± 0.325 |
| Random Forest | 2 — literature | ANCOM-BC2 taxa | 0.558 ± 0.213 | 0.404 ± 0.244 | 0.333 ± 0.236 | 0.687 ± 0.222 |
| Random Forest | 3 — data-driven | Full microbiome | 0.542 ± 0.280 | 0.508 ± 0.220 | 0.400 ± 0.279 | 0.507 ± 0.303 |
| Random Forest | 2 — literature | Full microbiome | 0.518 ± 0.225 | 0.381 ± 0.147 | 0.600 ± 0.279 | 0.640 ± 0.285 |
| Random Forest | 1 — DREAM (minimal) | Full microbiome | 0.496 ± 0.233 | 0.423 ± 0.149 | 0.267 ± 0.365 | 0.647 ± 0.282 |
| Elastic net (glmnet) | 1 — DREAM (minimal) | Full microbiome | 0.471 ± 0.251 | 0.453 ± 0.204 | 0.300 ± 0.298 | 0.760 ± 0.434 |
| Elastic net (glmnet) | 2 — literature | Full microbiome | 0.384 ± 0.158 | 0.356 ± 0.182 | 0.133 ± 0.183 | 0.727 ± 0.344 |
| Elastic net (glmnet) | 1 — DREAM (minimal) | ANCOM-BC2 taxa | 0.262 ± 0.086 | 0.225 ± 0.036 | 0.067 ± 0.149 | 0.853 ± 0.202 |
| Elastic net (glmnet) | 2 — literature | ANCOM-BC2 taxa | 0.256 ± 0.174 | 0.227 ± 0.048 | 0.133 ± 0.183 | 0.587 ± 0.382 |

Main analysis → sensitivity analysis (means):

| Model | Clinical variable set | Microbiome input | AUROC | ΔAUROC | PRAUC | Sensitivity | Specificity |
|---|---|---|---|---|---|---|---|
| Elastic net (glmnet) | 1 — DREAM (minimal) | ANCOM-BC2 taxa | 0.262 → 0.262 | +0.000 | 0.225 → 0.225 | 0.067 → 0.067 | 0.853 → 0.853 |
| Random Forest | 1 — DREAM (minimal) | ANCOM-BC2 taxa | 0.560 → 0.560 | +0.000 | 0.480 → 0.480 | 0.533 → 0.533 | 0.507 → 0.507 |
| Elastic net (glmnet) | 1 — DREAM (minimal) | Full microbiome | 0.471 → 0.471 | +0.000 | 0.453 → 0.453 | 0.300 → 0.300 | 0.760 → 0.760 |
| Random Forest | 1 — DREAM (minimal) | Full microbiome | 0.496 → 0.496 | +0.000 | 0.423 → 0.423 | 0.267 → 0.267 | 0.647 → 0.647 |
| Elastic net (glmnet) | 2 — literature | ANCOM-BC2 taxa | 0.389 → 0.256 | −0.133 | 0.406 → 0.227 | 0.133 → 0.133 | 0.827 → 0.587 |
| Random Forest | 2 — literature | ANCOM-BC2 taxa | 0.569 → 0.558 | −0.011 | 0.406 → 0.404 | 0.333 → 0.333 | 0.787 → 0.687 |
| Elastic net (glmnet) | 2 — literature | Full microbiome | 0.396 → 0.384 | −0.011 | 0.393 → 0.356 | 0.333 → 0.133 | 0.593 → 0.727 |
| Random Forest | 2 — literature | Full microbiome | 0.562 → 0.518 | −0.044 | 0.397 → 0.381 | 0.600 → 0.600 | 0.707 → 0.640 |
| Elastic net (glmnet) | 3 — data-driven | ANCOM-BC2 taxa | 0.758 → 0.658 | −0.100 | 0.714 → 0.584 | 0.600 → 0.633 | 0.520 → 0.480 |
| Random Forest | 3 — data-driven | ANCOM-BC2 taxa | 0.591 → 0.589 | −0.002 | 0.586 → 0.556 | 0.600 → 0.467 | 0.707 → 0.627 |
| Elastic net (glmnet) | 3 — data-driven | Full microbiome | 0.700 → 0.636 | −0.064 | 0.629 → 0.604 | 0.500 → 0.567 | 0.707 → 0.667 |
| Random Forest | 3 — data-driven | Full microbiome | 0.573 → 0.542 | −0.031 | 0.515 → 0.508 | 0.333 → 0.400 | 0.607 → 0.507 |

Adding the variables known after the visit lowers the AUROC of all eight Approach 2 and 3
combinations (−0.002 to −0.133); over the 12 combinations the mean AUROC goes from 0.527 to 0.494
and the mean PR-AUC from 0.469 to 0.433. The best combination goes from 0.758 to 0.658. Adding these
variables back does not improve on the main analysis. They are rare in this cohort: `oligohidramnios`
equals 1 in 3 participants (2 preterm), `compvaginf` in 2 (1 preterm), `comppreeclam` and `rciu` in
the same single preterm participant, and `diabetes_gest` in none
([`CLINICAL_VARIABLE_AVAILABILITY.md`](CLINICAL_VARIABLE_AVAILABILITY.md) §3). A variable that
equals 1 in two or three participants can pass the screening in the outer-training subjects without
carrying over to the test fold.

In the Approach 2 list of the sensitivity analysis, the slot that `rpm` held goes to
`diabetes_gest`, which is 0 for every participant; the preprocessing removes it as a zero-variance
column, so those models use nine clinical variables.

## 5. Which variables left each approach

### Approach 2: the final list of 10

The list is the same in the 5 outer folds and in the four Approach 2 combinations.

| Variable | Class | Before | Main | Sensitivity |
|---|---|---|---|---|
| `sdg_visita` | at_visit | ✓ | ✓ | ✓ |
| `edad_cronologicamujer` | at_visit | ✓ | ✓ | ✓ |
| `imc_pregestacional` | at_visit | ✓ | ✓ | ✓ |
| `peso_pregestacional_kg` | at_visit | ✓ | ✓ | ✓ |
| `hemoglobin_alti_adj` | at_visit | ✓ | ✓ | ✓ |
| `folic_ac_correg` | at_visit | ✓ | ✓ | ✓ |
| `comp1tribleed` | at_visit | ✓ | ✓ | ✓ |
| `comppreeclam` | after_visit | ✓ | — | ✓ |
| `rpm` | outcome_defined | ✓ | — | — |
| `oligohidramnios` | after_visit | ✓ | — | ✓ |
| `anemia_visita` | at_visit | — | ✓ | — |
| `ses_risk_score` | at_visit | — | ✓ | — |
| `dietsuppl3mon` | at_visit | — | ✓ | — |
| `diabetes_gest` | after_visit | — | — | ✓ |

In the main analysis, the three variables not known at the visit are replaced by the next
`at_visit` candidates in the evidence ranking: anemia at the visit, the socioeconomic risk score and
dietary supplement use.

### Approach 3: the screening pool and the per-fold selection

Excluded from the screening pool:

- **Main analysis (11):** `comppreeclam`, `compvaginf`, `diabetes_gest`, `oligohidramnios`, `rciu`,
  `sex_baby` (`after_visit`); `rpm`, `rpm_preterm`, `complication_count`, `any_complication`,
  `multiple_complications` (`outcome_defined`).
- **Sensitivity analysis (5):** the five `outcome_defined` variables above.

Number of outer folds (of 5) in which each variable was selected. The screening uses the clinical
data only, so the selection is the same in the four Approach 3 combinations.

| Variable | Class | Before | Main | Sensitivity |
|---|---|---|---|---|
| `extreme_bmi` | at_visit | 4 | 4 | 4 |
| `comp1tribleed` | at_visit | 4 | 4 | 4 |
| `workouthome` | at_visit | 4 | 4 | 4 |
| `extreme_age` | at_visit | 4 | 4 | 4 |
| `obese` | at_visit | 3 | 3 | 3 |
| `ses_risk_score` | at_visit | 3 | 3 | 3 |
| `maritalstat` | at_visit | 3 | 3 | 3 |
| `age_risk_category` | at_visit | 1 | 3 | 1 |
| `imc_pregestacional` | at_visit | 2 | 2 | 2 |
| `peso_kg` | at_visit | 2 | 2 | 2 |
| `dietsuppl3mon` | at_visit | 2 | 2 | 2 |
| `hemoglobin_alti_adj` | at_visit | 1 | 1 | 1 |
| `peso_pregestacional_kg` | at_visit | 1 | 1 | 1 |
| `sdg_visita` | at_visit | 1 | 1 | 1 |
| `folic_ac_correg` | at_visit | 1 | 1 | 1 |
| `hemoglobin_g_dl` | at_visit | 1 | 1 | 1 |
| `anemia_visita` | at_visit | 0 | 1 | 1 |
| `oligohidramnios` | after_visit | 3 | 0 | 3 |
| `compvaginf` | after_visit | 3 | 0 | 3 |
| `sex_baby` | after_visit | 0 | 0 | 1 |
| `any_complication` | outcome_defined | 1 | 0 | 0 |
| `complication_count` | outcome_defined | 1 | 0 | 0 |

Variables per fold: before 10, 8, 8, 9, 10; main analysis 9, 7, 7, 8, 9; sensitivity analysis 10, 8,
8, 9, 10. In the main analysis, `age_risk_category` (folds 1 and 5) and `anemia_visita` (fold 1)
entered where variables were excluded; no other `at_visit` variable met the screening criteria, so
each fold uses one variable fewer than before. In the sensitivity analysis the selection is the same
as before except in fold 1, where `anemia_visita` and `sex_baby` take the places of
`complication_count` and `any_complication`.

## 6. Not covered here

- **The permutation test** reported in the [README](../README.md#permutation-test) (observed AUROC
  0.760, p = 0.019) was computed on an earlier state of the pipeline. It has not been repeated on
  these results.
- **The participant whose pregnancy ended at 20+0 weeks** is retained as preterm
  ([`CLINICAL_VARIABLE_AVAILABILITY.md`](CLINICAL_VARIABLE_AVAILABILITY.md) §6). The analysis without
  this participant will be repeated on these results and reported separately.
- **The README results table** predates this restriction and other updates; it will be updated
  together with [`UPDATE_SINCE_PUBLICATION.md`](UPDATE_SINCE_PUBLICATION.md).

## 7. Files and reproducing these numbers

[`results/clinical_availability/`](../results/clinical_availability/) holds the numbers of the three
runs (`analysis` = `before`, `main`, `sensitivity`):

- `metrics_by_combination.csv` — one row per run and combination: mean and SD across the outer folds
  of AUROC, PR-AUC, sensitivity, specificity and balanced accuracy.
- `clinical_variables_by_fold.csv` — the clinical variables each combination used: the fixed list
  for Approaches 1 and 2 (`fold` = `all`) and the per-fold selection for Approach 3, with the
  availability class of each variable.

The runs need the cohort data in `data/raw/`. The admitted classes come from the config; the
environment variable `PTB_AVAILABILITY` overrides them for one run (admitting `outcome_defined`
prints a warning, because it is only meant to reproduce the earlier analysis):

```bash
PTB_AVAILABILITY=at_visit,after_visit,outcome_defined PTB_RUN_TAG=availctrl_2026-10-07 Rscript scripts/run_baseline.R
```

```bash
PTB_RUN_TAG=atvisit_2026-10-07 Rscript scripts/run_baseline.R
```

```bash
PTB_AVAILABILITY=at_visit,after_visit PTB_RUN_TAG=aftervisit_2026-10-07 Rscript scripts/run_baseline.R
```

```bash
PTB_RESULTS_DIR=results/clinical_availability Rscript scripts/compare_availability_runs.R
```

Each run took 4–8 minutes (Windows 11, R 4.4.2). The comparison script prints the checks of
section 1 and the tables of this document.
