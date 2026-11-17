# Clinical-only arm: what the microbiome features add

**Last updated:** 2026-10-16

The 12 combinations of the main analysis (2 models × 3 clinical variable sets × 2 microbiome inputs,
with clinical variables known at the sample visit:
[`CLINICAL_AVAILABILITY_RESULTS.md`](CLINICAL_AVAILABILITY_RESULTS.md),
[`SOCIOECONOMIC_INDICATORS.md`](SOCIOECONOMIC_INDICATORS.md)) always combine clinical variables with
microbiome features, so on their own they cannot show what the microbiome adds. This document adds
the **clinical-only arm**: the same 2 models × 3 clinical variable sets, run with the same engine on
the same outer folds, with the same seeds and the same clinical variables, and no microbiome
feature (no taxon and no Shannon diversity). The 12 combinations with microbiome are the current
main-analysis results, reproduced in the same run (section 5).

All values are the mean ± SD across the 5 outer test folds: AUROC with the orientation fixed in
advance (preterm = higher predicted probability) and PR-AUC with preterm birth as the event. Four
test folds hold 9 subjects (3 preterm) and one holds 7 (2 preterm), so one preterm–term pair ranked
the other way round moves a fold's AUROC by 0.056 (0.100 in the smaller fold). A classifier without
discrimination has an expected PR-AUC equal to the preterm prevalence: 0.326 (14 of 43), 0.324
averaged over the folds.

The engine seeds every fold by its number, whatever the microbiome input, so both arms see the same
subjects in inner training, inner validation and outer test: the differences by fold are paired.

---

## 1. Approach 1: gestational age at the visit and maternal age

Approach 1 uses two clinical variables, gestational age at the sample visit (`sdg_visita`) and
maternal age (`edad_cronologicamujer`), mirroring the minimal adjustment of the DREAM Preterm Birth
Prediction Challenge (Golob et al. 2024). Its clinical part is the only one of the three sets that
can be carried to the DREAM challenge data, so this is the comparison that bears on an external
validation of the combined model.

| Model | Input | AUROC | PR-AUC |
|---|---|---|---|
| Elastic net (glmnet) | **Clinical only** | **0.378 ± 0.240** | **0.410 ± 0.191** |
| Elastic net (glmnet) | + ANCOM-BC2 taxa | 0.262 ± 0.086 | 0.225 ± 0.036 |
| Elastic net (glmnet) | + full microbiome | 0.471 ± 0.251 | 0.453 ± 0.204 |
| Random Forest | **Clinical only** | **0.600 ± 0.114** | **0.517 ± 0.137** |
| Random Forest | + ANCOM-BC2 taxa | 0.560 ± 0.268 | 0.480 ± 0.181 |
| Random Forest | + full microbiome | 0.496 ± 0.233 | 0.423 ± 0.149 |
| Reference | Number of samples per participant | 0.771 ± 0.160 | 0.712 ± 0.184 |

- **Adding the microbiome does not raise the Random Forest's AUROC:** 0.600 ± 0.114 with the two
  clinical variables alone, 0.560 ± 0.268 with ANCOM-BC2 taxa and 0.496 ± 0.233 with the full
  microbiome. With either input the AUROC is higher in 2 folds and lower in 3.
- **With the two clinical variables alone, the elastic net ranks preterm participants below term
  participants on average** (0.378 ± 0.240; in fold 5 both preterm participants score below all five
  term participants). Adding the ANCOM-BC2 taxa lowers it to 0.262 ± 0.086. Adding the full
  microbiome raises it to 0.471 ± 0.251 (higher in 3 folds, lower in 2): this is the only one of the
  12 model–set pairs in which the microbiome raises the mean AUROC, and the result is still below 0.5.
- The highest Approach 1 value, 0.600 ± 0.114, is within one fold-to-fold SD of 0.5 and below the
  number of samples per participant on the same folds (0.771 ± 0.160,
  [`SAMPLE_COUNT_REFERENCE.md`](SAMPLE_COUNT_REFERENCE.md)). **In this cohort, no Approach 1 model
  discriminates clearly better than chance, with or without microbiome, and adding the microbiome
  features does not improve on the two clinical variables.**
- Gestational age at the visit is tied to how the samples were collected: the samples of preterm
  participants cover a shorter part of pregnancy ([`SAMPLE_COUNT_REFERENCE.md`](SAMPLE_COUNT_REFERENCE.md)).
  Whether either arm depends on that sampling pattern is not measured here.

AUROC by outer fold:

| Fold | Subjects (preterm) | glmnet: clinical only | glmnet: + ANCOM | glmnet: + full | RF: clinical only | RF: + ANCOM | RF: + full | Number of samples |
|---|---|---|---|---|---|---|---|---|
| 1 | 9 (3) | 0.611 | 0.389 | 0.556 | 0.611 | 0.611 | 0.389 | 0.667 |
| 2 | 9 (3) | 0.556 | 0.278 | 0.111 | 0.500 | 0.111 | 0.556 | 0.889 |
| 3 | 9 (3) | 0.389 | 0.278 | 0.500 | 0.611 | 0.722 | 0.833 | 0.556 |
| 4 | 9 (3) | 0.333 | 0.167 | 0.389 | 0.778 | 0.556 | 0.500 | 0.944 |
| 5 | 7 (2) | 0.000 | 0.200 | 0.800 | 0.500 | 0.800 | 0.200 | 0.800 |

## 2. The three clinical sets with and without microbiome

AUROC:

| Clinical variable set | Model | Clinical only | + ANCOM-BC2 taxa | + full microbiome |
|---|---|---|---|---|
| 1 — DREAM (minimal) | Elastic net (glmnet) | 0.378 ± 0.240 | 0.262 ± 0.086 | 0.471 ± 0.251 |
| 1 — DREAM (minimal) | Random Forest | 0.600 ± 0.114 | 0.560 ± 0.268 | 0.496 ± 0.233 |
| 2 — literature | Elastic net (glmnet) | 0.507 ± 0.234 | 0.344 ± 0.250 | 0.351 ± 0.221 |
| 2 — literature | Random Forest | 0.636 ± 0.156 | 0.558 ± 0.190 | 0.562 ± 0.236 |
| 3 — data-driven | Elastic net (glmnet) | **0.802 ± 0.233** | 0.738 ± 0.270 | 0.669 ± 0.184 |
| 3 — data-driven | Random Forest | 0.711 ± 0.144 | 0.580 ± 0.200 | 0.531 ± 0.284 |
| Reference | Number of samples per participant | 0.771 ± 0.160 | | |

PR-AUC:

| Clinical variable set | Model | Clinical only | + ANCOM-BC2 taxa | + full microbiome |
|---|---|---|---|---|
| 1 — DREAM (minimal) | Elastic net (glmnet) | 0.410 ± 0.191 | 0.225 ± 0.036 | 0.453 ± 0.204 |
| 1 — DREAM (minimal) | Random Forest | 0.517 ± 0.137 | 0.480 ± 0.181 | 0.423 ± 0.149 |
| 2 — literature | Elastic net (glmnet) | 0.403 ± 0.213 | 0.350 ± 0.201 | 0.350 ± 0.207 |
| 2 — literature | Random Forest | 0.480 ± 0.227 | 0.395 ± 0.242 | 0.404 ± 0.139 |
| 3 — data-driven | Elastic net (glmnet) | **0.751 ± 0.295** | 0.716 ± 0.297 | 0.627 ± 0.235 |
| 3 — data-driven | Random Forest | 0.604 ± 0.222 | 0.579 ± 0.244 | 0.481 ± 0.222 |
| Reference | Number of samples per participant | 0.712 ± 0.184 | | |

- Mean AUROC of the six combinations: **0.606 with clinical variables only**, 0.507 with ANCOM-BC2
  taxa and 0.513 with the full microbiome. Five, four and three of the six are above 0.5.
- **In 11 of the 12 model–set pairs, adding microbiome features lowers the mean AUROC** (by 0.040 to
  0.180) **and the mean PR-AUC** (by 0.025 to 0.185). The exception is the elastic net with the
  Approach 1 variables and the full microbiome (section 1).
- **The highest of the 18 combinations is the elastic net with the data-driven clinical set and no
  microbiome features: 0.802 ± 0.233** (PR-AUC 0.751 ± 0.295). In every fold it uses exactly the
  same clinical variables as the highest combination with microbiome, the same model and set with
  ANCOM-BC2 taxa (0.738 ± 0.270). It is the only one of the 18 above the number of samples per
  participant (0.771 ± 0.160): fold by fold it is higher in folds 1 and 3, equal in fold 2 and lower
  in folds 4 and 5.
- In Approach 3, the variables selected in four of the five folds are `workouthome`,
  `low_education`, `extreme_bmi`, `extreme_age` and `comp1tribleed`; gestational age at the visit is
  selected in one fold. The selection in each fold is listed in
  [`results/socioeconomic_indicators/clinical_variables_by_fold.csv`](../results/socioeconomic_indicators/clinical_variables_by_fold.csv)
  (rows with `run` = `after`; it is the same with and without microbiome).
- Sensitivity, specificity and balanced accuracy at the threshold chosen on the inner validation
  split are in [`results/clinical_only_arm/metrics_by_combination.csv`](../results/clinical_only_arm/metrics_by_combination.csv).

## 3. Paired differences by fold

With microbiome minus clinical only, on the same folds:

| Clinical variable set | Model | Microbiome input | ΔAUROC, mean (min to max) | Folds higher / lower / equal | ΔPR-AUC, mean |
|---|---|---|---|---|---|
| 1 — DREAM (minimal) | Elastic net (glmnet) | + ANCOM-BC2 taxa | −0.116 (−0.278 to +0.200) | 1 / 4 / 0 | −0.185 |
| 1 — DREAM (minimal) | Elastic net (glmnet) | + full microbiome | +0.093 (−0.444 to +0.800) | 3 / 2 / 0 | +0.043 |
| 1 — DREAM (minimal) | Random Forest | + ANCOM-BC2 taxa | −0.040 (−0.389 to +0.300) | 2 / 3 / 0 | −0.037 |
| 1 — DREAM (minimal) | Random Forest | + full microbiome | −0.104 (−0.300 to +0.222) | 2 / 3 / 0 | −0.093 |
| 2 — literature | Elastic net (glmnet) | + ANCOM-BC2 taxa | −0.162 (−0.389 to +0.000) | 0 / 4 / 1 | −0.053 |
| 2 — literature | Elastic net (glmnet) | + full microbiome | −0.156 (−0.500 to +0.000) | 0 / 3 / 2 | −0.053 |
| 2 — literature | Random Forest | + ANCOM-BC2 taxa | −0.078 (−0.333 to +0.056) | 1 / 3 / 1 | −0.085 |
| 2 — literature | Random Forest | + full microbiome | −0.073 (−0.200 to +0.056) | 1 / 4 / 0 | −0.076 |
| 3 — data-driven | Elastic net (glmnet) | + ANCOM-BC2 taxa | −0.064 (−0.222 to +0.056) | 1 / 3 / 1 | −0.035 |
| 3 — data-driven | Elastic net (glmnet) | + full microbiome | −0.133 (−0.278 to +0.000) | 0 / 4 / 1 | −0.124 |
| 3 — data-driven | Random Forest | + ANCOM-BC2 taxa | −0.131 (−0.222 to +0.000) | 0 / 4 / 1 | −0.025 |
| 3 — data-driven | Random Forest | + full microbiome | −0.180 (−0.400 to +0.000) | 1 / 4 / 0 | −0.123 |

Over the 60 fold pairs (12 model–set–input pairs × 5 folds), the AUROC with microbiome is lower in
41, higher in 12 and equal in 7.

## 4. How to read these results

- **They are descriptive.** No test of the difference between the arms was run. Each test fold
  holds 7–9 subjects, the fold-to-fold SDs of the AUROC range from 0.09 to 0.28, and the 12 pairs are
  not independent: they share the subjects, the folds and the clinical variables. That 11 of 12
  pairs go the same way is therefore not 11 independent replications.
- **They concern these models and inputs.** They do not show that the vaginal microbiome carries no
  information about preterm birth. They show that, with these two models and their fixed
  hyperparameters, adding these microbiome inputs (5 to 19 ANCOM-BC2 taxa per fold, or the 49
  filtered genera, each with Shannon diversity) to the clinical variables known at the visit did not
  raise the cross-validated discrimination in this cohort. Why the microbiome features lower it was
  not examined; each inner-training split holds 23 subjects.
- **Only clinical variables known at the visit** enter both arms. The sensitivity analysis that adds
  the variables known only after the visit was not repeated for the clinical-only arm.

**Not covered here:** a permutation test of any clinical-only combination; the same comparison
restricted to the samples taken before a fixed gestational week (landmark design).

> 📌 **Update (2026-10-23).** The permutation test has since been run across the 18 combinations,
> adjusted for the choice among them: selection-adjusted p = 0.060 for the clinical-only elastic net
> with Approach 3 (the best of the 18) and 0.190 for the main model; see
> [`results/permutation_test/`](../results/permutation_test/).

> 📌 **Update (2026-11-02).** How large a difference with and without microbiome this cohort could
> detect is now computed in [`AUROC_PRECISION.md`](AUROC_PRECISION.md).
>
> - **From a model at 0.70:** with the pipeline's metric, about 0.28 or more for 80% power.
> - **From 0.80:** no attainable difference reaches that power.
> - **Reading:** the results of this document are compatible both with no contribution of the
>   microbiome features and with a moderate one.

## 5. How the clinical-only arm was run

The engine, the analysis notebook and the configuration are unchanged.
[`scripts/clinical_only_arm.R`](../scripts/clinical_only_arm.R) runs the notebook's chunks by name
up to `run_nested_cv`, which runs the 12 combinations with microbiome as the pipeline does. It then
calls the engine, `train_with_nested_cv()`, for each model and clinical set defined in that chunk,
with a microbiome table that holds only the sample and subject keys. The engine treats as taxa every
column that is not a key, the outcome or a clinical variable, so it finds none and applies no CLR
transformation. The microbiome input is passed as `"Full_Microbiome"`, so ANCOM-BC2 is not run, and
the count table and genus list it would need are not passed, so the run would stop if it tried.

The script stops unless all of the following hold:

- **The 12 combinations with microbiome reproduce the current main analysis:** identical per-fold
  metrics, subject-level test predictions, ANCOM-BC2 taxa and clinical variables to the saved run of
  [`SOCIOECONOMIC_INDICATORS.md`](SOCIOECONOMIC_INDICATORS.md) (`sesfix_2026-10-16`, when its output
  is present), and the same values as the `after` rows of
  [`results/socioeconomic_indicators/metrics_by_combination.csv`](../results/socioeconomic_indicators/metrics_by_combination.csv).
- **No microbiome feature reaches the clinical-only arm:** the engine reports 0 taxa columns in each
  of the 5 folds and in the model refitted on all the data; the predictors of that model are exactly
  the clinical variables of the set, with no taxon, no Shannon diversity and no CLR column; and no
  object with the names the engine looks up for its per-fold microbiome table exists before any call.
- **The only difference is the microbiome:** each clinical-only combination has results in the 5
  folds, and uses the same clinical variables (the fixed list of Approaches 1 and 2, the same
  selection in each fold of Approach 3) and the same test subjects in each fold as its two
  counterparts with microbiome.
- **The reference reproduces [`SAMPLE_COUNT_REFERENCE.md`](SAMPLE_COUNT_REFERENCE.md):** the
  per-fold AUROC of the number of samples per participant (score = minus the number of samples,
  orientation fixed in advance) equals the published one. Its PR-AUC is reported here for the first
  time.

## 6. Files and reproducing these numbers

[`results/clinical_only_arm/`](../results/clinical_only_arm/) holds the numbers of this document:

- `metrics_by_combination.csv` — the 18 combinations (`Microbiome` = `None`, `ANCOM_Taxa`,
  `Full_Microbiome`): clinical variables, number of folds, and mean and SD of AUROC, PR-AUC,
  sensitivity, specificity and balanced accuracy;
- `metrics_by_fold.csv` — the same metrics by outer fold;
- `microbiome_vs_clinical_only.csv` — the paired differences of section 3;
- `sample_count_reference_by_fold.csv` — AUROC and PR-AUC of the number of samples per participant
  by outer fold.

The script needs the cohort data in `data/raw/`:

```bash
PTB_RUN_TAG=clinonly_2026-10-16 PTB_RESULTS_DIR=results/clinical_only_arm Rscript scripts/clinical_only_arm.R
```

It took 3.8 and 4.8 minutes in two runs with identical results (Windows 11, R 4.4.2), and prints the
checks of section 5 and the tables of this document. The comparison with the saved run uses the output of the "after" run of
[`SOCIOECONOMIC_INDICATORS.md`](SOCIOECONOMIC_INDICATORS.md) §6; without it, the script says that
this check was skipped and still compares with the versioned CSV.
