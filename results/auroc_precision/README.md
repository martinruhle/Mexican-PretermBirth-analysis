# Precision of the AUROC and the smallest detectable difference

Numbers behind [`docs/AUROC_PRECISION.md`](../../docs/AUROC_PRECISION.md).

## Terms used in the column names

- **AUROC values:** `A` is the model's true AUROC; `n1` and `n0` are the numbers of preterm and term
  participants.
- **`pooled`** refers to a fixed model with a single AUROC over all participants.
- **`engine` or `foldmean`** refers to the pipeline's metric: the mean of the 5 per-fold AUROCs. The
  `engine` values also include the refitting factor f.
- **`HM`** is the Hanley–McNeil variance.
- **`delta`** is the smallest difference in AUROC detected with a two-sided test at α = 0.05 and 80%
  power. `NA` means that no difference with A + Δ ≤ 1 reaches that power.

## Files

### Reproduction of the preliminary estimate

- **`prior_reproduction.csv`** — the preliminary planning estimate recomputed with its own method
  (Δ = 2.8 × SE_HM, the same SE for both models). For each row it gives:
  - `SE`, `CI_halfwidth` (1.96 × SE), `delta_prior` and `A_plus_delta`;
  - the value as stated (`stated`) and whether it is reproduced to three decimals (`reproduced`).
- **`prior_cohort_proportion.csv`** — the same estimate for N = 81–500 with this cohort's preterm
  share (14/43), with `n1` exact and rounded.

### Simulation

- **`simulation_check.csv`** — the binormal simulation for a fixed model, at A = 0.60, 0.70 and 0.80
  with r = 0.5:
  - the score correlation used and the empirical r;
  - the false-positive rate of the paired DeLong test with no true difference;
  - simulated and Hanley–McNeil SEs, for the pooled AUROC and the fold mean;
  - Δ at 80% power from the simulated test (`delta_simulated_delong`);
  - Δ from the Obuchowski form (`delta_analytic`, `delta_oracle_binormal`) and from the Wald form
    (`delta_analytic_wald`, `delta_oracle_binormal_wald`). `oracle` uses the variances of the
    simulated binormal scores instead of Hanley–McNeil;
  - the engine's Δ with Hanley–McNeil and binormal variances;
  - the largest Δ on the simulation grid.

### Refitting and the spread between folds

- **`permutation_null_sd.csv`** — for each of the 18 combinations, from the permutation test:
  - `null_mean` and `null_sd`: mean and SD of the permuted AUROCs;
  - `null_sd_test_sampling`: the SD with random scores on the same folds, without refitting;
  - `inflation`: their ratio. Its median over the 18 is f, in `parameters.csv`.
- **`pairs_permutation_null.csv`** — for each of the 12 microbiome pairs, under the same
  permutations: the correlation of the two models' AUROCs (`r_null`) and the SD of their difference.
- **`between_fold_sd.csv`** — for each of the 18 combinations of the clinical-only arm: the observed
  between-fold SD and variance of the AUROC, and the variance expected from sampling the test folds
  alone.

### Correlation between the paired models

- **`pairs_delong.csv`** — for each of the 12 microbiome pairs, from the saved out-of-fold
  predictions:
  - the pooled AUROC of both models;
  - the DeLong correlation, over all participants (`r_pooled`) and combined over folds as the mean of
    the folds (`r_folds`);
  - the SE of the difference for both metrics.

### Counts by gestational-week window

- **`window_counts.csv`** — participants with at least one sample in each window.
  - **Columns:** `eligible`, `participants`, cases and non-cases for the window's outcome, preterm
    and term, eligible participants without a sample in the window, and samples in the window.
  - **Low counts:** for the outcome delivery before 32+0, `cases` is reported as `<= 1`. The DREAM
    sub-challenge 2 rows give no further counts (`note` = `not evaluable (<= 1 case)`).

### Tables by cohort size

- **`precision_by_size.csv`** — for each cohort row and A = 0.60, 0.70 and 0.80:
  - N, `n1`, `n0` and the preterm share;
  - SE and 95% CI half-width for a fixed model (`pooled`) and for the pipeline's metric (`engine`);
  - the preliminary 2.8 × SE.
- **`detectable_delta_by_size.csv`** — Δ for the same rows, in long format:
  - **`r_label`:** 0.3, 0.5, 0.7, or one of the three medians of `parameters.csv`;
  - **`critical`:** `null` is the Obuchowski form, critical value from the variance under no
    difference; `alternative` is the Wald form, used in the document;
  - **`delta_pooled` and `delta_engine`:** the fixed model and the pipeline's metric.

### Parameters

- **`parameters.csv`** — α, power, the two z values and the number of simulation replicates.
  - **The refitting factor:** f with its bootstrap 95% CI.
  - **The null SDs** from sampling the test folds alone.
  - **The expected share of permutations with a skipped fold.**
  - **The values of r used.**

## How these files are generated

All files are written by [`scripts/auroc_precision.R`](../../scripts/auroc_precision.R). The
command and its data requirements are in the last section of the document.
