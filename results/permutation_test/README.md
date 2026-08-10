# Permutation test — saved artefacts

These files are the **saved output of the permutation test** for the model reported in
[`docs/UPDATE_SINCE_PUBLICATION.md`](../../docs/UPDATE_SINCE_PUBLICATION.md):

> **Elastic net (glmnet) · Approach 3 (data-driven) · ANCOM-BC2-selected taxa**

They are versioned deliberately. Regenerating the null distribution takes **≈3.6 hours on 8
cores**, so having the result in the repository lets a reader check the p-value, the shape of
the null and the per-fold numbers without re-running anything — and without installing R, in
the case of the figure.

## Result

| | Value |
|---|---|
| Observed AUROC (5 outer folds, real cohort) | **0.760** |
| Null distribution, 999 subject-level label permutations (998 valid) | **0.502 ± 0.121** (median 0.504, range 0.184–0.913) |
| Permutations reaching the observed value | 18 / 998 |
| **One-sided permutation p-value** | **0.019** |

The null is centred on 0.502 and 48.6% of its mass falls below 0.5, which is the check that
the null is legitimate (a permuted-label AUROC should average 0.5). With the statistic used
previously the same data gave p = 0.139; see
[`docs/UPDATE_SINCE_PUBLICATION.md` §2.4](../../docs/UPDATE_SINCE_PUBLICATION.md) for what changed
and why.

![Null distribution of AUROC under 999 subject-level label permutations](permutation_null_final999.png)

## Files

| File | What it is |
|---|---|
| `permutation_null_final999.png` | Histogram of the null distribution with the observed AUROC and the p-value annotated. Raster version, renders inline on GitHub. |
| `permutation_null_final999.pdf` | The same figure, vector format, for publication. |
| `permutation_summary_final999.csv` | One-row summary: target model/approach/microbiome input, observed AUROC, null mean/SD/median/min/max, `n_permutations` (999), `n_valid` (998), `n_exceeding` (18) and `p_value` (0.019). The final column carries 0.139, the value the automatic-orientation statistic used previously would have reported, kept for comparison. |
| `permutation_null_final999.csv` | The null distribution as plain text: 999 rows with `permutation`, `seed`, `auroc_fixed` (fixed-orientation AUROC, the statistic used), `auroc_auto` (automatic-orientation AUROC, kept for comparison) and `n_folds_valid`. |
| `perm_null_final999.rds` | The raw R object behind the two CSVs: the 999 null values (fixed and automatic orientation), valid folds per permutation, the seeds, the wall-clock runtime, and the observed run embedded for reference. |
| `perm_observed.rds` | The unpermuted run: observed AUROC (fixed and automatic orientation), per-fold metrics, per-fold fixed-orientation AUROC, and the out-of-fold subject-level predictions. |
| `supplementary_per_fold_metrics.csv` | Per-fold table for the reported model (one row per outer fold): `n_test`, `n_ptb`, `threshold`, `AUROC`, `AUROC_fixed`, `PRAUC`, `Sensitivity`, `Specificity`, `Balanced_Accuracy`, `Youden`. Prepared as a supplementary table in response to reviewer comments. |

## How to regenerate

All three modes are driven by [`scripts/permutation_test.R`](../../scripts/permutation_test.R)
and require the **restricted cohort data** in `data/raw/` (see
[`docs/DATA_ACCESS.md`](../../docs/DATA_ACCESS.md)). Outputs are written to
`analysis/_output/` (untracked scratch); the files here are copies of that output.

```bash
# 1. Observed value + per-fold supplementary table (~40 s)
#    -> perm_observed.rds, supplementary_per_fold_metrics.csv
PTB_PERM_MODE=verify Rscript scripts/permutation_test.R
```

```bash
# 2. Null distribution: 999 permutations, ~3.6 h on 8 cores
#    -> perm_null_final999.rds
PTB_PERM_MODE=run PTB_PERM_N=999 PTB_PERM_CORES=8 PTB_PERM_TAG=final999 Rscript scripts/permutation_test.R
```

```bash
# 3. p-value, summary tables and figure from the saved null (seconds)
#    -> permutation_null_final999.{csv,png,pdf}, permutation_summary_final999.csv
PTB_PERM_MODE=report PTB_PERM_TAG=final999 Rscript scripts/permutation_test.R
```

Step 3 only reads the `.rds` files, so the figure and the summary table can be rebuilt from
`perm_null_final999.rds` and `perm_observed.rds` in this directory without the cohort data —
copy them into `analysis/_output/` first.

The permutation seeds are fixed (`10001`–`10999`), so the null distribution and the p-value
are reproducible; this was confirmed across independent runs.
