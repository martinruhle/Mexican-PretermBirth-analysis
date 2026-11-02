# Results by clinical variable availability

Numbers behind [`docs/CLINICAL_AVAILABILITY_RESULTS.md`](../../docs/CLINICAL_AVAILABILITY_RESULTS.md):
the 12 combinations run three times, differing only in the availability classes admitted as clinical
predictors. Column `analysis`:

| `analysis` | Classes admitted |
|---|---|
| `before` | `at_visit`, `after_visit`, `outcome_defined` (reproduces the run made before the pipeline read the classification) |
| `main` | `at_visit` — the main analysis |
| `sensitivity` | `at_visit`, `after_visit` |

- `metrics_by_combination.csv` — one row per run and combination (`Model`, `Approach`,
  `Microbiome`): mean and SD across the 5 outer folds of AUROC, PR-AUC, sensitivity, specificity and
  balanced accuracy.
- `clinical_variables_by_fold.csv` — the clinical variables each combination used: the fixed list
  of Approaches 1 and 2 (`fold` = `all`) and the per-fold selection of Approach 3, with the
  `availability` class of each variable from [`config/data_dictionary.csv`](../../config/data_dictionary.csv).

Both files are written by [`scripts/compare_availability_runs.R`](../../scripts/compare_availability_runs.R)
from the three runs; the commands are in section 7 of the document.
