# Results with the education indicator flagging secondary education

Numbers behind [`docs/SOCIOECONOMIC_INDICATORS.md`](../../docs/SOCIOECONOMIC_INDICATORS.md): the 12
combinations of the main analysis (clinical variables known at the sample visit) run twice,
differing only in how `low_education`, and therefore `ses_risk_score`, is derived. Column `run`:

| `run` | `low_education` = 1 for |
|---|---|
| `before` | primary education only (the education level was compared as text) |
| `after` | primary and secondary education, as defined in [`config/data_dictionary.csv`](../../config/data_dictionary.csv) |

- `metrics_by_combination.csv` — one row per run and combination (`Model`, `Approach`,
  `Microbiome`): mean and SD across the 5 outer folds of AUROC, PR-AUC, sensitivity, specificity and
  balanced accuracy. `uses_ses_variables` names `low_education` or `ses_risk_score` when the
  combination used them (in its fixed list, or in at least one fold of Approach 3).
- `clinical_variables_by_fold.csv` — the clinical variables each combination used: the fixed list
  of Approaches 1 and 2 (`fold` = `all`) and the per-fold selection of Approach 3, with the
  `availability` class of each variable from the data dictionary.

Both files are written by [`scripts/compare_ses_runs.R`](../../scripts/compare_ses_runs.R) from the
two runs; the commands are in the last section of the document.
