# Clinical-only arm and the 12 combinations with microbiome

Numbers behind [`docs/CLINICAL_ONLY_ARM.md`](../../docs/CLINICAL_ONLY_ARM.md): the 2 models × 3
clinical variable sets of the main analysis (clinical variables known at the sample visit), run
with the same engine, outer folds and seeds with three microbiome inputs. Column `Microbiome`:

| `Microbiome` | Microbiome features in the model |
|---|---|
| `None` | none: clinical variables only |
| `ANCOM_Taxa` | the taxa selected by ANCOM-BC2 in each fold, plus Shannon diversity |
| `Full_Microbiome` | the 49 filtered genera, plus Shannon diversity |

- `metrics_by_combination.csv` — one row per combination (`Model`, `Approach`, `Microbiome`): the
  clinical variables (the fixed list of Approaches 1 and 2; Approach 3 selects them in each fold),
  the number of folds with results, and the mean and SD across the 5 outer folds of AUROC, PR-AUC,
  sensitivity, specificity and balanced accuracy.
- `metrics_by_fold.csv` — the same metrics for each outer fold, with the fold's number of test
  subjects (`n_test`), preterm subjects (`n_ptb`) and the classification threshold chosen on the
  inner validation split.
- `microbiome_vs_clinical_only.csv` — for each model and clinical set, the paired difference by
  fold between each microbiome input and the clinical-only arm (with microbiome minus clinical only):
  mean, SD, minimum and maximum of the AUROC difference, the number of folds in which it is higher,
  lower or equal, and the mean and SD of the PR-AUC difference.
- `sample_count_reference_by_fold.csv` — the reference of
  [`docs/SAMPLE_COUNT_REFERENCE.md`](../../docs/SAMPLE_COUNT_REFERENCE.md) (score = minus the number
  of samples per participant) on the same outer folds: AUROC and PR-AUC by fold.

The Approach 3 selection in each fold is the same in the three inputs; it is listed in
[`results/socioeconomic_indicators/clinical_variables_by_fold.csv`](../socioeconomic_indicators/clinical_variables_by_fold.csv)
(rows with `run` = `after`).

All four files are written by [`scripts/clinical_only_arm.R`](../../scripts/clinical_only_arm.R);
the command is in the last section of the document.
