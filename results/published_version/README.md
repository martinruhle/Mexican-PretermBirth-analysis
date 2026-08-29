# `results/published_version/` — output of the published analysis

> **📌 These are archived outputs. They correspond to the analysis as published and predate the
> methodological updates described in
> [`docs/UPDATE_SINCE_PUBLICATION.md`](../../docs/UPDATE_SINCE_PUBLICATION.md). Please cite the
> current results rather than the numbers in this directory.**

This directory is kept so that the state of the analysis at publication time remains inspectable
alongside the updated state. It is a record, not a current result.

## What is here

| Path | What it is |
|---|---|
| `integrated_preterm_prediction_workflow.html` | Rendered report of the pipeline, from a run made before these updates. |
| `figures/` | The 11 figures produced by that same run. |

## Why these numbers differ from the current ones

This render predates even the published article: it was produced **before the contaminant filtering
step was applied**, so it works from 59 genera instead of 49 and reports a best-model AUROC of
**0.849**. The published article reports **0.813** (post-filtering). Both predate the updates below.

The figure filenames also reflect that age — some are named after pipeline chunks
(`feature_importance_no_retrain`) that no longer exist in
`analysis/integrated_preterm_prediction_workflow.Rmd`.

Since this render was made, three methodological updates changed the results:

1. The CLR transformation now applies as intended (a pseudocount appropriate for counts had been
   added to relative abundances). This reorders the model ranking.
2. One genus (`Escherichia-Shigella`) is now retained in differential-abundance analysis.
3. The permutation test supporting the significance claim was rebuilt.

The current best model is **elastic net · Approach 3 · ANCOM-BC2 taxa, AUROC 0.760 ± 0.270**,
with a permutation p-value of **0.019**. Full detail, including what each update changed, is in
[`docs/UPDATE_SINCE_PUBLICATION.md`](../../docs/UPDATE_SINCE_PUBLICATION.md); the current results
table is in the [repository README](../../README.md#results).

## Where the current outputs live

- **Results table and headline numbers:** [`README.md`](../../README.md#results)
- **Permutation test artefacts (versioned):** [`results/permutation_test/`](../permutation_test/)
- **Regenerating the full report:** running
  `analysis/integrated_preterm_prediction_workflow.Rmd` writes a fresh HTML report to
  `analysis/_output/` (untracked). It requires the restricted cohort data; see
  [`docs/DATA_ACCESS.md`](../../docs/DATA_ACCESS.md).
