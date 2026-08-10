# `results/published_version/` — output of the published analysis

> **⚠ These are archived outputs. They correspond to the analysis as published and predate the
> methodological corrections described in
> [`docs/UPDATE_SINCE_PUBLICATION.md`](../../docs/UPDATE_SINCE_PUBLICATION.md). Do not cite the
> numbers in this directory as current results.**

This directory is kept so that the state of the analysis at publication time remains inspectable
alongside the corrected state. It is a record, not a result.

## What is here

| Path | What it is |
|---|---|
| `integrated_preterm_prediction_workflow.html` | Rendered report of the pipeline, from a run made before the corrections. |
| `figures/` | The 11 figures produced by that same run. |

## Why these numbers differ from the current ones

This render is even older than the published article: it was produced **before the contaminant
filtering step was applied**, so it works from 59 genera instead of 49 and reports a best-model
AUROC of **0.849**. The published article reports **0.813** (post-filtering). Both predate the
corrections below.

The figure filenames also reflect that age — some are named after pipeline chunks
(`feature_importance_no_retrain`) that no longer exist in
`analysis/integrated_preterm_prediction_workflow.Rmd`.

Since this render was made, three corrections changed the results:

1. The microbiome was never actually CLR-transformed (a pseudocount appropriate for counts was
   added to relative abundances). Fixing it reorders the model ranking.
2. One genus (`Escherichia-Shigella`) was silently dropped from differential-abundance analysis.
3. The permutation test supporting the significance claim was invalid and did not run.

The current best model is **elastic net · Approach 3 · ANCOM-BC2 taxa, AUROC 0.760 ± 0.270**,
with a permutation p-value of **0.019**. Full detail, including what each correction changed, is
in [`docs/UPDATE_SINCE_PUBLICATION.md`](../../docs/UPDATE_SINCE_PUBLICATION.md); the corrected
results table is in the [repository README](../../README.md#results-current).

## Where the current outputs live

- **Results table and headline numbers:** [`README.md`](../../README.md#results-current)
- **Permutation test artefacts (versioned):** [`results/permutation_test/`](../permutation_test/)
- **Regenerating the full report:** running
  `analysis/integrated_preterm_prediction_workflow.Rmd` writes a fresh HTML report to
  `analysis/_output/` (untracked). It requires the restricted cohort data; see
  [`docs/DATA_ACCESS.md`](../../docs/DATA_ACCESS.md).
