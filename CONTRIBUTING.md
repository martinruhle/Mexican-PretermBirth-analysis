# Contributing

Thanks for looking at this code. This is a research compendium for a published analysis, so the
priority is that results stay verifiable: a change that improves the code but moves a reported number
without saying so is worse than no change at all.

Please read [`docs/UPDATE_SINCE_PUBLICATION.md`](docs/UPDATE_SINCE_PUBLICATION.md) first. It
explains which published numbers were updated and why — most of the rules below exist because of one
of those updates.

---

## Setting up

```bash
git clone https://github.com/martinruhle/Mexican-PretermBirth-analysis.git
```

```r
install.packages("renv")
renv::restore()      # pinned R 4.4.2 library, including Bioconductor
```

`.Rprofile` activates `renv` for you. Load the engine for interactive work with:

```r
devtools::load_all()
```

You do **not** need the restricted cohort data to develop against this repository. The synthetic
`data/example/` dataset drives the test suite and can run the whole pipeline:

```bash
PTB_PROFILE=example Rscript scripts/run_baseline.R
```

---

## Running the tests

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

That skips the slow end-to-end test and finishes in a couple of minutes. Before opening a pull
request, run the full suite the way CI does:

```bash
PTB_RUN_SLOW_TESTS=1 Rscript -e 'testthat::test_dir("tests/testthat")'
```

Expect **56 passes, 0 failures, 0 skips**. Budget around 20 minutes: `test-leakage-permutation.R`
re-runs a reduced nested CV under 10 label permutations and accounted for ~17 of the 18.4 minutes in
a recent local run. (The estimate in that file's header comment is out of date and much too
optimistic.)

If anything skips with `PTB_RUN_SLOW_TESTS=1`, check your environment — a green
suite that skipped the leakage test is not green.

CI ([`.github/workflows/check-code.yml`](.github/workflows/check-code.yml)) restores the pinned
library, installs the package, runs the full suite with `PTB_RUN_SLOW_TESTS=1`, and smoke-runs the
example pipeline. All of it must pass.

Any new function in `R/` needs a `testthat` test.

---

## The invariants — please preserve these

These are the properties that keep the reported performance estimates trustworthy. Every one of them
is enforced by a test, and several were made explicit as part of the methodological updates
documented in `docs/UPDATE_SINCE_PUBLICATION.md`. A change that trips one of these is not a
refactor.

1. **Nothing is ever learned from the outer-test fold.** Within a fold: ANCOM-BC2 taxa selection
   and the Approach-3 univariate screening are fitted on the **outer-training** subjects; the
   preprocessing recipe and the model are fitted on **inner-train** and then *applied* to
   inner-validation and outer-test. Never the other way around, and never on the pooled data. If
   you change which subset a step is fitted on, say so explicitly — it changes what the reported
   numbers mean.
2. **The classification threshold is chosen on inner-validation and evaluated on outer-test.**
   There is no inner *k*-fold — the inner loop is a single stratified 70/30 split
   (`cv.inner_train_prop`). Do not "simplify" it into a k-fold; the numbers would change.
3. **Splits are at the subject level** (the `id` column), stratified by `preterm`, seed 123. A
   subject's longitudinal samples must never straddle a train/test boundary.
4. **Microbiome and clinical columns are separated by the explicit list in
   `config/data_dictionary.csv`** (`role == "microbiome"`). Never by positional indices
   (`data[, 1:99]`), never by a name regex. The positional version had overlapping ranges that put
   the key columns in both domains, and it is what prevented the pipeline being used on other
   datasets.
5. **Zero replacement is outcome-blind.** `fit_clr_zerorepl()` learns replacement levels once, over
   the full microbiome, without ever seeing the outcome, and they are passed in as the `zero_levels`
   argument. Do not make them depend on labels or on a fold.
6. **The CLR is per sample.** Each sample is centred by the geometric mean of its own taxa, so each
   transformed row sums to zero. Do not centre per feature, and do not add a pseudocount to
   relative abundances — the earlier version did, which prevented the compositional transform from
   taking effect (see update 2.1).
7. **Each figure and table comes from the same model, trained once.** Do not refit a model to draw
   a plot. Every builder in `R/plots.R` takes already-computed results and returns a ggplot.
8. **Permuting labels must reach every consumer of the outcome.** The outcome is read in three
   places (the model's outcome column, the stratification labels, and the metadata that ANCOM-BC2
   uses). Permuting fewer than all three contaminates the null — ANCOM keeps selecting taxa with
   the true labels, and the null comes out far too high.

---

## `R/` is the only implementation of the engine

**Do not copy engine code into a script, notebook or `.Rmd` in order to modify it.** If a script
needs different behaviour, add a parameter to the function in `R/` and call it.

This is not a style preference; it comes from experience on this project. The earlier permutation
test was a second, copied implementation of the nested CV loop. Because it was a copy, it drifted out
of sync with the engine as the engine evolved: it ended up pointing at a superseded reference value
and a different model, and using an AUROC orientation that placed the null distribution near 0.68
rather than 0.5. Keeping one implementation is what prevents that class of divergence. It has since
been replaced by [`scripts/permutation_test.R`](scripts/permutation_test.R), which re-executes the
real pipeline setup and calls the engine directly. Full account:
[`docs/UPDATE_SINCE_PUBLICATION.md` §2.4](docs/UPDATE_SINCE_PUBLICATION.md).

Related rules:

- **Values belong in `config/config.yml`, not in code.** Paths, column names, approach definitions,
  preprocessing parameters, CV settings and model hyperparameters are all configuration. Adapting
  the pipeline to a new dataset must never require editing `R/`.
- **Functions in `R/` are pure.** Everything enters through arguments; no reading of globals, no
  hidden side effects, no `setwd()`, no absolute paths — resolve from the project root with
  `here::here()`. One responsibility per function, roxygen on each.
- **Match the surrounding code.** Comment density, naming and idiom vary between files; follow the
  file you are editing.

---

## Changing the engine means re-checking the numbers

Refactoring must not change results. The workflow is:

1. Capture a baseline **before** your change:
   ```bash
   PTB_RUN_TAG=before Rscript scripts/run_baseline.R
   ```
2. Make the change.
3. Capture again with `PTB_RUN_TAG=after` and diff
   `analysis/_output/before_metrics.csv` against `after_metrics.csv` row by row.

Any difference is something to investigate and explain, not noise. If a change is *meant* to move
the numbers, say so explicitly in the pull request, quantify the shift, and update
[`docs/UPDATE_SINCE_PUBLICATION.md`](docs/UPDATE_SINCE_PUBLICATION.md) and the results table in the
[README](README.md) in the same change. Those two documents are what readers arriving from the
article rely on, so keeping them in step with the code is part of the change, not a follow-up.

Without the restricted data you cannot reproduce the real-cohort baseline. Use
`PTB_PROFILE=example` for a same-vs-same comparison instead: the absolute values are meaningless,
but a refactor that is result-neutral on the real data should be result-neutral on the example data
too.

---

## Proposing a change

1. Branch off `master`.
2. Keep the change focused — one concern per pull request.
3. Add or update tests.
4. Run `PTB_RUN_SLOW_TESTS=1 Rscript -e 'testthat::test_dir("tests/testthat")'`.
5. If you touched `R/` roxygen, regenerate `NAMESPACE`:
   ```r
   devtools::document()
   ```
6. In the pull request, state whether the change is result-neutral and how you verified it.

---

## Known limitations and open items

Already tracked, so you are not surprised and do not need to report them. Contributions welcome as
separate pull requests.

- **`scripts/sensitivity_nonindependence_weight.R` predates the current engine API:** it calls
  `apply_clr_transform(data, taxa_cols)` without the now-required `zero_levels` argument, so it needs
  migrating onto the current `R/` interface before it will run.
- **`R CMD check` emits NOTES.** The engine calls several analysis packages unqualified (a
  legacy of its `.Rmd` origins) and `PRROC` and `microbiome` are used at runtime without being
  declared in `DESCRIPTION`. It installs and tests cleanly; full check-cleanliness is still open.
- **`renv.lock` contains two hand-pinned packages** (`microbiome`, `Rtsne`). `microbiome` is a
  `Suggests` of ANCOMBC that this pipeline genuinely needs at runtime, and renv's implicit snapshot
  does not follow `Suggests`. **Do not run `renv::snapshot()` without re-pinning them** — a
  snapshot can drop them, and a clean `renv::restore()` would then leave the ANCOM tests without a
  dependency they need.
- **The pipeline prints `12 arguments not used by format` warnings** — a `sprintf` call in the
  completeness helper inside the analysis `.Rmd` whose format string needs consolidating. It affects
  the log only, not results.
- **AUROC uses `pROC`'s automatic orientation** in the engine (`R/nested_cv.R`, `R/threshold.R`).
  Immaterial with real signal, but it raises the AUROC of *random* labels above 0.5, which is why
  the end-to-end permutation test asserts strictly on balanced accuracy and only loosely on AUROC.
  Fixing the orientation would let that guard be tightened — but it is a results-affecting change,
  so it needs the baseline comparison above.
- **Alphanumeric subject/sample ids are not yet supported by the PCA figure**, and the exploratory
  zero replacement needs each taxon present in at least 2 samples. Both are documented under
  [Known reusability limits](README.md#known-reusability-limits).

---

## Reporting a problem

Open an issue with: what you ran (command and `PTB_PROFILE`), what you expected, what happened
(full console output), and the output of `sessionInfo()`. If it concerns results, say whether you were
using the real cohort data or `data/example/`.

Please be constructive, keep discussion on scientific and reproducibility merit, and remember that
this work involves sensitive health data from a small cohort.
