# ============================================================================
# Test setup (Chat 7) — attach the packages the engine calls UNQUALIFIED, then
# load R/ into scope. renv is not restored yet (Chat 8), so this runs against
# the system library where tidymodels / ANCOMBC / phyloseq / pROC / PRROC /
# zCompositions / vegan already live.
# ============================================================================

suppressWarnings(suppressMessages(suppressPackageStartupMessages({
  # Heavy / bioconductor deps first (guarded so the pure tests still run if a
  # slow test's package is absent); tidymodels LAST so dplyr's verbs win the
  # common conflicts (filter/select/slice) exactly as the engine expects.
  if (requireNamespace("phyloseq", quietly = TRUE)) library(phyloseq)
  if (requireNamespace("ANCOMBC",  quietly = TRUE)) library(ANCOMBC)
  library(vegan)
  library(zCompositions)
  library(pROC)
  if (requireNamespace("PRROC", quietly = TRUE)) library(PRROC)
  library(tidymodels)
})))

# Load the engine. The proper path is devtools::load_all(), but the package
# DESCRIPTION lists `config` under Imports and it is not installed (renv is
# unrestored), so load_all() aborts on the dependency check. Fall back to
# sourcing R/*.R into the global environment — functionally identical for the
# tests, and it also lets train_with_nested_cv() resolve its `subject_labels`
# global-dep via the search path. Once renv + NAMESPACE land in Chat 8,
# load_all() (and tests/testthat.R's test_check) become the real path.
local({
  loaded <- FALSE
  if (requireNamespace("devtools", quietly = TRUE)) {
    loaded <- tryCatch({
      suppressWarnings(suppressMessages(devtools::load_all(here::here(), quiet = TRUE)))
      TRUE
    }, error = function(e) FALSE)
  }
  if (!loaded) {
    for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) {
      source(f, local = FALSE, encoding = "UTF-8")
    }
  }
})
