# ============================================================================
# Tiny, in-memory fixtures for the pure-function tests (Chat 7).
# testthat auto-sources helper*.R before running the suite, so every make_*()
# below is available to all test files. Passing `dict` explicitly to the
# contract / io functions keeps the pure tests free of disk I/O.
# ============================================================================

# Taxa used by the contract + io fixtures. Deliberately includes a hyphenated
# name (`Escherichia-Shigella`) and family-rank prefixes (`f__...`) so the tests
# exercise name preservation, not just plain genus names.
.FIXTURE_TAXA <- c("Lactobacillus", "Gardnerella", "Escherichia-Shigella",
                   "f__Rhizobiaceae", "Prevotella")

# Minimal config list: the `columns` block the contract/io functions read, plus
# an active `models` set and a reserved `models_future` set (which
# build_model_specs() must ignore).
make_cfg <- function() {
  list(
    columns = list(
      subject_id       = "id",
      sample_id        = "index",
      outcome          = "preterm",
      outcome_positive = "1",
      gestational_age  = "sdg_parto"),
    models = list(
      rf_base     = list(engine = "ranger", mtry = 4, trees = 500, min_n = 10,
                         importance = "impurity", num_threads = 1),
      glmnet_base = list(engine = "glmnet", penalty = 0.01, mixture = 0.5)),
    models_future = list(
      xgb_base  = list(engine = "xgboost", trees = 100, tree_depth = 3),
      tree_base = list(engine = "rpart",   cost_complexity = 0.01, tree_depth = 5)))
}

# Data dictionary: keys + outcome + gestational-age + a couple of clinical vars
# + the microbiome taxa (role == "microbiome"), in a fixed order.
make_dict <- function() {
  data.frame(
    variable = c("id", "index", "preterm", "sdg_parto",
                 "edad_cronologicamujer", "imc_pregestacional",
                 .FIXTURE_TAXA),
    role = c("subject_id", "sample_id", "outcome", "gestational_age",
             "clinical", "clinical",
             rep("microbiome", length(.FIXTURE_TAXA))),
    stringsAsFactors = FALSE)
}

# Wide matrix: 8 samples from 4 subjects (2 samples each), subject-level 0/1
# outcome (2 preterm / 2 term), a couple of clinical columns, and the taxa as
# relative abundances in [0, 1]. Valid against every contract check; individual
# tests mutate a copy to trigger a specific violation.
make_matrix <- function() {
  set.seed(42)
  taxa_mat <- matrix(round(runif(8 * length(.FIXTURE_TAXA), 0, 0.2), 4),
                     nrow = 8, dimnames = list(NULL, .FIXTURE_TAXA))
  cbind(
    data.frame(
      id                    = c("P1","P1","P2","P2","P3","P3","P4","P4"),
      index                 = paste0("S", 1:8),
      preterm               = c(1, 1, 0, 0, 1, 1, 0, 0),
      sdg_parto             = c(35, 35, 39, 39, 34, 34, 40, 40),
      edad_cronologicamujer = c(31, 31, 28, 28, 35, 35, 26, 26),
      imc_pregestacional    = c(24, 24, 22, 22, 29, 29, 21, 21),
      stringsAsFactors      = FALSE),
    as.data.frame(taxa_mat, check.names = FALSE))
}

# Longitudinal metadata whose sample_id / subject_id keys match make_matrix().
make_metadata <- function() {
  data.frame(
    id        = c("P1","P1","P2","P2","P3","P3","P4","P4"),
    index     = paste0("S", 1:8),
    preterm   = c(1, 1, 0, 0, 1, 1, 0, 0),
    sdg_parto = c(35, 35, 39, 39, 34, 34, 40, 40),
    stringsAsFactors = FALSE)
}

# ---- CLR fixture -----------------------------------------------------------
# 6 samples x 4 taxa. Every taxon has >= 2 positive values (a hard requirement
# of zCompositions::cmultRepl's GBM estimator) and there are scattered zeros.
# Row 6 (0.5, 0.3, 0.1, 0.1) is deliberately zero-free -> used for the
# scale-invariance check. Row 1 has a zero -> used for the row-locality check.
CLR_TAXA <- c("Lactobacillus", "Gardnerella", "Escherichia-Shigella", "f__Rhizobiaceae")

make_clr_df <- function() {
  mat <- matrix(c(
    0.5, 0.3, 0.0, 0.2,
    0.6, 0.0, 0.3, 0.1,
    0.4, 0.4, 0.2, 0.0,
    0.7, 0.1, 0.0, 0.2,
    0.2, 0.5, 0.3, 0.0,
    0.5, 0.3, 0.1, 0.1),
    ncol = 4, byrow = TRUE, dimnames = list(NULL, CLR_TAXA))
  cbind(data.frame(index = paste0("S", 1:6), id = paste0("P", 1:6),
                   stringsAsFactors = FALSE),
        as.data.frame(mat, check.names = FALSE))
}

# fit_clr_zerorepl() prints cmultRepl's "No. adjusted imputations:" line to
# stdout; swallow it so test output stays clean.
fit_clr_quiet <- function(data, taxa_cols) {
  invisible(utils::capture.output(zl <- fit_clr_zerorepl(data, taxa_cols)))
  zl
}
