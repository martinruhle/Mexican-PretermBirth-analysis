# ============================================================================
# Centered Log-Ratio (CLR) transformation — COMPOSITIONAL, per SAMPLE
# ----------------------------------------------------------------------------
# Extracted verbatim from the chunk `nested_cv_helper_functions` of
# analysis/integrated_preterm_prediction_workflow.Rmd (Chat 2). The methodology
# is unchanged; the only edits are (a) roxygen docs, (b) `zero_levels` is now a
# required argument of apply_clr_transform() instead of defaulting to the global
# `clr_zero_levels`, and (c) the matrix product is fully qualified as `base::%*%`
# (see note inside apply_clr_transform).
#
# CLR closes each SAMPLE (row) by the geometric mean of that sample's own taxa:
#   clr(x_ij) = log(x_ij) - mean_k log(x_ik)   (each row then sums to 0).
# The previous version subtracted each TAXON's across-sample geometric mean (a
# per-feature centering, NOT compositional CLR): row-sums were not 0, and once the
# recipe's step_normalize ran it was mathematically identical to a plain
# log(x + 0.65) transform — i.e. the microbiome was never CLR-transformed.
#
# Zeros are replaced first with Bayesian-multiplicative estimates
# (zCompositions::cmultRepl, GBM). This matrix is 68–81% zeros, so cmultRepl's
# default deletion of high-zero columns/rows is switched off (z.delete = FALSE,
# z.warning = 1) to keep every taxon and sample. The replacement is UNSUPERVISED
# (it never sees `preterm`), so the per-taxon levels are learned once on the full
# microbiome (fit_clr_zerorepl -> clr_zero_levels): per-fold refitting is
# infeasible here because some taxa are all-zero inside a training split. The
# per-sample CLR step is stateless (row-local), so applying it to train/val/test
# separately is leakage-free by construction — resolving the invariant that the
# old per-split recomputation of the centering violated.
# ============================================================================

#' Learn per-taxon zero-replacement levels from a training matrix
#'
#' Learns Bayesian-multiplicative (GBM) zero-replacement levels with
#' [zCompositions::cmultRepl()] on the closed relative-abundance matrix. GBM
#' replacement on closed relative data is ~constant per taxon, so the result is
#' summarised as one replacement level per taxon (the median of the replaced
#' values for the zeros of that taxon). This is **outcome-blind** and must be fit
#' on TRAINING data only; the returned object is then passed to
#' [apply_clr_transform()] for train/val/test alike.
#'
#' @param train_data A data.frame containing the taxa columns (relative
#'   abundances). Only the columns named in `taxa_cols` are used.
#' @param taxa_cols Character vector of taxa column names to learn levels for.
#'
#' @return A named numeric vector of length `length(taxa_cols)` giving the
#'   per-taxon zero-replacement level, carrying an attribute `"floor"` (the
#'   smallest positive level, used as a fallback for taxa absent from `train_data`
#'   and later dropped by `step_zv`).
#'
#' @details Taxa with no non-zero value in `train_data` cannot get a data-driven
#'   level and are assigned the `"floor"` value; they are removed downstream by
#'   the recipe's zero-variance step.
#'
#' @seealso [apply_clr_transform()]
#' @export
fit_clr_zerorepl <- function(train_data, taxa_cols) {
  M <- as.matrix(train_data[, taxa_cols, drop = FALSE]); storage.mode(M) <- "numeric"
  present <- colSums(M) > 0
  levels  <- setNames(rep(NA_real_, length(taxa_cols)), taxa_cols)
  if (any(present)) {
    R <- as.matrix(suppressWarnings(zCompositions::cmultRepl(
           M[, present, drop = FALSE], label = 0, method = "GBM",
           output = "prop", z.delete = FALSE, z.warning = 1)))
    pres_idx <- which(present)
    for (jj in seq_along(pres_idx)) {
      z <- which(M[, pres_idx[jj]] == 0)
      if (length(z)) levels[pres_idx[jj]] <- stats::median(R[z, jj])
    }
  }
  floor_val <- suppressWarnings(min(levels[levels > 0], na.rm = TRUE))
  if (!is.finite(floor_val)) floor_val <- 1e-6
  levels[is.na(levels)] <- floor_val   # taxa absent in train (dropped later by step_zv)
  attr(levels, "floor") <- floor_val
  levels
}

#' Apply zero-replacement and per-sample CLR to a data.frame
#'
#' Replaces zeros in the taxa columns with the (training-derived) `zero_levels`
#' and then applies a per-**sample** centered log-ratio transform (each row is
#' centred by the geometric mean of its own taxa, so rows sum to 0). The step is
#' stateless per row, so calling it separately on inner-train / inner-val /
#' outer-test is leakage-free provided `zero_levels` was learned on training data
#' only (see [fit_clr_zerorepl()]).
#'
#' @param data A data.frame containing the taxa columns to transform (other
#'   columns are returned untouched).
#' @param taxa_cols Character vector of taxa column names to transform.
#' @param zero_levels Named numeric vector of per-taxon replacement levels from
#'   [fit_clr_zerorepl()] (carrying its `"floor"` attribute). Required — it is
#'   passed in explicitly rather than read from a global (Chat 2 change).
#'
#' @return `data` with the `taxa_cols` columns replaced by their CLR-transformed
#'   values; all other columns unchanged.
#'
#' @seealso [fit_clr_zerorepl()]
#' @export
apply_clr_transform <- function(data, taxa_cols, zero_levels) {
  M  <- as.matrix(data[, taxa_cols, drop = FALSE]); storage.mode(M) <- "numeric"
  rj <- zero_levels[taxa_cols]
  floor_val <- attr(zero_levels, "floor"); if (is.null(floor_val)) floor_val <- 1e-6
  rj[is.na(rj)] <- floor_val
  Z <- M == 0
  M <- M * (1 - as.numeric(base::`%*%`(Z, rj)))             # shrink observed parts (multiplicative)
  # NB: base:: is required. `compositions` is on the search path and masks %*%
  # with an S4 method that mishandles (matrix %*% vector) here; with `conflicted`
  # active (via tidymodels_prefer) the bare call errors and blocks headless
  # render. This is plain matrix-vector multiply, so base:: is both the intended
  # operator and numerically inert.
  rjmat <- matrix(rj, nrow(M), ncol(M), byrow = TRUE); M[Z] <- rjmat[Z]
  M[M <= 0] <- floor_val
  logM <- log(M)
  clr  <- logM - rowMeans(logM)                             # per-SAMPLE centering (rows sum to 0)
  colnames(clr) <- taxa_cols
  data[taxa_cols] <- as.data.frame(clr)
  data
}
