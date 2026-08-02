# fit_clr_zerorepl() + apply_clr_transform(): compositional, per-sample CLR with
# training-derived zero replacement. The no-leakage guarantee here is row-local:
# each sample is transformed using only its own taxa + the (train-learned)
# zero_levels, so the transform can be applied to inner-train / inner-val /
# outer-test separately without leakage.

test_that("apply_clr_transform centres each sample (rows sum to ~0), stays finite, keeps keys", {
  df <- make_clr_df()
  zl <- fit_clr_quiet(df, CLR_TAXA)
  res <- apply_clr_transform(df, CLR_TAXA, zl)

  expect_equal(unname(rowSums(res[, CLR_TAXA])), rep(0, nrow(df)), tolerance = 1e-8)
  expect_true(all(is.finite(as.matrix(res[, CLR_TAXA]))))   # zeros replaced -> no -Inf/NaN
  expect_identical(res$index, df$index)                     # non-taxa columns untouched
  expect_identical(res$id, df$id)
})

test_that("apply_clr_transform is scale-invariant on a zero-free sample: CLR(x) == CLR(c*x)", {
  df <- make_clr_df()
  zl <- fit_clr_quiet(df, CLR_TAXA)
  res <- apply_clr_transform(df, CLR_TAXA, zl)

  df_scaled <- df
  df_scaled[6, CLR_TAXA] <- df_scaled[6, CLR_TAXA] * 10     # row 6 has no zeros
  res_scaled <- apply_clr_transform(df_scaled, CLR_TAXA, zl)

  expect_equal(as.numeric(res_scaled[6, CLR_TAXA]),
               as.numeric(res[6, CLR_TAXA]), tolerance = 1e-8)
})

test_that("apply_clr_transform is row-local: perturbing OTHER samples cannot change a target row", {
  df <- make_clr_df()
  zl <- fit_clr_quiet(df, CLR_TAXA)
  res <- apply_clr_transform(df, CLR_TAXA, zl)

  df_perturbed <- df
  df_perturbed[2:6, CLR_TAXA] <- df_perturbed[2:6, CLR_TAXA] * 3 + 0.01   # rows 2..6 changed
  res_perturbed <- apply_clr_transform(df_perturbed, CLR_TAXA, zl)        # SAME zero_levels

  # row 1 (kept identical) transforms identically -> no cross-sample coupling
  expect_equal(as.numeric(res_perturbed[1, CLR_TAXA]), as.numeric(res[1, CLR_TAXA]))
})

test_that("fit_clr_zerorepl is deterministic and a pure function of train_data", {
  df <- make_clr_df()
  expect_identical(fit_clr_quiet(df, CLR_TAXA), fit_clr_quiet(df, CLR_TAXA))
})

test_that("fit_clr_zerorepl assigns the floor level to taxa absent from train_data", {
  df <- make_clr_df()
  df_absent <- cbind(df, Absent = 0)                       # a taxon with no positive value
  zl <- fit_clr_quiet(df_absent, c(CLR_TAXA, "Absent"))

  expect_equal(unname(zl["Absent"]), attr(zl, "floor"))    # falls back to the floor
  res <- apply_clr_transform(df_absent, c(CLR_TAXA, "Absent"), zl)
  expect_true(all(is.finite(res[["Absent"]])))             # still finite (dropped later by step_zv)
})
