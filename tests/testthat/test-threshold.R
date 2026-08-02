# optimize_threshold_cv(): pick the probability threshold that maximises Youden's
# J (sensitivity + specificity - 1) on a set of validation predictions. Leakage-
# wise it is a pure function of its input (it is called on inner-validation preds
# and never sees the outer-test fold).

test_that("optimize_threshold_cv returns the Youden-optimal operating point on a separable set", {
  vp <- data.frame(
    true_class = factor(c("0", "0", "0", "1", "1", "1"), levels = c("0", "1")),
    pred_prob  = c(0.1, 0.2, 0.3, 0.7, 0.8, 0.9))          # perfectly separable at ~0.5

  res <- optimize_threshold_cv(vp, method = "youden")

  expect_equal(res$sensitivity, 1)                          # J = 1 + 1 - 1 = 1
  expect_equal(res$specificity, 1)
  expect_equal(res$criterion_value, 1)
  expect_gt(res$threshold, 0.3)                             # lands in the (0.3, 0.7) gap
  expect_lt(res$threshold, 0.7)
})

test_that("optimize_threshold_cv is deterministic (pure function of its input)", {
  vp <- data.frame(
    true_class = factor(c("0", "0", "1", "0", "1", "1"), levels = c("0", "1")),
    pred_prob  = c(0.20, 0.35, 0.40, 0.55, 0.75, 0.85))
  expect_identical(optimize_threshold_cv(vp), optimize_threshold_cv(vp))
})
