# optimize_threshold_cv(): pick the probability threshold that maximises Youden's
# J (sensitivity + specificity - 1) on a set of validation predictions. Leakage-
# wise it is a pure function of its input (it is called on inner-validation preds
# and never sees the outer-test fold). The J it maximises must be the J of the
# rule the engine applies afterwards, pred_prob >= threshold -> "1", including
# when preterm subjects score LOWER on the validation set (the case in which
# pROC's direction = "auto" would orient the curve the other way).

# Sensitivity, specificity and J of the applied rule (pred_prob >= t -> "1").
applied_rule <- function(vp, t) {
  y <- as.character(vp$true_class)
  sens <- mean(vp$pred_prob[y == "1"] >= t)
  spec <- mean(vp$pred_prob[y == "0"] <  t)
  c(sensitivity = sens, specificity = spec, J = sens + spec - 1)
}

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

test_that("optimize_threshold_cv scores the applied rule when preterm subjects score lower", {
  # Controls 0.30/0.50/0.60/0.80, cases 0.10/0.20/0.70: median(controls) 0.55 >
  # median(cases) 0.20, so direction = "auto" would orient the curve as ">".
  vp <- data.frame(
    true_class = factor(c("0", "0", "0", "0", "1", "1", "1"), levels = c("0", "1")),
    pred_prob  = c(0.30, 0.50, 0.60, 0.80, 0.10, 0.20, 0.70))
  r_auto <- pROC::roc(vp$true_class, vp$pred_prob, levels = c("0", "1"),
                      direction = "auto", quiet = TRUE)
  expect_identical(r_auto$direction, ">")                   # the fixture is the inverted case

  # Rule pred_prob >= t over the midpoints 0.15 ... 0.75 gives J = -1/3, -2/3,
  # -5/12, -1/6, 1/12, -1/4: the maximum is at t = 0.65 (sens 1/3, spec 3/4).
  # Under "auto" the maximum would be J = 2/3 at t = 0.25, whose J for the
  # applied rule is -2/3.
  res <- optimize_threshold_cv(vp, method = "youden")
  expect_equal(res$threshold, 0.65)
  expect_equal(res$sensitivity, 1 / 3)
  expect_equal(res$specificity, 3 / 4)
  expect_equal(res$criterion_value, 1 / 12)

  a <- applied_rule(vp, res$threshold)
  expect_equal(unname(a["sensitivity"]), res$sensitivity)
  expect_equal(unname(a["specificity"]), res$specificity)
  expect_equal(unname(a["J"]), res$criterion_value)
})

test_that("optimize_threshold_cv returns the maximum J of the applied rule (inner_val-sized sets)", {
  # 11 subjects with 4 preterm, as in the real inner-validation splits; scores
  # drawn independently of the labels, so about half the sets are inverted.
  runs <- do.call(rbind, lapply(1:40, function(s) {
    set.seed(1000 + s)
    vp <- data.frame(
      true_class = factor(sample(rep(c("0", "1"), c(7, 4))), levels = c("0", "1")),
      pred_prob  = round(runif(11), 3))
    sc <- sort(unique(vp$pred_prob))
    mids <- (head(sc, -1) + tail(sc, -1)) / 2                # pROC's finite thresholds
    res <- optimize_threshold_cv(vp, method = "youden")
    a <- applied_rule(vp, res$threshold)
    data.frame(
      best_J    = max(vapply(mids, function(t) applied_rule(vp, t)[["J"]], numeric(1))),
      criterion = res$criterion_value,
      J_applied = a[["J"]],
      sens = res$sensitivity, sens_applied = a[["sensitivity"]],
      spec = res$specificity, spec_applied = a[["specificity"]],
      inverted  = pROC::roc(vp$true_class, vp$pred_prob, levels = c("0", "1"),
                            direction = "auto", quiet = TRUE)$direction == ">")
  }))

  expect_equal(runs$criterion, runs$best_J)                 # the maximum over all thresholds
  expect_equal(runs$J_applied, runs$criterion)              # ... of the rule actually applied
  expect_equal(runs$sens_applied, runs$sens)
  expect_equal(runs$spec_applied, runs$spec)
  expect_gt(sum(runs$inverted), 5)                          # the inverted case is exercised
})
