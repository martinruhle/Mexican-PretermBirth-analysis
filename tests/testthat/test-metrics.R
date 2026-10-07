# Outer-test metrics: roc_ptb() (AUROC with a fixed orientation, preterm = higher
# score) and prauc_ptb() (PR-AUC with preterm as the event). The hand-built cases
# below are chosen so that the alternative definitions — pROC direction = "auto"
# and yardstick's default event_level = "first" — give a DIFFERENT number, so a
# regression to either default fails here. Expected values are worked out by hand
# in the comments.

test_that("roc_ptb keeps the orientation fixed when the model ranks the classes backwards", {
  # controls ("0") score 0.6, 0.7, 0.4; cases ("1") score 0.3, 0.5, 0.2.
  # Of the 9 (case, control) pairs only (0.5, 0.4) has the case above the
  # control -> AUROC = 1/9 with preterm = higher score.
  truth <- factor(c("0", "0", "0", "1", "1", "1"), levels = c("0", "1"))
  score <- c(0.6, 0.7, 0.4, 0.3, 0.5, 0.2)

  r <- roc_ptb(truth, score)
  expect_identical(r$direction, "<")
  expect_equal(as.numeric(pROC::auc(r)), 1 / 9)

  # direction = "auto" sees median(controls) = 0.6 > median(cases) = 0.3, flips
  # to ">" and reports the complement, 8/9. This is the behaviour roc_ptb()
  # replaces; asserting it here documents why the two must not be confused.
  r_auto <- pROC::roc(truth, score, levels = c("0", "1"), direction = "auto", quiet = TRUE)
  expect_identical(r_auto$direction, ">")
  expect_equal(as.numeric(pROC::auc(r_auto)), 8 / 9)
})

test_that("roc_ptb matches direction = 'auto' when the ranking is already the right way round", {
  truth <- factor(c("0", "0", "0", "1", "1", "1"), levels = c("0", "1"))
  score <- 1 - c(0.6, 0.7, 0.4, 0.3, 0.5, 0.2)           # same data, ranking reversed
  expect_equal(as.numeric(pROC::auc(roc_ptb(truth, score))), 8 / 9)
})

test_that("roc_ptb accepts a character outcome and ignores its order of appearance", {
  # levels are imposed as c("0", "1") regardless of which label appears first
  truth <- c("1", "0", "1", "0")
  score <- c(0.8, 0.6, 0.4, 0.1)
  r <- roc_ptb(truth, score)
  expect_identical(r$direction, "<")
  expect_equal(as.numeric(pROC::auc(r)), 3 / 4)         # 3 of 4 pairs ranked right
})

test_that("prauc_ptb scores the preterm class, not the term class", {
  # Ranked by score: 0.8 (PTB), 0.6 (term), 0.4 (PTB), 0.1 (term).
  # PR curve for PTB (recall, precision): (0, 1) (0.5, 1) (0.5, 0.5) (1, 2/3) (1, 0.5)
  # trapezoids: 0.5 * 1 + 0.5 * (0.5 + 2/3) / 2 = 1/2 + 7/24 = 19/24.
  truth <- factor(c("1", "0", "1", "0"), levels = c("0", "1"))
  score <- c(0.8, 0.6, 0.4, 0.1)
  expect_equal(prauc_ptb(truth, score), 19 / 24)

  # yardstick's default (event_level = "first") treats term ("0") as the event
  # while still ranking by P(preterm): (0, 1) (0, 0) (0.5, 0.5) (0.5, 1/3) (1, 0.5)
  # -> 0.5 * 0.5 / 2 + 0.5 * (1/3 + 0.5) / 2 = 1/8 + 5/24 = 1/3. A different number.
  d <- data.frame(true_class = truth, .pred_1 = score)
  expect_equal(yardstick::pr_auc(d, truth = true_class, .pred_1)$.estimate, 1 / 3)
})

test_that("prauc_ptb is 1 for a perfect preterm ranking and falls when the ranking flips", {
  truth <- factor(c("0", "0", "0", "1", "1", "1"), levels = c("0", "1"))
  expect_equal(prauc_ptb(truth, c(0.1, 0.2, 0.3, 0.7, 0.8, 0.9)), 1)
  expect_lt(prauc_ptb(truth, c(0.9, 0.8, 0.7, 0.3, 0.2, 0.1)), 0.5)  # prevalence is 0.5
})

# ----------------------------------------------------------------------------
# Wiring: train_with_nested_cv() must report exactly these metrics per fold, and
# the ROC curve it stores must be the curve behind the reported AUROC.
# ----------------------------------------------------------------------------
test_that("train_with_nested_cv reports roc_ptb / prauc_ptb and stores the matching ROC curve", {
  skip_if(Sys.getenv("PTB_RUN_SLOW_TESTS") == "",
          "runs the engine on the example data — set PTB_RUN_SLOW_TESTS=1 to run")
  skip_if_not_installed("PRROC")

  cfg <- read_config("example")
  ds  <- load_dataset(cfg)
  sp  <- split_domains(ds$matrix, cfg)
  micro_data <- sp$microbiome
  clin       <- sp$clinical

  # scaffolding mirrored from the .Rmd chunks (see test-leakage-permutation.R)
  taxa_all    <- setdiff(names(micro_data), c("index", "id"))
  min_samples <- ceiling(nrow(micro_data) * cfg$preprocessing$prevalence_filter)
  present     <- vapply(micro_data[taxa_all], function(x) sum(x > 0), numeric(1))
  genera_keep  <- names(present)[present >= min_samples]
  genera_clean <- genera_keep[!genera_keep %in% cfg$preprocessing$contaminant_genera]

  micro_mat <- as.matrix(micro_data[, taxa_all]); storage.mode(micro_mat) <- "numeric"
  diversity_data <- data.frame(index = micro_data$index,
                               shannon_diversity = vegan::diversity(micro_mat, index = "shannon"),
                               stringsAsFactors = FALSE)
  micro_genus_full <- dplyr::left_join(
    micro_data[, c("index", "id", genera_clean)], diversity_data, by = "index")

  clinical_approach1 <- clin[, c("index", "id", cfg$approaches$approach1_dream, "preterm")]
  labels <- clinical_approach1 %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(preterm = ifelse(any(preterm == "1"), "1", "0"), .groups = "drop")

  clr_taxa <- setdiff(names(micro_genus_full), c("index", "id", "shannon_diversity"))
  invisible(capture.output(
    clr_zero_levels <- fit_clr_zerorepl(micro_genus_full, clr_taxa)))

  # the engine reads `subject_labels` from the global env (known global-dep)
  had_sl <- exists("subject_labels", envir = globalenv(), inherits = FALSE)
  old_sl <- if (had_sl) get("subject_labels", envir = globalenv()) else NULL
  on.exit({
    if (had_sl) assign("subject_labels", old_sl, envir = globalenv())
    else if (exists("subject_labels", envir = globalenv(), inherits = FALSE))
      rm("subject_labels", envir = globalenv())
  }, add = TRUE)
  assign("subject_labels", labels, envir = globalenv())

  set.seed(123)
  cv_folds <- rsample::vfold_cv(labels, v = 3, strata = preterm)

  # Full_Microbiome: no ANCOM-BC2, so this stays a few seconds
  invisible(capture.output(res <- train_with_nested_cv(
    model_name = "glmnet_base", model_spec = build_model_specs(cfg)$glmnet_base,
    clinical_data_all = clinical_approach1, microbiome_data_all = micro_genus_full,
    approach_name = "Approach1_DREAM", microbiome_option = "Full_Microbiome",
    cv_folds = cv_folds, clr_zero_levels = clr_zero_levels,
    dict = read_data_dictionary(cfg), allowed_availability = "at_visit")))

  expect_false(is.null(res))
  fr <- res$fold_results
  tp <- res$test_predictions
  expect_setequal(unique(tp$fold), fr$fold)

  for (f in fr$fold) {
    d <- tp[tp$fold == f, ]
    expect_equal(fr$AUROC[fr$fold == f], as.numeric(pROC::auc(roc_ptb(d$true_class, d$pred_prob))))
    expect_equal(fr$PRAUC[fr$fold == f], prauc_ptb(d$true_class, d$pred_prob))

    # stored ROC curve (pROC order: threshold ascending) -> trapezoidal AUC must
    # reproduce the reported fold AUROC, i.e. the curve has the same orientation
    rc <- res$roc_curves[res$roc_curves$fold == f, ]
    auc_curve <- sum(diff(rc$specificity) *
                       (head(rc$sensitivity, -1) + tail(rc$sensitivity, -1)) / 2)
    expect_equal(auc_curve, fr$AUROC[fr$fold == f])
  }
})
