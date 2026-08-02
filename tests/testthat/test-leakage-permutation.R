# ============================================================================
# End-to-end no-leakage: subject-level label permutation (#8).
# ----------------------------------------------------------------------------
# Why this exists: a pipeline WITH leakage still runs and still validates — the
# smoke test (#6) cannot detect it, it only inflates metrics. Leakage was a real
# bug in this project (threshold optimised on the test fold; feature selection
# outside the folds). So the no-leakage invariant is proven end-to-end here, not
# merely inferred from "the parts are pure + the pipeline runs".
#
# Mechanic: permute `preterm` across SUBJECTS (every sample of a subject shares
# its permuted label, preserving the longitudinal structure), run a REDUCED
# pipeline (Approach1_DREAM x ANCOM_Taxa x glmnet x 3 outer folds) and average
# over permutations. With the feature<->outcome relation destroyed, performance
# must collapse to chance.
#
# GATED on PTB_RUN_SLOW_TESTS (slow, ~1.5-2 min). This no-leakage test MUST run
# in CI: Chat 8's CI workflow HAS TO export PTB_RUN_SLOW_TESTS=1, otherwise the
# test is dormant — a skipped no-leakage test is green by omission, not by
# passing. It is deliberately NOT gated with skip_on_ci() for that exact reason.
# Off by default locally (set the env var to run it on demand).
# Runtime here ~9-10s per permutation; 10 permutations ~1.5-2 min.
# ============================================================================

test_that("permuted subject labels collapse the reduced pipeline to chance (end-to-end no leakage)", {
  # NB: intentionally NOT skip_on_ci() — CI is exactly where this must run
  # (Chat 8's workflow must set PTB_RUN_SLOW_TESTS=1). See the header note.
  skip_if(Sys.getenv("PTB_RUN_SLOW_TESTS") == "",
          "slow permutation test — set PTB_RUN_SLOW_TESTS=1 to run (Chat 8 CI must set this)")
  skip_if_not_installed("ANCOMBC")
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("PRROC")

  cfg <- read_config("example")
  ds  <- load_dataset(cfg)
  sp  <- split_domains(ds$matrix, cfg)
  micro_data <- sp$microbiome
  clin       <- sp$clinical

  # ---- minimal scaffolding, mirrored from the .Rmd chunks --------------------
  # prevalence >=5% + decontaminant filter -> genera_clean
  taxa_all    <- setdiff(names(micro_data), c("index", "id"))
  min_samples <- ceiling(nrow(micro_data) * cfg$preprocessing$prevalence_filter)
  present     <- vapply(micro_data[taxa_all], function(x) sum(x > 0), numeric(1))
  genera_keep  <- names(present)[present >= min_samples]
  genera_clean <- genera_keep[!genera_keep %in% cfg$preprocessing$contaminant_genera]

  # Shannon diversity (untransformed rel. abundances) -> micro_genus_full
  micro_mat <- as.matrix(micro_data[, taxa_all]); storage.mode(micro_mat) <- "numeric"
  diversity_data <- data.frame(index = micro_data$index,
                               shannon_diversity = vegan::diversity(micro_mat, index = "shannon"),
                               stringsAsFactors = FALSE)
  micro_genus_full <- dplyr::left_join(
    micro_data[, c("index", "id", genera_clean)], diversity_data, by = "index")

  # Approach 1 (DREAM): index/id/sdg_visita/edad_cronologicamujer/preterm
  clinical_approach1 <- clin[, c("index", "id", cfg$approaches$approach1_dream, "preterm")]

  # subject-level outcome table (unpermuted baseline)
  base_labels <- clinical_approach1 %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(preterm = ifelse(any(preterm == "1"), "1", "0"), .groups = "drop")

  # CLR zero-levels: outcome-blind, learned once on the full microbiome
  clr_taxa <- setdiff(names(micro_genus_full), c("index", "id", "shannon_diversity"))
  invisible(capture.output(
    clr_zero_levels <- fit_clr_zerorepl(micro_genus_full, clr_taxa)))

  abs_data    <- load_abs_matrix(cfg)
  glmnet_spec <- build_model_specs(cfg)$glmnet_base

  # train_with_nested_cv() reads `subject_labels` from the global env (a known
  # engine global-dep, flagged in R/nested_cv.R). Provide it per permutation and
  # restore whatever was there on exit.
  had_sl <- exists("subject_labels", envir = globalenv(), inherits = FALSE)
  old_sl <- if (had_sl) get("subject_labels", envir = globalenv()) else NULL
  on.exit({
    if (had_sl) assign("subject_labels", old_sl, envir = globalenv())
    else if (exists("subject_labels", envir = globalenv(), inherits = FALSE))
      rm("subject_labels", envir = globalenv())
  }, add = TRUE)

  run_one_permutation <- function(seed) {
    set.seed(seed)
    perm <- base_labels
    perm$preterm <- sample(perm$preterm)                 # shuffle labels across SUBJECTS
    lab <- stats::setNames(perm$preterm, perm$id)        # id -> permuted label

    # propagate the permuted label everywhere it is consumed (keeps it consistent
    # across the model outcome, the inner split, and ANCOM's selection input)
    clin_perm <- clinical_approach1
    clin_perm$preterm <- unname(lab[as.character(clin_perm$id)])
    abs_perm  <- abs_data
    abs_perm$meta$preterm <- unname(lab[as.character(abs_perm$meta$id)])

    assign("subject_labels", perm, envir = globalenv())  # engine reads this
    set.seed(1000 + seed)
    cv_folds <- rsample::vfold_cv(perm, v = 3, strata = preterm)

    res <- train_with_nested_cv(
      model_name = "glmnet_base", model_spec = glmnet_spec,
      clinical_data_all = clin_perm, microbiome_data_all = micro_genus_full,
      approach_name = "Approach1_DREAM", microbiome_option = "ANCOM_Taxa",
      cv_folds = cv_folds, clr_zero_levels = clr_zero_levels,
      genera_clean = genera_clean, abs_data = abs_perm)

    c(auroc  = as.numeric(res$summary$AUROC_mean),
      balacc = as.numeric(res$summary$Balanced_Accuracy_mean))
  }

  n_perm <- 10   # average over 10 perms to tame small-fold variance (fixed seeds -> reproducible)
  M <- matrix(NA_real_, nrow = n_perm, ncol = 2, dimnames = list(NULL, c("auroc", "balacc")))
  for (i in seq_len(n_perm)) {
    invisible(capture.output(M[i, ] <- run_one_permutation(200 + i)))
  }
  mean_auroc  <- mean(M[, "auroc"])
  mean_balacc <- mean(M[, "balacc"])
  message(sprintf("permutation null over %d perms: mean AUROC=%.3f  mean BalAcc=%.3f",
                  n_perm, mean_auroc, mean_balacc))

  # ---- PRIMARY (strict): a THRESHOLD-DEPENDENT metric ------------------------
  # With the labels destroyed, the Youden-on-inner-validation threshold applied
  # to the independent outer-test fold yields chance balanced accuracy (~0.5,
  # empirically ~0.52). This is the arm that would stay HIGH under the real
  # historical bug (threshold optimised on the test fold), because a leaked
  # threshold inflates sens/spec/accuracy even when AUROC does not.
  expect_lt(mean_balacc, 0.65)

  # ---- LOOSE GUARD: AUROC ----------------------------------------------------
  # The bound is 0.75, NOT ~0.5, on purpose — and this is NOT a threshold
  # loosened merely to pass. The engine computes AUROC with
  # pROC::roc(direction = "auto"), which orients each fold to whichever direction
  # gives AUC >= 0.5. Under random labels on these small outer folds that auto-
  # flip folds the NULL AUROC into [0.5, 1] with mean ~0.63 (a known property of
  # the metric, not leakage — it cannot collapse to 0.5 here). The real
  # threshold-leakage verification is carried by Balanced_Accuracy above, and
  # feature-selection leakage is covered deterministically at the component level
  # by test-ancom-leakage.R. This guard still trips on a GROSS feature-selection
  # leak, which would push the permuted AUROC well past 0.75.
  expect_lt(mean_auroc, 0.75)
})
