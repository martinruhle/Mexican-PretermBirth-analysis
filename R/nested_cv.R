# ============================================================================
# Nested cross-validation engine (subject-level, leakage-safe)
# ----------------------------------------------------------------------------
# Extracted from the chunk `nested_cv_helper_functions` of
# analysis/integrated_preterm_prediction_workflow.Rmd (Chat 2). The methodological
# logic is unchanged. Edits relative to the source:
#   * roxygen docs added.
#   * train_with_nested_cv() gains an explicit `clr_zero_levels` argument, threaded
#     into every apply_clr_transform() call, so the CLR zero-replacement levels are
#     passed in (learned once on the full microbiome) instead of being read from a
#     global default. (Chat 2 decision: "thread as arg, compute in driver".)
#
# ✓ RESOLVED in Chat 3 (now explicit arguments of train_with_nested_cv()):
#     - `genera_clean`   contaminant-filtered genus list -> `genera_clean` argument.
#     - the hardcoded absolute `abs_file` + positional indexing (row.names = 98,
#       [, 1:97], [, 98:166]) -> `abs_data` argument, produced once by
#       io.R::load_abs_matrix(cfg) (config abs_matrix_path, split by NAME).
#
# ⚠ GLOBAL-STATE DEPENDENCIES still present in train_with_nested_cv() — these are
#   read from the calling (global) environment and are flagged inline with
#   "GLOBAL DEP" comments to be resolved in later chats (config / models):
#     - `subject_labels`            (subject-level outcome table for the inner split)
#     - univariate-screening helpers used in the Approach-3 branch:
#         calculate_completeness, filter_by_completeness,
#         calculate_univariate_association, prioritize_continuous,
#         detect_collinear_vars, resolve_collinearity
#     - accumulator objects grown via exists()/<- inside the fold loop
#         (selected_taxa_by_fold, roc_curves_storage, pr_curves_storage)
#
# Depends on the analysis packages being attached (tidymodels/rsample: analysis,
# assessment, initial_split, training, testing, recipe/workflow/fit steps;
# yardstick: accuracy/sens/spec/pr_auc; pROC: roc/auc; PRROC: pr.curve; dplyr/tidyr/
# purrr), exactly as in the source .Rmd.
# ============================================================================

#' Detect the positive-class probability column in a predictions tibble
#'
#' Returns the name of the predicted-probability column for the positive class,
#' handling the tidymodels `.pred_<level>` naming and falling back to the first
#' `.pred_*` column (excluding `.pred_class`).
#'
#' @param df A predictions data.frame/tibble (e.g. from `predict(type = "prob")`).
#' @param pos_class The positive class label. Default `"1"`.
#'
#' @return The name of the probability column (character scalar).
#' @export
detect_prob_col <- function(df, pos_class = "1") {
  if(".pred_1" %in% names(df)) return(".pred_1")
  pcol <- paste0(".pred_", pos_class)
  if(pcol %in% names(df)) return(pcol)
  pcols <- grep("^\\.pred_", names(df), value = TRUE)
  pcols <- setdiff(pcols, ".pred_class")
  if (length(pcols) > 0) {
    return(pcols[1])
  } else {
    stop("No probability column found in predictions.")
  }
}

#' Train and evaluate one model configuration with nested cross-validation
#'
#' Runs the outer 5-fold, subject-level CV loop for a single
#' (model x approach x microbiome) combination: leakage-safe feature selection
#' inside each fold (ANCOM-BC2 taxa and/or Approach-3 univariate screening),
#' per-sample CLR, an inner 70/30 split to optimise the classification threshold,
#' evaluation on the independent outer-test fold, and finally a model refit on the
#' full data. Returns fold-level and summary metrics plus stored ROC/PR curves and
#' the final fitted model.
#'
#' @param model_name Character label for the model (e.g. `"rf_base"`).
#' @param model_spec A parsnip model specification.
#' @param clinical_data_all Clinical data for the approach (all subjects). For
#'   Approach 3 this is the full clinical table; variable selection happens inside
#'   the loop.
#' @param microbiome_data_all Microbiome data for all subjects (`index`, `id`,
#'   `shannon_diversity`, taxa columns).
#' @param approach_name One of `"Approach1_DREAM"`, `"Approach2_Literature"`,
#'   `"Approach3_DataDriven"`.
#' @param microbiome_option `"ANCOM_Taxa"` (taxa selected per fold) or
#'   `"Full_Microbiome"`.
#' @param cv_folds An rsample `vfold_cv` object of subject-level folds.
#' @param clr_zero_levels Named numeric vector of per-taxon CLR zero-replacement
#'   levels from [fit_clr_zerorepl()], learned once on the full microbiome and
#'   passed to every [apply_clr_transform()] call (threaded in explicitly rather
#'   than read from a global).
#' @param genera_clean Character vector of contaminant-filtered genus names
#'   (post prevalence + decontaminant filter). Required when
#'   `microbiome_option == "ANCOM_Taxa"`; the ANCOM OTU table is restricted to
#'   these genera. Passed in explicitly (Chat 3) instead of read from a global.
#' @param abs_data Absolute-count data for ANCOM-BC2 as a list with `otu`
#'   (samples x taxa) and `meta` (subject_id + clinical), from
#'   [load_abs_matrix()]. Required when `microbiome_option == "ANCOM_Taxa"`.
#'   Replaces the old hardcoded absolute path + positional indexing.
#' @param use_pca Logical; if `TRUE` and there are >10 taxa, adds a PCA step to the
#'   recipe. Default `FALSE`.
#'
#' @return A named list with `model_name`, `approach`, `microbiome`,
#'   `fold_results`, `summary`, and (when available) `test_predictions`,
#'   `roc_curves`, `pr_curves`, `selected_variables`, `selected_taxa` and
#'   `final_model`; or `NULL` if no fold produced results. `test_predictions` holds
#'   the subject-level out-of-fold predictions (`fold`, `id`, `pred_prob`,
#'   `true_class`, `pred_class`) used to build `fold_results`, so downstream code
#'   can recompute alternative metrics (e.g. a direction-fixed AUROC) without
#'   re-running the pipeline.
#' @export
train_with_nested_cv <- function(model_name, model_spec,
                                  clinical_data_all, microbiome_data_all,
                                  approach_name, microbiome_option,
                                  cv_folds, clr_zero_levels,
                                  genera_clean = NULL, abs_data = NULL,
                                  use_pca = FALSE) {

  cat(sprintf("\n========================================\n"))
  cat(sprintf("MODEL: %s | %s | %s\n",
              model_name, approach_name, microbiome_option))
  cat(sprintf("========================================\n"))

  # Store results from each fold
  fold_results <- list()

  # Subject-level out-of-fold test predictions, one entry per fold. ADDITIVE: no
  # metric computed below reads this; it is only packaged into the returned list
  # so downstream analyses (e.g. the permutation test) can recompute a
  # direction-fixed AUROC without re-running the pipeline. Explicit local list —
  # deliberately NOT the exists()/<- accumulator pattern used elsewhere here.
  test_predictions_by_fold <- list()

  # For Approach 3, we'll also store which variables were selected in each fold
  if(approach_name == "Approach3_DataDriven") {
    selected_variables_by_fold <- list()
  }

  # Iterate through outer CV folds
  for(fold_idx in 1:nrow(cv_folds)) {

    set.seed(123 + fold_idx)

    cat(sprintf("\n--- Outer Fold %d/%d ---\n", fold_idx, nrow(cv_folds)))

    # Get outer fold split (subject IDs)
    outer_train_subjects <- analysis(cv_folds$splits[[fold_idx]])$id
    outer_test_subjects <- assessment(cv_folds$splits[[fold_idx]])$id

    # ========================================================================
    # IMPROVED: UNIVARIATE SCREENING FOR APPROACH 3
    # ========================================================================
    if(approach_name == "Approach3_DataDriven") {

      cat("  Performing improved univariate screening...\n")

      # GLOBAL DEP: this branch calls univariate-screening helpers defined in
      # other .Rmd chunks (calculate_completeness, filter_by_completeness,
      # calculate_univariate_association, prioritize_continuous,
      # detect_collinear_vars, resolve_collinearity). Resolve/extract later.

      # Filter to OUTER TRAINING subjects only
      clinical_for_screening <- clinical_data_all %>%
        filter(id %in% outer_train_subjects)

      # Get candidate variables
      candidate_vars_all <- clinical_for_screening %>%
        select(-index, -id, -preterm) %>%
        names()

      # STEP 1: Filter by data completeness
      completeness_stats <- calculate_completeness(
        clinical_for_screening,
        candidate_vars_all
      )

      candidate_vars_complete <- filter_by_completeness(
        candidate_vars_all,
        completeness_stats,
        subject_threshold = 80,  # ≥80% subjects
        sample_threshold = 70    # ≥70% samples
      )

      cat(sprintf("  After completeness filter: %d variables\n",
                  length(candidate_vars_complete)))

      # STEP 2: Calculate univariate associations
      univariate_results_fold <- tibble(variable = candidate_vars_complete) %>%
        mutate(
          association = map(
            variable,
            ~calculate_univariate_association(.x, clinical_for_screening)
          )
        ) %>%
        unnest_wider(association) %>%
        filter(!is.na(p_value)) %>%
        arrange(p_value)

      # STEP 3: Prioritize continuous variables
      # Order by: continuous → numeric → categorical, then by p-value
      vars_by_type <- prioritize_continuous(
        univariate_results_fold$variable,
        clinical_for_screening
      )

      univariate_results_fold <- univariate_results_fold %>%
        mutate(priority_rank = match(variable, vars_by_type)) %>%
        arrange(priority_rank, p_value)

      # STEP 4: Select top variables (liberal p-value threshold)
      # Take top 15 candidates to allow for collinearity removal
      p_threshold <- 0.30  # Liberal threshold
      max_candidates <- 15

      top_candidates <- univariate_results_fold %>%
        filter(p_value < p_threshold) %>%
        head(max_candidates) %>%
        pull(variable)

      cat(sprintf("  Top candidates (p<%.2f): %d variables\n",
                  p_threshold, length(top_candidates)))

      # STEP 5: Detect and resolve collinearity
      collinear_pairs <- detect_collinear_vars(
        clinical_for_screening,
        top_candidates,
        threshold = 0.95
      )

      if(nrow(collinear_pairs) > 0) {
        cat(sprintf("  Found %d collinear pairs (|r|>0.85)\n",
                    nrow(collinear_pairs)))

        # Resolve by priority (continuous > others, lower p-value)
        candidate_vars_ordered <- univariate_results_fold %>%
          filter(variable %in% top_candidates) %>%
          arrange(priority_rank, p_value) %>%
          pull(variable)

        top_candidates <- resolve_collinearity(
          candidate_vars_ordered,
          collinear_pairs
        )
      }

      # STEP 6: Select final 10 variables
      max_features <- 10
      approach3_vars_fold <- head(top_candidates, max_features)

      cat(sprintf("  Final selection: %d variables\n", length(approach3_vars_fold)))

      # Show selected variables with details
      selected_info <- univariate_results_fold %>%
        filter(variable %in% approach3_vars_fold) %>%
        arrange(match(variable, approach3_vars_fold))

      cat("\n  Selected variables:\n")
      for(i in 1:nrow(selected_info)) {
        var_type <- if(is.numeric(clinical_for_screening[[selected_info$variable[i]]])) {
          n_unique <- length(unique(na.omit(clinical_for_screening[[selected_info$variable[i]]])))
          if(n_unique > 10) "continuous" else "numeric"
        } else "categorical"

        cat(sprintf("    %2d. %-30s p=%.4f  (%s)\n",
                    i,
                    selected_info$variable[i],
                    selected_info$p_value[i],
                    var_type))
      }
      cat("\n")

      # Store selected variables for this fold
      selected_variables_by_fold[[fold_idx]] <- selected_info %>%
       select(variable, p_value, coefficient, n_complete) %>%
       mutate(fold = fold_idx)

      # Create clinical dataset with ONLY selected variables
      clinical_data_fold <- clinical_data_all %>%
        select(index, id, all_of(approach3_vars_fold), preterm)

    } else {
      # For Approaches 1 and 2, use pre-defined clinical data
      clinical_data_fold <- clinical_data_all
    }

    # ========================================================================
    # ANCOM-BC2 ANALYSIS FOR MICROBIOME TAXA SELECTION
    # ========================================================================
    if(microbiome_option == "ANCOM_Taxa") {

      cat("  Performing ANCOM-BC2 on outer training subjects...\n")

      # Chat 3: la matriz de conteos absolutos llega PRECARGADA vía io.R
      # load_abs_matrix(cfg) (config abs_matrix_path, split por NOMBRE), en lugar
      # de una ruta absoluta hardcodeada + indexado posicional
      # (read.csv(row.names = 98), [, 1:97] OTU, [, 98:166] metadata).
      if(is.null(abs_data)) {
        stop("microbiome_option = 'ANCOM_Taxa' requiere abs_data = ",
             "load_abs_matrix(cfg); se recibió NULL.", call. = FALSE)
      }
      otu_table_raw <- abs_data$otu    # samples x 97 taxa (rownames = sample_id)
      meta_data_raw <- abs_data$meta   # subject_id + clínicas (mismas filas/orden)

      # ── Filtrar OTU table al mismo conjunto de géneros limpios (genera_clean) ────
      # genera_clean ahora entra por argumento (antes se leía del entorno global).
      if(is.null(genera_clean)) {
        stop("microbiome_option = 'ANCOM_Taxa' requiere genera_clean ",
             "(lista de géneros post-decontaminante).", call. = FALSE)
      }
      cols_to_keep <- intersect(genera_clean, colnames(otu_table_raw))
      otu_table_raw <- otu_table_raw[, cols_to_keep, drop = FALSE]
      cat(sprintf("  OTU table filtered to %d clean genera\n", ncol(otu_table_raw)))
      # ─────────────────────────────────────────────────────────────────────────────

      # Run ANCOM-BC2 on this fold's training subjects
      selected_taxa_fold <- tryCatch({
        run_ancombc_on_fold(otu_table_raw, meta_data_raw, outer_train_subjects)
      }, error = function(e) {
        cat("  ERROR in ANCOM-BC2:", e$message, "\n")
        return(character(0))
      })

      # Handle case where no taxa are significant
      if(length(selected_taxa_fold) == 0) {
        cat("  WARNING: No significant taxa found in this fold\n")
        cat("  Using top 5 taxa by mean abundance as fallback\n")

        # Fallback: use top abundant taxa in training set
        train_micro_data <- microbiome_data_all %>%
          filter(id %in% outer_train_subjects)

        taxa_cols_all <- setdiff(names(microbiome_data_all),
                                  c("index", "id", "shannon_diversity"))

        taxa_abundance <- train_micro_data %>%
          select(all_of(taxa_cols_all)) %>%
          summarize(across(everything(), ~mean(., na.rm = TRUE))) %>%
          pivot_longer(everything(), names_to = "taxon", values_to = "mean_abundance") %>%
          filter(taxon != "o__Chloroplast") %>%  # Exclude chloroplast
          arrange(desc(mean_abundance))

        selected_taxa_fold <- head(taxa_abundance$taxon, 5)
      }

      cat(sprintf("  Selected %d ANCOM taxa for this fold\n", length(selected_taxa_fold)))
      cat("  Taxa:", paste(selected_taxa_fold, collapse = ", "), "\n")

      # CRITICAL: Create microbiome_data_fold with ONLY selected taxa
      microbiome_data_fold <- microbiome_data_all %>%
        select(index, id, shannon_diversity, all_of(selected_taxa_fold))

      cat(sprintf("  Created microbiome_data_fold: %d samples × %d taxa\n",
                  nrow(microbiome_data_fold),
                  length(selected_taxa_fold)))

      # GLOBAL DEP: accumulator grown via exists()/<- across fold iterations
      # (function-local, but exists() searches enclosing envs). Refactor to an
      # explicit list initialised before the loop.
      # Store selected taxa for this fold (for later analysis)
      if(!exists("selected_taxa_by_fold")) {
        selected_taxa_by_fold <- list()
      }
      selected_taxa_by_fold[[fold_idx]] <- tibble(
        fold = fold_idx,
        taxon = selected_taxa_fold
      )

    } else {
      # Full Microbiome - clean up microbiome_data_fold if it exists
      if(exists("microbiome_data_fold")) {
        rm(microbiome_data_fold)
      }
    }
    # ========================================================================

    # Further split outer training into inner train/validation (70/30)
    # GLOBAL DEP: `subject_labels` (subject-level outcome table, built in the
    # nested_cv_setup chunk). Pass in as an argument later.
    outer_train_labels <- subject_labels %>%
      filter(id %in% outer_train_subjects)

    set.seed(123 + fold_idx + 1000)
    inner_split <- initial_split(outer_train_labels, prop = 0.70, strata = preterm)
    inner_train_subjects <- training(inner_split)$id
    inner_val_subjects <- testing(inner_split)$id

    cat(sprintf("  Inner train: %d subjects | Inner val: %d subjects | Outer test: %d subjects\n",
                length(inner_train_subjects),
                length(inner_val_subjects),
                length(outer_test_subjects)))

    # ========================================================================
    # PREPARE DATA (COMBINE CLINICAL + MICROBIOME)
    # ========================================================================
    # Use microbiome_data_fold if it exists (ANCOM path), otherwise microbiome_data_all
    if(exists("microbiome_data_fold")) {
      microbiome_to_use <- microbiome_data_fold
      cat("  Using ANCOM-selected taxa:", ncol(microbiome_data_fold) - 2, "genera\n")
    } else {
      microbiome_to_use <- microbiome_data_all
      cat("  Using full microbiome:", ncol(microbiome_data_all) - 2, "genera\n")
    }

    full_data <- clinical_data_fold %>%
      left_join(microbiome_to_use, by = c("index", "id")) %>%
      mutate(preterm = factor(preterm, levels = c("0", "1")))

    inner_train_data <- full_data %>% filter(id %in% inner_train_subjects)
    inner_val_data <- full_data %>% filter(id %in% inner_val_subjects)
    outer_test_data <- full_data %>% filter(id %in% outer_test_subjects)

    # ========================================================================
    # IDENTIFY TAXA COLUMNS (CRITICAL: BASED ON ACTUAL DATA)
    # ========================================================================
    # Get clinical variable names from clinical_data_fold
    clinical_var_names <- setdiff(names(clinical_data_fold),
                                   c("index", "id", "preterm"))

    # Taxa cols = all columns EXCEPT: index, id, shannon_diversity, preterm, clinical vars
    taxa_cols <- setdiff(
      names(full_data),
      c("index", "id", "shannon_diversity", "preterm", clinical_var_names)
    )

    cat(sprintf("  Taxa columns identified for CLR: %d\n", length(taxa_cols)))

    # ========================================================================
    # APPLY CLR TRANSFORMATION
    # ========================================================================

    if(length(taxa_cols) > 0) {
      inner_train_data <- apply_clr_transform(inner_train_data, taxa_cols, zero_levels = clr_zero_levels)
      inner_val_data <- apply_clr_transform(inner_val_data, taxa_cols, zero_levels = clr_zero_levels)
      outer_test_data <- apply_clr_transform(outer_test_data, taxa_cols, zero_levels = clr_zero_levels)
    }

    # Create recipe
    if(use_pca && length(taxa_cols) > 10) {
      rec <- recipe(preterm ~ ., data = inner_train_data) %>%
        update_role(index, id, new_role = "ID") %>%
        step_zv(all_predictors()) %>%
        step_corr(all_numeric_predictors(), threshold = 0.95) %>%
        step_novel(all_nominal_predictors(), new_level = "(new)") %>%
        step_other(all_nominal_predictors(), threshold = 0.01) %>%
        step_impute_median(all_numeric_predictors()) %>%
        step_impute_mode(all_nominal_predictors()) %>%
        step_normalize(all_numeric_predictors()) %>%
        step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
        step_pca(starts_with("CLR_") | starts_with("Genus_"),
                 num_comp = min(20, length(taxa_cols)), prefix = "PC_genus_") %>%
        step_nzv(all_predictors())
    } else {
      rec <- recipe(preterm ~ ., data = inner_train_data) %>%
        update_role(index, id, new_role = "ID") %>%
        step_zv(all_predictors()) %>%
        step_corr(all_numeric_predictors(), threshold = 0.95) %>%
        step_novel(all_nominal_predictors(), new_level = "(new)") %>%
        step_other(all_nominal_predictors(), threshold = 0.01) %>%
        step_impute_median(all_numeric_predictors()) %>%
        step_impute_mode(all_nominal_predictors()) %>%
        step_normalize(all_numeric_predictors()) %>%
        step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
        step_nzv(all_predictors())
    }

    # Create workflow and fit on inner training
    wf <- workflow() %>%
      add_recipe(rec) %>%
      add_model(model_spec)

    set.seed(123 + fold_idx + 2000)

    final_fit <- tryCatch({
      fit(wf, data = inner_train_data)
    }, error = function(e) {
      cat("  ERROR in fitting:", e$message, "\n")
      return(NULL)
    })

    if(is.null(final_fit)) {
      fold_results[[fold_idx]] <- NULL
      next
    }

    # ========================================================================
    # STEP 1: Optimize threshold on INNER VALIDATION set
    # ========================================================================

    val_preds_prob <- predict(final_fit, new_data = inner_val_data, type = "prob")
    prob_col <- detect_prob_col(val_preds_prob)

    val_preds_subject <- inner_val_data %>%
      select(id, preterm) %>%
      bind_cols(val_preds_prob) %>%
      group_by(id) %>%
      summarize(
        pred_prob = mean(!!sym(prob_col), na.rm = TRUE),
        true_class = first(preterm),
        .groups = "drop"
      ) %>%
      mutate(true_class = factor(true_class, levels = c("0", "1")))

    threshold_opt <- optimize_threshold_cv(val_preds_subject, method = "youden")

    cat(sprintf("  Optimal threshold (inner validation): %.4f (Youden=%.3f)\n",
                threshold_opt$threshold, threshold_opt$criterion_value))

    # ========================================================================
    # STEP 2: Apply optimized threshold to OUTER TEST fold
    # ========================================================================

    test_preds_prob <- predict(final_fit, new_data = outer_test_data, type = "prob")
    test_prob_col <- detect_prob_col(test_preds_prob)

    if(nrow(test_preds_prob) == 0) {
      cat("  ERROR: No predictions generated for test set\n")
      fold_results[[fold_idx]] <- NULL
      next
    }

    test_combined <- outer_test_data %>%
      select(id, preterm) %>%
      bind_cols(test_preds_prob)

    if(nrow(test_combined) == 0) {
      cat("  ERROR: bind_cols produced 0 rows\n")
      fold_results[[fold_idx]] <- NULL
      next
    }

    test_preds_subject <- test_combined %>%
      group_by(id) %>%
      summarize(
        pred_prob = mean(!!sym(test_prob_col), na.rm = TRUE),
        true_class = first(preterm),
        .groups = "drop"
      )

    if(nrow(test_preds_subject) == 0) {
      cat("  ERROR: Subject aggregation produced 0 subjects\n")
      fold_results[[fold_idx]] <- NULL
      next
    }

    test_preds_subject <- test_preds_subject %>%
      mutate(
        true_class = factor(true_class, levels = c("0", "1")),
        pred_class = factor(
          ifelse(pred_prob >= threshold_opt$threshold, "1", "0"),
          levels = c("0", "1")
        )
      )

    # Degenerate outer-test fold: a single outcome class present -> AUROC/sens/spec
    # are undefined and pROC aborts with "No case observation". Skip the fold, in
    # line with the other guards above, instead of killing the whole model run.
    # UNREACHABLE with the real labels (outer folds are stratified on `preterm`);
    # it only arises under subject-level label permutation, where stratification
    # no longer holds. Leaves every non-degenerate fold untouched.
    if(length(unique(as.character(test_preds_subject$true_class))) < 2) {
      cat("  SKIP: outer-test fold has a single outcome class - metrics undefined\n")
      fold_results[[fold_idx]] <- NULL
      next
    }

    # Calculate metrics
    acc <- accuracy(test_preds_subject, truth = true_class, estimate = pred_class)$.estimate
    sens <- sens(test_preds_subject, truth = true_class, estimate = pred_class,
                 event_level = "second")$.estimate
    spec <- spec(test_preds_subject, truth = true_class, estimate = pred_class,
                 event_level = "second")$.estimate

    bal_acc <- (sens + spec) / 2

    roc_test <- roc(test_preds_subject$true_class,
                    test_preds_subject$pred_prob,
                    levels = c("0", "1"),
                    direction = "auto",
                    quiet = TRUE)
    auroc <- auc(roc_test)

    pr_data_input <- test_preds_subject %>%
      rename(.pred_1 = pred_prob)
    prauc <- pr_auc(pr_data_input, truth = true_class, .pred_1)$.estimate

    cat(sprintf("  Outer test metrics: AUROC=%.3f, Sens=%.3f, Spec=%.3f, Acc=%.3f\n",
                auroc, sens, spec, acc))

    # Store fold results
    fold_results[[fold_idx]] <- tibble(
      fold = fold_idx,
      threshold = threshold_opt$threshold,
      n_test = nrow(test_preds_subject),
      n_ptb = sum(test_preds_subject$true_class == "1"),
      AUROC = as.numeric(auroc),
      PRAUC = prauc,
      Accuracy = acc,
      Balanced_Accuracy = bal_acc,
      Sensitivity = sens,
      Specificity = spec,
      Youden = sens + spec - 1
    )

    # ADDITIVE: keep this fold's subject-level test predictions (see the list
    # initialisation above). Stored AFTER the metrics so it can never influence
    # them; a fold skipped by any guard above contributes nothing here either.
    test_predictions_by_fold[[fold_idx]] <- test_preds_subject %>%
      mutate(fold = fold_idx)

        # ========================================================================
    # STORE ROC AND PR CURVES
    # ========================================================================

    # ROC curve
    roc_obj <- roc(test_preds_subject$true_class,
                   test_preds_subject$pred_prob,
                   levels = c("0", "1"),
                   direction = "auto")

    # GLOBAL DEP: accumulator grown via exists()/<- (see note above).
    if(!exists("roc_curves_storage")) {
      roc_curves_storage <- list()
    }

    roc_curves_storage[[fold_idx]] <- tibble(
      fold = fold_idx,
      sensitivity = roc_obj$sensitivities,
      specificity = roc_obj$specificities,
      threshold = roc_obj$thresholds
    )

    # PR curve
    pr_obj <- pr.curve(
      scores.class0 = test_preds_subject$pred_prob[test_preds_subject$true_class == "1"],
      scores.class1 = test_preds_subject$pred_prob[test_preds_subject$true_class == "0"],
      curve = TRUE
    )

    # GLOBAL DEP: accumulator grown via exists()/<- (see note above).
    if(!exists("pr_curves_storage")) {
      pr_curves_storage <- list()
    }

    pr_curves_storage[[fold_idx]] <- tibble(
      fold = fold_idx,
      recall = pr_obj$curve[, 1],
      precision = pr_obj$curve[, 2],
      threshold = pr_obj$curve[, 3]
    )

    cat(sprintf("  ✓ ROC curve: %d points | PR curve: %d points\n",
                length(roc_obj$sensitivities),
                nrow(pr_obj$curve)))

    # ========================================================================
    # CLEANUP: Remove fold-specific variables
    # ========================================================================
    if(exists("microbiome_data_fold")) {
      rm(microbiome_data_fold)
    }
    if(exists("microbiome_to_use")) {
      rm(microbiome_to_use)
    }
    # ========================================================================

  }  # End of fold loop

  # Aggregate results across folds
  fold_results_df <- bind_rows(fold_results)

  if(nrow(fold_results_df) == 0) {
    return(NULL)
  }

  # Calculate mean and SD across folds
  summary_results <- fold_results_df %>%
    summarize(
      n_folds = n(),
      across(c(threshold, AUROC, PRAUC, Accuracy, Balanced_Accuracy,
               Sensitivity, Specificity, Youden),
             list(mean = ~mean(., na.rm = TRUE),
                  sd = ~sd(., na.rm = TRUE)),
             .names = "{.col}_{.fn}")
    )

  cat("\n=== CROSS-VALIDATION SUMMARY ===\n")
  cat(sprintf("  AUROC: %.3f ± %.3f\n",
              summary_results$AUROC_mean, summary_results$AUROC_sd))
  cat(sprintf("  Sensitivity: %.3f ± %.3f\n",
              summary_results$Sensitivity_mean, summary_results$Sensitivity_sd))
  cat(sprintf("  Specificity: %.3f ± %.3f\n",
              summary_results$Specificity_mean, summary_results$Specificity_sd))
  cat(sprintf("  Accuracy: %.3f ± %.3f\n",
              summary_results$Accuracy_mean, summary_results$Accuracy_sd))
  cat(sprintf("  Threshold: %.3f ± %.3f\n",
              summary_results$threshold_mean, summary_results$threshold_sd))

    # ========================================================================
  # TRAIN FINAL MODEL ON FULL DATA
  # ========================================================================

  cat("\n╔════════════════════════════════════════╗\n")
  cat("║  TRAINING FINAL MODEL (Full Dataset)  ║\n")
  cat("╚════════════════════════════════════════╝\n\n")

  # Use appropriate clinical data
  if(approach_name == "Approach3_DataDriven") {
    clinical_data_full <- clinical_data_all  # All variables for Approach 3
  } else {
    clinical_data_full <- clinical_data_all
  }

  # ========================================================================
  # CRITICAL FIX: Use correct microbiome based on option
  # ========================================================================

  if(microbiome_option == "ANCOM_Taxa") {
    cat("  Microbiome: ANCOM Taxa (using fold-selected taxa)\n")

    # Get taxa selected across folds
    if(exists("selected_taxa_by_fold") && length(selected_taxa_by_fold) > 0) {

      # =================================================================
      # CRITICAL FIX: Convert list to tibble first
      # =================================================================
      selected_taxa_df <- bind_rows(selected_taxa_by_fold)
      # =================================================================

      # Count frequency across folds
      taxa_freq <- selected_taxa_df %>%
        group_by(taxon) %>%
        summarize(n_folds = n(), .groups = "drop") %>%
        arrange(desc(n_folds))

      cat(sprintf("  Taxa selected across %d folds:\n",
                  length(unique(selected_taxa_df$fold))))
      print(taxa_freq)

      # Use taxa in ≥3 folds (majority)
      n_total_folds <- length(unique(selected_taxa_df$fold))
      min_folds <- max(3, ceiling(n_total_folds / 2))

      selected_taxa_final <- taxa_freq %>%
        filter(n_folds >= min_folds) %>%
        pull(taxon)

      if(length(selected_taxa_final) == 0) {
        cat("  WARNING: No taxa in majority of folds, using all\n")
        selected_taxa_final <- unique(selected_taxa_df$taxon)
      }

      cat(sprintf("  Using %d ANCOM taxa (≥%d folds):\n",
                  length(selected_taxa_final), min_folds))
      cat(paste("   ", selected_taxa_final, collapse = "\n"), "\n")

      # Filter microbiome to only these taxa
      available_taxa <- intersect(selected_taxa_final, names(microbiome_data_all))

      if(length(available_taxa) == 0) {
        stop("ERROR: No ANCOM taxa found in microbiome_data_all")
      }

      microbiome_data_full <- microbiome_data_all %>%
        select(index, id, shannon_diversity, all_of(available_taxa))

      cat(sprintf("  ✓ Final model: %d ANCOM taxa\n\n", length(available_taxa)))

    } else {
      cat("  WARNING: No selected_taxa_by_fold found\n")
      cat("  Using full microbiome as fallback\n")
      microbiome_data_full <- microbiome_data_all
    }

  } else {
    # Full microbiome
    cat("  Microbiome: Full Microbiome (all taxa)\n")
    microbiome_data_full <- microbiome_data_all
  }
  # ========================================================================

  # Combine
  full_data_final <- clinical_data_full %>%
    left_join(microbiome_data_full, by = c("index", "id")) %>%
    mutate(preterm = factor(preterm, levels = c("0", "1")))

  # Identify taxa columns
  clinical_var_names_full <- setdiff(names(clinical_data_full),
                                      c("index", "id", "preterm"))

  taxa_cols_full <- setdiff(
    names(full_data_final),
    c("index", "id", "shannon_diversity", "preterm", clinical_var_names_full)
  )

  cat(sprintf("  Clinical variables: %d\n", length(clinical_var_names_full)))
  cat(sprintf("  Microbiome taxa: %d\n", length(taxa_cols_full)))

  # Apply CLR
  if(length(taxa_cols_full) > 0) {
    full_data_final <- apply_clr_transform(full_data_final, taxa_cols_full, zero_levels = clr_zero_levels)
  }

  # Create recipe (same structure as in folds)
  if(use_pca && length(taxa_cols_full) > 10) {
    rec_final <- recipe(preterm ~ ., data = full_data_final) %>%
      update_role(index, id, new_role = "ID") %>%
      step_zv(all_predictors()) %>%
      step_corr(all_numeric_predictors(), threshold = 0.95) %>%
      step_novel(all_nominal_predictors(), new_level = "(new)") %>%
      step_other(all_nominal_predictors(), threshold = 0.01) %>%
      step_impute_median(all_numeric_predictors()) %>%
      step_impute_mode(all_nominal_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
      step_pca(starts_with("CLR_") | starts_with("Genus_"),
               num_comp = min(20, length(taxa_cols_full)),
               prefix = "PC_genus_") %>%
      step_nzv(all_predictors())
  } else {
    rec_final <- recipe(preterm ~ ., data = full_data_final) %>%
      update_role(index, id, new_role = "ID") %>%
      step_zv(all_predictors()) %>%
      step_corr(all_numeric_predictors(), threshold = 0.95) %>%
      step_novel(all_nominal_predictors(), new_level = "(new)") %>%
      step_other(all_nominal_predictors(), threshold = 0.01) %>%
      step_impute_median(all_numeric_predictors()) %>%
      step_impute_mode(all_nominal_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
      step_nzv(all_predictors())
  }

  # Fit final model
  wf_final <- workflow() %>%
    add_recipe(rec_final) %>%
    add_model(model_spec)

  set.seed(123)  # Reproducibility (unificado a 123 con la semilla del pipeline; antes 42)

  final_model_fit <- tryCatch({
    fit(wf_final, data = full_data_final)
  }, error = function(e) {
    cat("  ✗ ERROR:", e$message, "\n")
    return(NULL)
  })

  if(!is.null(final_model_fit)) {
    cat("  ✓ Final model trained successfully\n")
  } else {
    cat("  ✗ Final model training failed\n")
  }

  # ========================================================================
  # RETURN COMPREHENSIVE RESULTS
  # ========================================================================

  cat("\n╔════════════════════════════════════════╗\n")
  cat("║      PACKAGING RESULTS                 ║\n")
  cat("╚════════════════════════════════════════╝\n\n")

  result_list <- list(
    model_name = model_name,
    approach = approach_name,
    microbiome = microbiome_option,
    fold_results = fold_results_df,
    summary = summary_results
  )

  # Add subject-level out-of-fold test predictions (ADDITIVE; see above)
  if(length(test_predictions_by_fold) > 0) {
    result_list$test_predictions <- bind_rows(test_predictions_by_fold)
    cat(sprintf("  ✓ Test predictions: %d folds, %d subject-level rows\n",
                sum(!vapply(test_predictions_by_fold, is.null, logical(1))),
                nrow(result_list$test_predictions)))
  }

  # Add ROC curves
  if(exists("roc_curves_storage") && length(roc_curves_storage) > 0) {
    result_list$roc_curves <- bind_rows(roc_curves_storage)
    cat(sprintf("  ✓ ROC curves: %d folds, %d total points\n",
                length(roc_curves_storage),
                nrow(result_list$roc_curves)))
  } else {
    cat("  ✗ WARNING: No ROC curves stored\n")
  }

  # Add PR curves
  if(exists("pr_curves_storage") && length(pr_curves_storage) > 0) {
    result_list$pr_curves <- bind_rows(pr_curves_storage)
    cat(sprintf("  ✓ PR curves: %d folds, %d total points\n",
                length(pr_curves_storage),
                nrow(result_list$pr_curves)))
  } else {
    cat("  ✗ WARNING: No PR curves stored\n")
  }

  # Add selected variables (Approach 3)
  if(approach_name == "Approach3_DataDriven" && exists("selected_variables_by_fold")) {
    result_list$selected_variables <- bind_rows(selected_variables_by_fold)
    cat(sprintf("  ✓ Selected variables: %d folds\n",
                length(selected_variables_by_fold)))
  }

  # Add selected taxa (ANCOM)
  if(microbiome_option == "ANCOM_Taxa" && exists("selected_taxa_by_fold")) {
    result_list$selected_taxa <- bind_rows(selected_taxa_by_fold)
    cat(sprintf("  ✓ Selected taxa: %d folds\n",
                length(selected_taxa_by_fold)))
  }

  # Add final model
  if(exists("final_model_fit") && !is.null(final_model_fit)) {
    result_list$final_model <- final_model_fit
    cat("  ✓ Final model object stored\n")
  } else {
    cat("  ✗ WARNING: No final model stored\n")
  }

  cat("\n")
  return(result_list)

}  # End of train_with_nested_cv function
