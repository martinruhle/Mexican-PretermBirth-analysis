#!/usr/bin/env Rscript
# =============================================================================
# Sensitivity Analysis: Non-Independence of Longitudinal Samples
# =============================================================================
#
# Project:  Vaginal Microbiome-Based Preterm Birth Prediction in Mexican Women
# Authors:  Ruhle M, et al.
# Purpose:  Address reviewer concern that multiple samples per subject violate
#           independence assumptions and that differential sampling depth
#           (term: 3.0 samples/subject vs PTB: 1.7) may bias feature
#           identification and model performance.
#
# Analysis:
#   Inverse-frequency sample weighting: Retain all 110 samples but weight
#   each sample by 1/n_samples_per_subject, downweighting heavily-sampled
#   individuals during model training so every subject contributes equally (weight=1).
#
# Prerequisites:
#   - Run the main analysis pipeline first (integrated_preterm_prediction_workflow.Rmd)
#   - Required objects in environment:
#       complete_data, approach3_clinical_all, micro_genus_full, 
#       subject_labels, cv_folds, and all helper functions
#
# Last updated: 2026-03-11
# =============================================================================

# --- Setup -------------------------------------------------------------------

library(tidyverse)
library(tidymodels)
library(ranger)
library(pROC)
library(yardstick)

cat(
  "=============================================================\n",
  "  SENSITIVITY ANALYSIS: NON-INDEPENDENCE OF SAMPLES\n",
  "  (Inverse-Frequency Sample Weighting)\n",
  "=============================================================\n\n"
)

# Output directory
output_dir <- "sensitivity_nonindependence"
if (!dir.exists(output_dir)) dir.create(output_dir)


# =============================================================================
# SECTION 1: CHARACTERIZE THE SAMPLING IMBALANCE
# =============================================================================

cat("--- Section 1: Sampling structure characterization ---\n\n")

# Document the differential sampling that motivates this analysis
sampling_summary <- complete_data %>%
  group_by(id) %>%
  summarize(
    n_samples   = n(),
    preterm     = first(preterm),
    .groups     = "drop"
  ) %>%
  mutate(
    outcome = ifelse(preterm == 1, "PTB", "Term")
  )

# Per-group summary
group_sampling <- sampling_summary %>%
  group_by(outcome) %>%
  summarize(
    n_subjects      = n(),
    total_samples   = sum(n_samples),
    mean_samples    = mean(n_samples),
    median_samples  = median(n_samples),
    min_samples     = min(n_samples),
    max_samples     = max(n_samples),
    .groups = "drop"
  )

cat("  Sampling structure by outcome group:\n\n")
print(group_sampling)

cat(sprintf(
  "\n  Sample-level class ratio: %d term : %d PTB = %.1f:1\n",
  group_sampling$total_samples[group_sampling$outcome == "Term"],
  group_sampling$total_samples[group_sampling$outcome == "PTB"],
  group_sampling$total_samples[group_sampling$outcome == "Term"] /
    group_sampling$total_samples[group_sampling$outcome == "PTB"]
))
cat(sprintf(
  "  Subject-level class ratio: %d term : %d PTB = %.1f:1\n\n",
  group_sampling$n_subjects[group_sampling$outcome == "Term"],
  group_sampling$n_subjects[group_sampling$outcome == "PTB"],
  group_sampling$n_subjects[group_sampling$outcome == "Term"] /
    group_sampling$n_subjects[group_sampling$outcome == "PTB"]
))


# =============================================================================
# SECTION 2: INVERSE-FREQUENCY SAMPLE WEIGHTING
# =============================================================================

cat("\n--- Section 2: Inverse-frequency sample weighting ---\n\n")

# Strategy: Each sample receives weight = 1 / n_samples_from_that_subject
# This ensures that each subject contributes equally regardless of the number
# of longitudinal samples.

# ---------------------------------------------------------------------------
# 2a. Compute sample weights
# ---------------------------------------------------------------------------

sample_weights_df <- complete_data %>%
  group_by(id) %>%
  mutate(
    n_samples_subject = n(),
    sample_weight = 1 / n_samples_subject
  ) %>%
  ungroup() %>%
  select(index, id, n_samples_subject, sample_weight)

cat("  Sample weight distribution:\n")
weight_summary <- sample_weights_df %>%
  group_by(n_samples_subject) %>%
  summarize(
    n_samples = n(),
    weight = first(sample_weight),
    .groups = "drop"
  )
for (i in 1:nrow(weight_summary)) {
  cat(sprintf("    %d sample(s)/subject → weight = %.3f (%d samples)\n",
              weight_summary$n_samples_subject[i],
              weight_summary$weight[i],
              weight_summary$n_samples[i]))
}

# Verify: sum of weights per subject should equal 1
weight_check <- sample_weights_df %>%
  group_by(id) %>%
  summarize(total_weight = sum(sample_weight), .groups = "drop")
cat(sprintf("\n  Weight verification: all subjects contribute total weight = %.1f ✓\n\n",
            unique(weight_check$total_weight)))

# ---------------------------------------------------------------------------
# 2b. Modified training function with case weights
# ---------------------------------------------------------------------------
# We create a wrapper that adds case weights to the ranger model.
# The weights are computed from training data only (within each fold).

train_with_nested_cv_weighted <- function(model_name, model_spec,
                                          clinical_data_all, microbiome_data_all,
                                          approach_name, microbiome_option,
                                          cv_folds, weights_df, use_pca = FALSE) {
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("MODEL: %s | %s | %s (WEIGHTED)\n",
              model_name, approach_name, microbiome_option))
  cat(sprintf("========================================\n"))
  
  fold_results <- list()
  
  if (approach_name == "Approach3_DataDriven") {
    selected_variables_by_fold <- list()
  }
  
  for (fold_idx in 1:nrow(cv_folds)) {
    
    set.seed(27 + fold_idx)
    
    cat(sprintf("\n--- Outer Fold %d/%d ---\n", fold_idx, nrow(cv_folds)))
    
    outer_train_subjects <- analysis(cv_folds$splits[[fold_idx]])$id
    outer_test_subjects  <- assessment(cv_folds$splits[[fold_idx]])$id
    
    # ---- Approach 3 univariate screening ----
    if (approach_name == "Approach3_DataDriven") {
      
      cat("  Performing univariate screening...\n")
      clinical_for_screening <- clinical_data_all %>%
        filter(id %in% outer_train_subjects)
      
      candidate_vars_all <- clinical_for_screening %>%
        select(-index, -id, -preterm) %>%
        names()
      
      completeness_stats <- calculate_completeness(
        clinical_for_screening, candidate_vars_all
      )
      candidate_vars_complete <- filter_by_completeness(
        candidate_vars_all, completeness_stats,
        subject_threshold = 80, sample_threshold = 70
      )
      
      univariate_results_fold <- tibble(variable = candidate_vars_complete) %>%
        mutate(
          association = map(variable,
                            ~calculate_univariate_association(.x, clinical_for_screening))
        ) %>%
        unnest_wider(association) %>%
        filter(!is.na(p_value)) %>%
        arrange(p_value)
      
      vars_by_type <- prioritize_continuous(
        univariate_results_fold$variable, clinical_for_screening
      )
      univariate_results_fold <- univariate_results_fold %>%
        mutate(priority_rank = match(variable, vars_by_type)) %>%
        arrange(priority_rank, p_value)
      
      p_threshold <- 0.30
      max_candidates <- 15
      top_candidates <- univariate_results_fold %>%
        filter(p_value < p_threshold) %>%
        head(max_candidates) %>%
        pull(variable)
      
      collinear_pairs <- detect_collinear_vars(
        clinical_for_screening, top_candidates, threshold = 0.95
      )
      if (nrow(collinear_pairs) > 0) {
        candidate_vars_ordered <- univariate_results_fold %>%
          filter(variable %in% top_candidates) %>%
          arrange(priority_rank, p_value) %>%
          pull(variable)
        top_candidates <- resolve_collinearity(
          candidate_vars_ordered, collinear_pairs
        )
      }
      
      approach3_vars_fold <- head(top_candidates, 10)
      cat(sprintf("  Selected %d variables\n", length(approach3_vars_fold)))
      
      selected_variables_by_fold[[fold_idx]] <- univariate_results_fold %>%
        filter(variable %in% approach3_vars_fold) %>%
        select(variable, p_value, coefficient, n_complete) %>%
        mutate(fold = fold_idx)
      
      clinical_data_fold <- clinical_data_all %>%
        select(index, id, all_of(approach3_vars_fold), preterm)
    } else {
      clinical_data_fold <- clinical_data_all
    }
    
    # ---- Prepare combined data ----
    full_data <- clinical_data_fold %>%
      left_join(microbiome_data_all, by = c("index", "id")) %>%
      mutate(preterm = factor(preterm, levels = c("0", "1")))
    
    # ---- Add case weights ----
    full_data <- full_data %>%
      left_join(weights_df %>% select(index, sample_weight), by = "index")
    
    # ---- Inner split ----
    outer_train_labels <- subject_labels %>%
      filter(id %in% outer_train_subjects)
    
    set.seed(27 + fold_idx + 1000)
    inner_split <- initial_split(outer_train_labels, prop = 0.70, strata = preterm)
    inner_train_subjects <- training(inner_split)$id
    inner_val_subjects   <- testing(inner_split)$id
    
    inner_train_data <- full_data %>% filter(id %in% inner_train_subjects)
    inner_val_data   <- full_data %>% filter(id %in% inner_val_subjects)
    outer_test_data  <- full_data %>% filter(id %in% outer_test_subjects)
    
    cat(sprintf("  Inner train: %d samples | Inner val: %d | Outer test: %d\n",
                nrow(inner_train_data), nrow(inner_val_data), nrow(outer_test_data)))
    
    # ---- CLR transformation ----
    clinical_var_names <- setdiff(names(clinical_data_fold),
                                  c("index", "id", "preterm"))
    taxa_cols <- setdiff(
      names(full_data),
      c("index", "id", "shannon_diversity", "preterm",
        clinical_var_names, "sample_weight")
    )
    
    if (length(taxa_cols) > 0) {
      inner_train_data <- apply_clr_transform(inner_train_data, taxa_cols)
      inner_val_data   <- apply_clr_transform(inner_val_data, taxa_cols)
      outer_test_data  <- apply_clr_transform(outer_test_data, taxa_cols)
    }
    
    # ---- Extract weights for training samples ----
    train_weights <- inner_train_data$sample_weight
    
    # ---- Remove weight column before recipe ----
    inner_train_data <- inner_train_data %>% select(-sample_weight)
    inner_val_data   <- inner_val_data %>% select(-sample_weight)
    outer_test_data  <- outer_test_data %>% select(-sample_weight)
    
    # ---- Recipe ----
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
    
    # ---- Model spec with case weights passed to ranger ----
    weighted_rf_spec <- rand_forest(
      mtry = 4, trees = 500, min_n = 10
    ) %>%
      set_mode("classification") %>%
      set_engine("ranger",
                 importance = "impurity",
                 num.threads = 1,
                 case.weights = train_weights)
    
    # ---- Workflow and fit ----
    wf <- workflow() %>%
      add_recipe(rec) %>%
      add_model(weighted_rf_spec)
    
    set.seed(27 + fold_idx + 2000)
    
    final_fit <- tryCatch({
      fit(wf, data = inner_train_data)
    }, error = function(e) {
      cat("  ERROR:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(final_fit)) {
      fold_results[[fold_idx]] <- NULL
      next
    }
    
    # ---- Threshold optimization on inner validation ----
    val_preds_prob <- predict(final_fit, new_data = inner_val_data, type = "prob")
    prob_col <- detect_prob_col(val_preds_prob)
    
    val_preds_subject <- inner_val_data %>%
      select(id, preterm) %>%
      bind_cols(val_preds_prob) %>%
      group_by(id) %>%
      summarize(
        pred_prob  = mean(!!sym(prob_col), na.rm = TRUE),
        true_class = first(preterm),
        .groups = "drop"
      ) %>%
      mutate(true_class = factor(true_class, levels = c("0", "1")))
    
    threshold_opt <- optimize_threshold_cv(val_preds_subject, method = "youden")
    
    cat(sprintf("  Threshold: %.4f (Youden=%.3f)\n",
                threshold_opt$threshold, threshold_opt$criterion_value))
    
    # ---- Evaluate on outer test ----
    test_preds_prob <- predict(final_fit, new_data = outer_test_data, type = "prob")
    test_prob_col <- detect_prob_col(test_preds_prob)
    
    test_preds_subject <- outer_test_data %>%
      select(id, preterm) %>%
      bind_cols(test_preds_prob) %>%
      group_by(id) %>%
      summarize(
        pred_prob  = mean(!!sym(test_prob_col), na.rm = TRUE),
        true_class = first(preterm),
        .groups = "drop"
      ) %>%
      mutate(
        true_class = factor(true_class, levels = c("0", "1")),
        pred_class = factor(
          ifelse(pred_prob >= threshold_opt$threshold, "1", "0"),
          levels = c("0", "1")
        )
      )
    
    # Compute metrics
    acc  <- accuracy(test_preds_subject, truth = true_class, estimate = pred_class)$.estimate
    sens_val <- sens(test_preds_subject, truth = true_class, estimate = pred_class,
                     event_level = "second")$.estimate
    spec_val <- spec(test_preds_subject, truth = true_class, estimate = pred_class,
                     event_level = "second")$.estimate
    bal_acc <- (sens_val + spec_val) / 2
    
    roc_test <- pROC::roc(test_preds_subject$true_class,
                          test_preds_subject$pred_prob,
                          levels = c("0", "1"), direction = "auto", quiet = TRUE)
    auroc_val <- as.numeric(pROC::auc(roc_test))
    
    pr_data_input <- test_preds_subject %>% rename(.pred_1 = pred_prob)
    prauc_val <- pr_auc(pr_data_input, truth = true_class, .pred_1)$.estimate
    
    fold_results[[fold_idx]] <- tibble(
      fold              = fold_idx,
      threshold         = threshold_opt$threshold,
      n_test            = nrow(test_preds_subject),
      n_ptb             = sum(test_preds_subject$true_class == "1"),
      AUROC             = auroc_val,
      PRAUC             = prauc_val,
      Accuracy          = acc,
      Balanced_Accuracy = bal_acc,
      Sensitivity       = sens_val,
      Specificity       = spec_val,
      Youden            = sens_val + spec_val - 1
    )
    
    cat(sprintf("  AUROC: %.3f | Sens: %.3f | Spec: %.3f\n",
                auroc_val, sens_val, spec_val))
  }
  
  # Aggregate results
  valid_folds <- compact(fold_results)
  fold_df <- bind_rows(valid_folds)
  
  if (nrow(fold_df) == 0) {
    cat("  ERROR: No valid fold results\n")
    return(NULL)
  }
  
  summary_stats <- tibble(
    n_folds              = nrow(fold_df),
    threshold_mean       = mean(fold_df$threshold, na.rm = TRUE),
    threshold_sd         = sd(fold_df$threshold, na.rm = TRUE),
    AUROC_mean           = mean(fold_df$AUROC, na.rm = TRUE),
    AUROC_sd             = sd(fold_df$AUROC, na.rm = TRUE),
    PRAUC_mean           = mean(fold_df$PRAUC, na.rm = TRUE),
    PRAUC_sd             = sd(fold_df$PRAUC, na.rm = TRUE),
    Accuracy_mean        = mean(fold_df$Accuracy, na.rm = TRUE),
    Accuracy_sd          = sd(fold_df$Accuracy, na.rm = TRUE),
    Balanced_Accuracy_mean = mean(fold_df$Balanced_Accuracy, na.rm = TRUE),
    Balanced_Accuracy_sd   = sd(fold_df$Balanced_Accuracy, na.rm = TRUE),
    Sensitivity_mean     = mean(fold_df$Sensitivity, na.rm = TRUE),
    Sensitivity_sd       = sd(fold_df$Sensitivity, na.rm = TRUE),
    Specificity_mean     = mean(fold_df$Specificity, na.rm = TRUE),
    Specificity_sd       = sd(fold_df$Specificity, na.rm = TRUE),
    Youden_mean          = mean(fold_df$Youden, na.rm = TRUE),
    Youden_sd            = sd(fold_df$Youden, na.rm = TRUE)
  )
  
  result <- list(
    model_name  = model_name,
    approach    = approach_name,
    microbiome  = microbiome_option,
    method      = "inverse_frequency_weighted",
    fold_results = fold_df,
    summary      = summary_stats
  )
  
  if (approach_name == "Approach3_DataDriven") {
    result$selected_variables <- bind_rows(selected_variables_by_fold)
  }
  
  return(result)
}

# ---------------------------------------------------------------------------
# 2c. Run weighted analysis
# ---------------------------------------------------------------------------

cat("  Running nested CV with inverse-frequency weights...\n\n")

# Model specification
rf_spec <- rand_forest(
  mtry = 4, trees = 500, min_n = 10
) %>%
  set_mode("classification") %>%
  set_engine("ranger", importance = "impurity", num.threads = 1)

result_weighted <- train_with_nested_cv_weighted(
  model_name      = "rf_base",
  model_spec      = rf_spec,
  clinical_data_all   = approach3_clinical_all,
  microbiome_data_all = micro_genus_full,
  approach_name   = "Approach3_DataDriven",
  microbiome_option = "Full_Microbiome",
  cv_folds        = cv_folds,
  weights_df      = sample_weights_df,
  use_pca         = FALSE
)

if (!is.null(result_weighted)) {
  cat("\n  === INVERSE-FREQUENCY WEIGHTING RESULTS ===\n\n")
  cat(sprintf("  AUROC:       %.3f ± %.3f\n",
              result_weighted$summary$AUROC_mean,
              result_weighted$summary$AUROC_sd))
  cat(sprintf("  PRAUC:       %.3f ± %.3f\n",
              result_weighted$summary$PRAUC_mean,
              result_weighted$summary$PRAUC_sd))
}


# =============================================================================
# SECTION 3: COMPARISON TABLE
# =============================================================================

cat("\n--- Section 3: Comparison of approaches ---\n\n")

# Primary analysis values (from decontaminated dataset)
primary_auroc_mean <- 0.813
primary_auroc_sd   <- 0.110
primary_prauc_mean <- 0.592
primary_prauc_sd   <- 0.209

comparison_table <- tibble(
  Analysis = c(
    "Primary analysis (sample-level, unweighted)",
    "Sensitivity: Inverse-frequency weighting"
  ),
  N_observations = c(110, 110),
  N_subjects     = c(43, 43),
  AUROC = c(
    sprintf("%.3f ± %.3f", primary_auroc_mean, primary_auroc_sd),
    if (!is.null(result_weighted))
      sprintf("%.3f ± %.3f",
              result_weighted$summary$AUROC_mean,
              result_weighted$summary$AUROC_sd)
    else "FAILED"
  ),
  PRAUC = c(
    sprintf("%.3f ± %.3f", primary_prauc_mean, primary_prauc_sd),
    if (!is.null(result_weighted))
      sprintf("%.3f ± %.3f",
              result_weighted$summary$PRAUC_mean,
              result_weighted$summary$PRAUC_sd)
    else "FAILED"
  )
)

cat("  ┌───────────────────────────────────────────────────────────────┐\n")
cat("  │  NON-INDEPENDENCE SENSITIVITY ANALYSIS: COMPARISON TABLE     │\n")
cat("  └───────────────────────────────────────────────────────────────┘\n\n")
print(comparison_table, n = Inf, width = Inf)

# Save comparison table
write_csv(comparison_table,
          file.path(output_dir, "nonindependence_comparison.csv"))


# =============================================================================
# SECTION 4: SAVE ALL RESULTS
# =============================================================================

cat("\n--- Section 4: Saving results ---\n\n")

if (!is.null(result_weighted)) {
  write_csv(
    result_weighted$fold_results,
    file.path(output_dir, "weighted_fold_results.csv")
  )
}

# Save sampling structure documentation
write_csv(
  sampling_summary,
  file.path(output_dir, "sampling_structure.csv")
)

cat(sprintf("  All results saved to: %s/\n", output_dir))
cat("\n=== Non-independence sensitivity analysis complete ===\n")

sessionInfo()