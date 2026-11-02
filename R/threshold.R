# ============================================================================
# Classification-threshold optimisation (Youden / F1)
# ----------------------------------------------------------------------------
# Extracted from the chunk `nested_cv_helper_functions` of
# analysis/integrated_preterm_prediction_workflow.Rmd (Chat 2). Edits relative to
# the source:
#   * roxygen docs added.
#   * threshold search (2026-09): the validation ROC curve is built with
#     roc_ptb() (nested_cv.R; pROC direction "<", preterm = higher score) instead
#     of direction = "auto", so the operating points it scores are those of the
#     rule the engine applies to the outer-test fold. Rationale in the roxygen
#     block below.
#
# Leakage invariant (see nested_cv.R): this is called on the INNER-VALIDATION
# predictions to pick the threshold, which is then applied to the independent
# OUTER-TEST fold. It must never see outer-test data.
#
# Depends on roc_ptb() (nested_cv.R) and on the analysis packages being attached
# (pROC::coords, and dplyr for %>%/mutate/filter/slice), exactly as in the source
# .Rmd.
# ============================================================================

#' Optimise a classification threshold on validation predictions
#'
#' Builds the ROC curve on a set of validation predictions and returns the
#' probability threshold that maximises the chosen criterion (Youden's J by
#' default, or F1).
#'
#' The curve is built with [roc_ptb()], i.e. with the orientation fixed a
#' priori (`direction = "<"`): at every candidate threshold `t`, the
#' sensitivity and specificity are those of the rule
#' `pred_prob >= t -> "1"` (preterm), which is the rule
#' [train_with_nested_cv()] applies to the outer-test fold. The selected
#' operating point therefore describes the classifier that is actually used,
#' and it shares its orientation with the reported AUROC.
#'
#' Why not `direction = "auto"`: `"auto"` orients the curve by comparing the
#' class medians of the validation predictions. When preterm subjects scored
#' lower than term subjects there, pROC oriented the curve the other way
#' (`">"`: positive when `pred_prob <= t`), so the maximum of J belonged to the
#' opposite rule. Because pROC thresholds lie between observed scores, the two
#' rules are exact complements on the validation set, and the selected
#' threshold was the operating point with the LOWEST J for the rule applied
#' afterwards. With the fixed orientation the threshold is always the best
#' operating point of the rule that is applied; when the model ranks preterm
#' subjects lower on the validation set, that best J can be close to or below
#' 0, which is the faithful outcome for a model that does not discriminate
#' there.
#'
#' @param val_predictions A tibble/data.frame with columns `true_class` (factor
#'   with levels `c("0", "1")`) and `pred_prob` (predicted probability of the
#'   positive class).
#' @param method Optimisation criterion: `"youden"` (maximise
#'   sensitivity + specificity - 1) or `"f1"` (maximise F1).
#'
#' @return A list with elements `threshold`, `sensitivity`, `specificity` and
#'   `criterion_value` at the selected operating point, all for the rule
#'   `pred_prob >= threshold -> "1"`.
#'
#' @export
optimize_threshold_cv <- function(val_predictions, method = "youden") {
  # Optimizes classification threshold on validation set
  #
  # Args:
  #   val_predictions: tibble with columns 'true_class' and 'pred_prob'
  #   method: optimization criterion ('youden' or 'f1')
  #
  # Returns:
  #   list with optimal threshold and associated metrics

  # ROC curve on validation set, with the same fixed orientation as the reported
  # AUROC and as the outer-test rule (pred_prob >= threshold -> "1").
  roc_obj <- roc_ptb(val_predictions$true_class, val_predictions$pred_prob)

  # Get all possible thresholds with their metrics
  coords_all <- coords(roc_obj, "all",
                       ret = c("threshold", "sensitivity", "specificity"))

  if(method == "youden") {
    # Maximize Youden Index (J = Sensitivity + Specificity - 1)
    coords_all <- coords_all %>%
      mutate(criterion = sensitivity + specificity - 1)
  } else if(method == "f1") {
    # Alternative: Maximize F1 score
    coords_all <- coords_all %>%
      mutate(
        precision = sensitivity / (sensitivity + (1 - specificity) + 1e-10),
        criterion = 2 * (precision * sensitivity) / (precision + sensitivity + 1e-10)
      )
  }

  # Select optimal threshold
  optimal <- coords_all %>%
    filter(!is.na(criterion), !is.infinite(threshold)) %>%
    filter(criterion == max(criterion, na.rm = TRUE)) %>%
    slice(1)

  return(list(
    threshold = optimal$threshold,
    sensitivity = optimal$sensitivity,
    specificity = optimal$specificity,
    criterion_value = optimal$criterion
  ))
}
