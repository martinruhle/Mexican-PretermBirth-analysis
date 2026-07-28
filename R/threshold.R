# ============================================================================
# Classification-threshold optimisation (Youden / F1)
# ----------------------------------------------------------------------------
# Extracted verbatim from the chunk `nested_cv_helper_functions` of
# analysis/integrated_preterm_prediction_workflow.Rmd (Chat 2). Logic unchanged;
# only roxygen docs were added.
#
# Leakage invariant (see nested_cv.R): this is called on the INNER-VALIDATION
# predictions to pick the threshold, which is then applied to the independent
# OUTER-TEST fold. It must never see outer-test data.
#
# Depends on the analysis packages being attached (pROC::roc, pROC::coords, and
# dplyr for %>%/mutate/filter/slice), exactly as in the source .Rmd.
# ============================================================================

#' Optimise a classification threshold on validation predictions
#'
#' Builds the ROC curve on a set of validation predictions and returns the
#' probability threshold that maximises the chosen criterion (Youden's J by
#' default, or F1).
#'
#' @param val_predictions A tibble/data.frame with columns `true_class` (factor
#'   with levels `c("0", "1")`) and `pred_prob` (predicted probability of the
#'   positive class).
#' @param method Optimisation criterion: `"youden"` (maximise
#'   sensitivity + specificity - 1) or `"f1"` (maximise F1).
#'
#' @return A list with elements `threshold`, `sensitivity`, `specificity` and
#'   `criterion_value` at the selected operating point.
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

  # ROC curve on validation set
  roc_obj <- roc(val_predictions$true_class,
                 val_predictions$pred_prob,
                 levels = c("0", "1"),
                 direction = "auto",
                 quiet = TRUE)

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
