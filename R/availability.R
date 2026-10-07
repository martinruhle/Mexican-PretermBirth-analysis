# ============================================================================
# Clinical predictors restricted by availability class
# ----------------------------------------------------------------------------
# config/data_dictionary.csv records, in its `availability` column, when each
# variable is known relative to the visit at which the sample was taken:
#   at_visit        known before or at the visit
#   after_visit     becomes known after the visit (undated diagnoses, birth data)
#   outcome_defined defined by the delivery itself (the outcome and its cut-offs)
# Classification rules and evidence: docs/CLINICAL_VARIABLE_AVAILABILITY.md.
#
# The classes admitted as predictors are a parameter of the analysis
# (config `clinical_availability$allowed`, overridable with PTB_AVAILABILITY),
# and the engine receives them explicitly. Pure functions: everything enters by
# argument except the documented environment-variable override.
# ============================================================================

#' Availability classes of the data dictionary, earliest to latest
#'
#' @return Character vector `c("at_visit", "after_visit", "outcome_defined")`.
#' @export
availability_levels <- function() {
  c("at_visit", "after_visit", "outcome_defined")
}

#' Check a set of admitted availability classes
#'
#' @param allowed Character vector of classes.
#'
#' @return `allowed`, invisibly, if valid. Errors when it is empty, has `NA` or
#'   names an unknown class; warns when it admits `outcome_defined`, whose
#'   variables are defined by the outcome itself.
#' @keywords internal
check_availability_classes <- function(allowed) {
  if (!is.character(allowed) || length(allowed) == 0 || anyNA(allowed) ||
      any(!nzchar(allowed))) {
    stop("allowed_availability must be a non-empty character vector of classes: ",
         paste(availability_levels(), collapse = ", "), call. = FALSE)
  }
  unknown <- setdiff(allowed, availability_levels())
  if (length(unknown) > 0) {
    stop("Unknown availability class(es): ", paste(unknown, collapse = ", "),
         ". Valid classes: ", paste(availability_levels(), collapse = ", "),
         call. = FALSE)
  }
  if ("outcome_defined" %in% allowed) {
    warning("allowed_availability admits 'outcome_defined': these variables are ",
            "defined by the delivery (the outcome itself). Use it only to reproduce ",
            "analyses that predate the availability classification.", call. = FALSE)
  }
  invisible(allowed)
}

#' Availability classes admitted as clinical predictors for this run
#'
#' Reads `cfg$clinical_availability$allowed`. The environment variable
#' `PTB_AVAILABILITY` (comma-separated classes, e.g. `at_visit,after_visit`)
#' overrides it, so the sensitivity analysis runs on the same configuration
#' profile as the main analysis.
#'
#' @param cfg Configuration list from [read_config()].
#' @param override Comma-separated classes; an empty string means "use the
#'   config". Defaults to `Sys.getenv("PTB_AVAILABILITY")`.
#'
#' @return Character vector of admitted classes, in [availability_levels()]
#'   order.
#' @export
resolve_availability <- function(cfg, override = Sys.getenv("PTB_AVAILABILITY", "")) {
  allowed <- if (nzchar(trimws(override))) {
    trimws(strsplit(override, ",", fixed = TRUE)[[1]])
  } else {
    unlist(cfg$clinical_availability$allowed)
  }
  if (is.null(allowed)) {
    stop("The config has no clinical_availability$allowed entry and PTB_AVAILABILITY ",
         "is not set: declare which availability classes are admitted.", call. = FALSE)
  }
  check_availability_classes(allowed)
  intersect(availability_levels(), allowed)
}

#' Variables whose availability class is admitted
#'
#' @param dict Data dictionary (from [read_data_dictionary()]) with columns
#'   `variable` and `availability`.
#' @param allowed Character vector of admitted classes.
#'
#' @return Character vector of dictionary variables whose `availability` is in
#'   `allowed`, in dictionary order.
#' @export
admissible_variables <- function(dict, allowed) {
  check_availability_classes(allowed)
  if (!all(c("variable", "availability") %in% names(dict))) {
    stop("The data dictionary needs the columns `variable` and `availability`.",
         call. = FALSE)
  }
  as.character(dict$variable[dict$availability %in% allowed])
}

#' Restrict a clinical table to the admitted availability classes
#'
#' Every non-key column of `clinical_data` is looked up in the dictionary. A
#' column with no dictionary row is an error in both modes: its availability is
#' unknown, so it can neither be admitted nor removed silently. Columns of a
#' class not in `allowed` are removed (`action = "drop"`) or make the call fail
#' (`action = "error"`). Use `"drop"` for a pool from which variables are still
#' to be selected, and `"error"` for a variable list that is already final:
#' removing a variable from a final list would bypass the selection that built it.
#'
#' @param clinical_data Data frame with key columns and clinical variables.
#' @param dict Data dictionary with columns `variable` and `availability`.
#' @param allowed Character vector of admitted classes.
#' @param keys Key columns that are not predictors (identifiers and outcome).
#' @param action `"drop"` or `"error"`.
#'
#' @return A list with `data` (`clinical_data` without the excluded columns) and
#'   `excluded` (data frame with the excluded `variable` and its `availability`).
#' @export
restrict_to_availability <- function(clinical_data, dict, allowed,
                                     keys = c("index", "id", "preterm"),
                                     action = c("drop", "error")) {
  action <- match.arg(action)
  admitted <- admissible_variables(dict, allowed)

  vars <- setdiff(names(clinical_data), keys)
  unknown <- setdiff(vars, dict$variable)
  if (length(unknown) > 0) {
    stop("Clinical variable(s) without a row in the data dictionary (availability ",
         "unknown): ", paste(unknown, collapse = ", "),
         ". Add them to config/data_dictionary.csv with an availability class.",
         call. = FALSE)
  }

  out <- setdiff(vars, admitted)
  excluded <- data.frame(
    variable     = out,
    availability = as.character(dict$availability[match(out, dict$variable)]),
    stringsAsFactors = FALSE)

  if (length(out) > 0 && action == "error") {
    stop("Clinical variable(s) of a class not admitted (",
         paste(allowed, collapse = ", "), "): ",
         paste(sprintf("%s [%s]", excluded$variable, excluded$availability),
               collapse = ", "),
         ". Restrict the candidate list before selecting the variables.",
         call. = FALSE)
  }

  list(data = clinical_data[, setdiff(names(clinical_data), out), drop = FALSE],
       excluded = excluded)
}
