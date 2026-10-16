# ============================================================================
# Socioeconomic indicators derived from the clinical table
# ----------------------------------------------------------------------------
# The chunk `feature_engineering` of the analysis .Rmd derives three
# socioeconomic indicators from educational attainment and marital status:
#   low_education   1 if the highest level is secondary education or lower
#   unmarried       1 if the marital status is other than married
#   ses_risk_score  low_education + unmarried (0 = lower risk, 2 = higher risk)
# Their dictionary rows (config/data_dictionary.csv, role clinical_derived)
# record the definitions and the availability class they inherit.
# Pure functions: everything enters by argument.
# ============================================================================

#' Numeric code of an educational-attainment value
#'
#' The cohort table stores the level as text with the code first
#' (`"2.- Secondary education"`); the synthetic example stores the bare integer.
#' Comparing the text with a number compares strings (`"2.- Secondary
#' education" <= 2` is `FALSE`), so the code is read before any comparison.
#'
#' @param x Character or numeric vector of education levels.
#' @param column Column name, used in the error message.
#'
#' @return Integer vector of codes (`NA` where `x` is missing). Errors when a
#'   non-missing text value does not start with a numeric code.
#' @keywords internal
education_level_code <- function(x, column = "education") {
  if (is.numeric(x)) return(as.integer(x))
  x <- trimws(as.character(x))
  has_code <- grepl("^[0-9]+", x)
  bad <- !is.na(x) & !has_code
  if (any(bad)) {
    stop("Column '", column, "' has values without a leading numeric code: ",
         paste(unique(x[bad]), collapse = ", "), call. = FALSE)
  }
  code <- rep(NA_integer_, length(x))
  code[has_code] <- as.integer(regmatches(x[has_code], regexpr("^[0-9]+", x[has_code])))
  code
}

#' Derive the socioeconomic indicators
#'
#' Appends `low_education`, `unmarried` and `ses_risk_score` (in that order) to
#' `data`.
#'
#' @param data Data frame with the education and marital-status columns.
#' @param education Name of the educational-attainment column (codes 1 = primary
#'   or less, 2 = secondary, 3 = vocational/technical, 4 = associate degree or
#'   higher), as integers or as text that starts with the code.
#' @param marital Name of the marital-status column (1 = married).
#'
#' @return `data` with the three indicator columns appended. Missing inputs give
#'   missing indicators (imputed later inside the model recipe).
#' @export
derive_ses_features <- function(data, education = "nivel_academico",
                                marital = "maritalstat") {
  education_code      <- education_level_code(data[[education]], column = education)
  data$low_education  <- ifelse(education_code <= 2, 1, 0)
  data$unmarried      <- ifelse(data[[marital]] != 1, 1, 0)
  data$ses_risk_score <- data$low_education + data$unmarried
  data
}
