# low_education and ses_risk_score are defined as "secondary education or lower"
# (config/data_dictionary.csv). In the cohort table nivel_academico is text with
# the code first ("2.- Secondary education"), while the synthetic example stores
# the bare integer. These tests use the cohort's text values, so they fail if the
# level is compared as a string: "2.- Secondary education" <= 2 is FALSE in R,
# and only primary education would be flagged.

cohort_levels <- c("1.- Primary education or less", "2.- Secondary education",
                   "3.- Vocational/Technical", "4.- Associates degree or higher")

ses_fixture <- function(education) {
  data.frame(id = seq_along(education), nivel_academico = education,
             maritalstat = c(1, 1, 1, 1, 5, 2, 1)[seq_along(education)],
             stringsAsFactors = FALSE)
}

test_that("low_education flags primary and secondary education given as text", {
  out <- derive_ses_features(ses_fixture(c(cohort_levels, cohort_levels[2], cohort_levels[3], NA)))
  expect_identical(out$low_education, c(1, 1, 0, 0, 1, 0, NA))
  expect_identical(out$unmarried,      c(0, 0, 0, 0, 1, 1, 0))
  expect_identical(out$ses_risk_score, c(1, 1, 0, 0, 2, 1, NA))
})

test_that("text and integer codes give the same indicators", {
  as_text <- derive_ses_features(ses_fixture(cohort_levels))
  as_code <- derive_ses_features(ses_fixture(1:4))
  cols <- c("low_education", "unmarried", "ses_risk_score")
  expect_identical(as_text[cols], as_code[cols])
})

test_that("the three indicators are appended after the input columns", {
  out <- derive_ses_features(ses_fixture(cohort_levels))
  expect_identical(names(out), c("id", "nivel_academico", "maritalstat",
                                 "low_education", "unmarried", "ses_risk_score"))
})

test_that("an education value without a leading code is an error", {
  expect_error(derive_ses_features(ses_fixture(c(cohort_levels[1:3], "Secondary"))),
               "nivel_academico")
})
