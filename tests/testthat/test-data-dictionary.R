# config/data_dictionary.csv records, for every variable, WHEN its value is known
# relative to the sample collection visit (`availability`), the evidence for that
# classification (`availability_source`) and, for derived variables, their inputs
# (`derived_from`). A model may only use variables known at the visit; this file is
# where that is decided. Classification rules: docs/CLINICAL_VARIABLE_AVAILABILITY.md.

availability_levels <- c("at_visit", "after_visit", "outcome_defined")   # earliest -> latest

test_that("every dictionary row has an availability class and a source", {
  dict <- read_data_dictionary(read_config("example"))
  expect_true(all(c("availability", "availability_source", "derived_from") %in% names(dict)))
  expect_false(anyNA(dict$availability))
  expect_true(all(dict$availability %in% availability_levels))
  expect_false(anyNA(dict$availability_source))
  expect_true(all(nzchar(trimws(dict$availability_source))))
})

test_that("the outcome, gestational age at delivery and rupture-of-membranes flags are outcome_defined", {
  cfg  <- read_config("example")
  dict <- read_data_dictionary(cfg)
  avail <- stats::setNames(dict$availability, dict$variable)
  expect_identical(unname(avail[cfg$columns$outcome]), "outcome_defined")
  expect_identical(unname(avail[cfg$columns$gestational_age]), "outcome_defined")
  # rpm = rupture after week 37 (a 1 implies a term delivery); rpm_preterm = rupture
  # before week 37 (same cut-off as the outcome).
  expect_identical(unname(avail["rpm"]), "outcome_defined")
  expect_identical(unname(avail["rpm_preterm"]), "outcome_defined")
})

test_that("every clinical variable named in the approach definitions has a dictionary row", {
  cfg  <- read_config("example")
  dict <- read_data_dictionary(cfg)
  a2   <- cfg$approaches$approach2_literature
  used <- c(unlist(cfg$approaches$approach1_dream), names(a2$core_vars), names(a2$candidates))
  expect_true(length(used) > 0)
  expect_identical(setdiff(used, dict$variable), character(0))
})

test_that("derived variables list existing inputs and inherit the latest availability among them", {
  dict <- read_data_dictionary(read_config("example"))
  avail <- stats::setNames(dict$availability, dict$variable)

  # pipeline-derived rows must declare their inputs
  expect_false(anyNA(dict$derived_from[dict$role == "clinical_derived"]))

  derived <- dict[!is.na(dict$derived_from) & nzchar(dict$derived_from), ]
  expect_true(nrow(derived) > 0)
  inputs <- strsplit(derived$derived_from, ";", fixed = TRUE)
  expect_identical(setdiff(unlist(inputs), dict$variable), character(0))

  latest <- vapply(inputs, function(x) availability_levels[max(match(avail[x], availability_levels))],
                   character(1))
  # names of derived variables whose class differs from the latest of their inputs
  expect_identical(derived$variable[derived$availability != latest], character(0))
})
