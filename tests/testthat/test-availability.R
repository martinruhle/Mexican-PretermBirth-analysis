# Clinical predictors are restricted by the `availability` class of the data
# dictionary (R/availability.R). The admitted classes are a parameter of the
# engine: the main analysis admits `at_visit`, the sensitivity analysis adds
# `after_visit`, and `outcome_defined` enters neither. These tests fail if a
# variable of a class that is not admitted can reach a model: through the
# Approach 3 screening pool, or through a fixed Approach 1/2 list.

main_classes        <- "at_visit"
sensitivity_classes <- c("at_visit", "after_visit")

# Variables of the Approach 3 pool that are not known at the visit
# (docs/CLINICAL_VARIABLE_AVAILABILITY.md, section 4).
late_after   <- c("comppreeclam", "compvaginf", "diabetes_gest", "oligohidramnios", "rciu",
                  "sex_baby")
late_outcome <- c("rpm", "rpm_preterm", "complication_count", "any_complication",
                  "multiple_complications")

test_that("admissible_variables returns exactly the dictionary rows of the admitted classes", {
  dict <- read_data_dictionary(read_config("example"))
  avail <- stats::setNames(dict$availability, dict$variable)

  main <- admissible_variables(dict, main_classes)
  expect_true(all(avail[main] == "at_visit"))
  expect_setequal(main, dict$variable[dict$availability == "at_visit"])
  expect_true(all(c("sdg_visita", "edad_cronologicamujer") %in% main))
  expect_length(intersect(main, c(late_after, late_outcome)), 0)

  sens <- admissible_variables(dict, sensitivity_classes)
  expect_true(all(late_after %in% sens))
  expect_length(intersect(sens, late_outcome), 0)
})

test_that("resolve_availability reads the config, accepts an override and rejects unknown classes", {
  cfg <- read_config("example")
  expect_identical(resolve_availability(cfg, override = ""), main_classes)
  # order follows availability_levels(), whatever the order of the override
  expect_identical(resolve_availability(cfg, override = "after_visit, at_visit"),
                   sensitivity_classes)
  expect_error(resolve_availability(cfg, override = "at_visit,visit_day"), "visit_day")
  expect_warning(resolve_availability(cfg, override = "at_visit,outcome_defined"),
                 "outcome_defined")
  cfg$clinical_availability <- NULL
  expect_error(resolve_availability(cfg, override = ""), "clinical_availability")
})

test_that("restrict_to_availability drops a pool, rejects a final list, and refuses undeclared variables", {
  dict <- read_data_dictionary(read_config("example"))
  clin <- data.frame(index = paste0("S", 1:4), id = paste0("P", 1:4), preterm = c(0, 1, 0, 1),
                     sdg_visita = c(20, 22, 25, 30), comppreeclam = c(0, 1, 0, 0),
                     rpm = c(1, 0, 0, 0), stringsAsFactors = FALSE)

  pool <- restrict_to_availability(clin, dict, main_classes, action = "drop")
  expect_identical(names(pool$data), c("index", "id", "preterm", "sdg_visita"))
  expect_identical(pool$excluded$variable, c("comppreeclam", "rpm"))
  expect_identical(pool$excluded$availability, c("after_visit", "outcome_defined"))

  pool_sens <- restrict_to_availability(clin, dict, sensitivity_classes, action = "drop")
  expect_identical(names(pool_sens$data), c("index", "id", "preterm", "sdg_visita", "comppreeclam"))

  expect_error(restrict_to_availability(clin, dict, main_classes, action = "error"),
               "rpm \\[outcome_defined\\]")

  clin$undeclared_score <- 1:4
  expect_error(restrict_to_availability(clin, dict, sensitivity_classes, action = "drop"),
               "undeclared_score")
})

# ---- the engine ------------------------------------------------------------
# Approach 3: the pool the univariate screening receives. The screening helper
# (calculate_completeness, defined in the analysis .Rmd and resolved by the engine
# from the calling environment) is replaced by a spy that records its input and
# stops the run, so the test checks the pool without fitting anything.
spy_on_screening <- function(code) {
  seen <- new.env()
  had <- exists("calculate_completeness", envir = globalenv(), inherits = FALSE)
  old <- if (had) get("calculate_completeness", envir = globalenv()) else NULL
  on.exit({
    if (had) assign("calculate_completeness", old, envir = globalenv())
    else rm("calculate_completeness", envir = globalenv())
  }, add = TRUE)
  assign("calculate_completeness", function(data, vars) {
    seen$vars <- vars
    seen$cols <- names(data)
    stop(structure(class = c("screening_reached", "error", "condition"),
                   list(message = "screening reached", call = NULL)))
  }, envir = globalenv())
  invisible(utils::capture.output(
    tryCatch(force(code), screening_reached = function(e) seen$reached <- TRUE)))
  seen
}

engine_inputs <- function() {
  set.seed(1)
  n <- 12
  clin <- data.frame(
    index = paste0("S", seq_len(n)), id = paste0("P", seq_len(n)),
    preterm = rep(c(0, 0, 1), length.out = n),
    sdg_visita = runif(n, 12, 30), edad_cronologicamujer = runif(n, 18, 40),
    hemoglobin_g_dl = runif(n, 10, 14),
    comppreeclam = rbinom(n, 1, 0.2), oligohidramnios = rbinom(n, 1, 0.2),
    sex_baby = sample(c("F", "M"), n, replace = TRUE),
    rpm = rbinom(n, 1, 0.2), rpm_preterm = rbinom(n, 1, 0.2),
    stringsAsFactors = FALSE)
  clin$complication_count <- clin$comppreeclam + clin$oligohidramnios + clin$rpm + clin$rpm_preterm
  labels <- data.frame(id = clin$id, preterm = as.character(clin$preterm), stringsAsFactors = FALSE)
  list(clin = clin, folds = rsample::vfold_cv(labels, v = 3),
       dict = read_data_dictionary(read_config("example")))
}

run_engine <- function(inp, approach, allowed, clin = inp$clin) {
  train_with_nested_cv(
    model_name = "glmnet_base", model_spec = NULL,
    clinical_data_all = clin, microbiome_data_all = NULL,
    approach_name = approach, microbiome_option = "Full_Microbiome",
    cv_folds = inp$folds, clr_zero_levels = NULL,
    dict = inp$dict, allowed_availability = allowed)
}

test_that("the Approach 3 screening pool holds only admitted classes", {
  inp <- engine_inputs()
  pool_all <- setdiff(names(inp$clin), c("index", "id", "preterm"))

  for (allowed in list(main_classes, sensitivity_classes)) {
    seen <- spy_on_screening(run_engine(inp, "Approach3_DataDriven", allowed))
    expect_true(isTRUE(seen$reached))
    expected <- intersect(pool_all, admissible_variables(inp$dict, allowed))
    expect_identical(seen$vars, expected)
    expect_identical(setdiff(seen$cols, c("index", "id", "preterm")), expected)
  }
})

test_that("a fixed Approach 1/2 list with a class that is not admitted stops the engine", {
  inp <- engine_inputs()
  a2 <- inp$clin[, c("index", "id", "sdg_visita", "edad_cronologicamujer", "comppreeclam",
                     "rpm", "preterm")]
  expect_error(invisible(utils::capture.output(
    run_engine(inp, "Approach2_Literature", main_classes, clin = a2))),
    "comppreeclam \\[after_visit\\], rpm \\[outcome_defined\\]")
  expect_error(invisible(utils::capture.output(
    run_engine(inp, "Approach1_DREAM", sensitivity_classes, clin = a2))),
    "rpm \\[outcome_defined\\]")
})
