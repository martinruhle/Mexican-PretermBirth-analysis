# build_model_specs(): construct parsnip specs for the ACTIVE models only
# (cfg$models), dispatching on each model's engine. cfg$models_future is reserved
# and must be ignored, keeping the active set at rf_base + glmnet_base.

test_that("build_model_specs builds exactly the active models and ignores models_future", {
  specs <- build_model_specs(make_cfg())

  expect_named(specs, c("rf_base", "glmnet_base"))     # models_future not included
  expect_s3_class(specs$rf_base, "model_spec")
  expect_s3_class(specs$glmnet_base, "model_spec")
  expect_identical(specs$rf_base$engine, "ranger")
  expect_identical(specs$glmnet_base$engine, "glmnet")
  expect_identical(specs$rf_base$mode, "classification")
  expect_identical(specs$glmnet_base$mode, "classification")
})

test_that("build_model_specs preserves the declared order of cfg$models", {
  cfg <- make_cfg()
  cfg$models <- cfg$models[c("glmnet_base", "rf_base")]   # reverse the declaration order
  expect_named(build_model_specs(cfg), c("glmnet_base", "rf_base"))
})

test_that("build_model_specs errors clearly on an unsupported engine", {
  cfg <- make_cfg()
  cfg$models <- list(weird = list(engine = "svm"))
  expect_error(build_model_specs(cfg), "no soportado")
})

test_that("build_model_specs errors when cfg$models is empty or absent", {
  expect_error(build_model_specs(list(models = list())), "esta vacio o ausente")
  expect_error(build_model_specs(list()),                "esta vacio o ausente")
})
