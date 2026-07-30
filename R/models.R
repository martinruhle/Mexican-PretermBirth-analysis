# ============================================================================
# Model specifications — config-driven parsnip specs (Chat 4)
# ----------------------------------------------------------------------------
# Extrae el chunk `model_specifications` del .Rmd a una funcion pura. Construye
# las especificaciones parsnip UNICAMENTE de los modelos activos declarados en
# config$models (rf_base, glmnet_base -> 12 combinaciones = 2 modelos x 3
# approaches x 2 microbioma). Los specs de config$models_future (xgb_base,
# tree_base, rf_pca) se IGNORAN por diseno: reservados para un dataset ampliado
# (ver roadmap en CLAUDE.md); reactivarlos = moverlos a config$models +
# revalidar convergencia + nuevo baseline.
#
# Refactor puro: los specs producidos son identicos a los que estaban
# hardcodeados en el .Rmd (mismos hiperparametros, mismo engine, mismo mode).
#
# Depends on the analysis packages being attached (parsnip via tidymodels, and
# magrittr %>% via tidyverse), exactly as in the source .Rmd.
# ============================================================================

#' Build parsnip model specifications from config
#'
#' Constructs the fixed-hyperparameter parsnip specs for the ACTIVE models only,
#' i.e. those declared under `cfg$models`. Dispatches on each model's `engine`
#' field. `cfg$models_future` (reserved specs) is intentionally ignored, so the
#' active set — and hence the 12 nested-CV combinations — is unchanged.
#'
#' @param cfg Configuration list from [read_config()]. Must contain `cfg$models`,
#'   a named list where each element carries an `engine` field plus that engine's
#'   hyperparameters. Supported engines:
#'   * `"ranger"`  — reads `mtry`, `trees`, `min_n`, `importance`, `num_threads`.
#'   * `"glmnet"`  — reads `penalty`, `mixture`.
#'
#' @return A named list of parsnip model specifications, one per `cfg$models`
#'   entry, preserving the order in which they are declared in the config.
#'
#' @export
build_model_specs <- function(cfg) {
  models <- cfg$models
  if (is.null(models) || length(models) == 0) {
    stop("build_model_specs(): cfg$models esta vacio o ausente.", call. = FALSE)
  }

  specs <- lapply(names(models), function(nm) {
    m <- models[[nm]]
    engine <- m$engine

    if (identical(engine, "ranger")) {
      rand_forest(mtry = m$mtry, trees = m$trees, min_n = m$min_n) %>%
        set_mode("classification") %>%
        set_engine("ranger", importance = m$importance, num.threads = m$num_threads)

    } else if (identical(engine, "glmnet")) {
      logistic_reg(penalty = m$penalty, mixture = m$mixture) %>%
        set_engine("glmnet") %>%
        set_mode("classification")

    } else {
      eng_disp <- if (is.null(engine)) "NULL" else engine
      stop(sprintf(
        paste0("build_model_specs(): engine '%s' no soportado (modelo '%s'). ",
               "Los specs de models_future (xgboost/rpart) requieren reactivacion ",
               "explicita — ver roadmap en CLAUDE.md."),
        eng_disp, nm), call. = FALSE)
    }
  })

  names(specs) <- names(models)
  specs
}
