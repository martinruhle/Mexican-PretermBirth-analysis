#!/usr/bin/env Rscript
# =============================================================================
# run_baseline.R  —  Captura el BASELINE post-fix (antes de refactorizar en Chat 2)
# -----------------------------------------------------------------------------
# POR QUE: el Chat 2 extrae funciones del .Rmd a R/ SIN cambiar resultados. Para
# poder probarlo, hay que registrar las metricas AUROC/PRAUC de las 12 combinaciones
# ANTES de tocar nada, ya con el CLR corregido (per-sample + cmultRepl) y el seed
# unificado a 123. Ese numero es la "verdad" contra la que se compara tras cada fase.
#
# QUE HACE: renderiza el .Rmd tal cual y guarda el HTML en analysis/_output/. El HTML
# incluye la tabla de las 12 combinaciones (chunk `cv_results_table`). Tarda ~2.5-3 h
# (corre el nested CV completo + ANCOM-BC2 por fold).
#
# COMO CORRERLO (renv aun NO esta restaurado; hay que saltarlo para usar la libreria
# base de R, donde SI estan tidymodels/ranger/glmnet/ANCOMBC/zCompositions/...):
#
#   RENV_CONFIG_AUTOLOADER_ENABLED=FALSE R_PROFILE_USER=/dev/null \
#     "/c/Program Files/R/R-4.4.2/bin/Rscript.exe" scripts/run_baseline.R
#
# REQUISITO: los datos reales en C:/Users/marti/Documents/Datos_mexicanos/ (el .Rmd
# los lee de ahi via setwd(); eso se elimina en el Chat 3).
# =============================================================================

stopifnot(requireNamespace("rmarkdown", quietly = TRUE))

rmd    <- normalizePath("analysis/integrated_preterm_prediction_workflow.Rmd", mustWork = TRUE)
outdir <- "analysis/_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message(">> Renderizando el pipeline completo... (nested CV + ANCOM-BC2; tarda horas)")
t0 <- Sys.time()
rmarkdown::render(
  input       = rmd,
  output_dir  = normalizePath(outdir),
  output_file = "baseline_postfix.html",
  envir       = new.env()
)
message(sprintf(">> Listo en %.1f min. Baseline: %s",
                as.numeric(difftime(Sys.time(), t0, units = "mins")),
                file.path(outdir, "baseline_postfix.html")))
message(">> Guarda la tabla de 12 combinaciones (chunk cv_results_table) como referencia del Chat 2.")
