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

# rmarkdown necesita pandoc; el Rscript "pelado" no lo ve (RStudio lo trae bundled y
# solo lo expone dentro de la IDE). Si no esta disponible, localizamos el pandoc de
# RStudio/Quarto y apuntamos RSTUDIO_PANDOC para que el render headless funcione igual.
if (!rmarkdown::pandoc_available()) {
  cand <- c(
    Sys.getenv("RSTUDIO_PANDOC"),
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
    "C:/Program Files/RStudio/bin/quarto/bin/tools",
    "C:/Program Files/RStudio/bin/pandoc",
    "C:/Program Files/Quarto/bin/tools"
  )
  cand <- cand[nzchar(cand)]
  exe  <- if (.Platform$OS.type == "windows") "pandoc.exe" else "pandoc"
  hit  <- cand[file.exists(file.path(cand, exe))]
  if (length(hit)) Sys.setenv(RSTUDIO_PANDOC = hit[1])
}
stopifnot("pandoc no encontrado (instala RStudio/Quarto o define RSTUDIO_PANDOC)" =
            rmarkdown::pandoc_available())

# El MISMO script sirve para el baseline pre-refactor (Paso 0) y la verificacion
# post-refactor (Paso final): basta cambiar el tag. Los dos *_metrics.csv se comparan
# fila a fila. Tag por defecto = baseline_postfix.
tag    <- Sys.getenv("PTB_RUN_TAG", "baseline_postfix")
rmd    <- normalizePath("analysis/integrated_preterm_prediction_workflow.Rmd", mustWork = TRUE)
outdir <- "analysis/_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Conservamos el entorno de render (no new.env() anonimo) para poder extraer despues
# la tabla de las 12 combinaciones (cv_results_all) y guardarla como referencia.
render_env <- new.env()

message(sprintf(">> [%s] Renderizando el pipeline completo... (nested CV + ANCOM-BC2; tarda horas)", tag))
t0 <- Sys.time()
# render() puede fallar en un chunk de FIGURA posterior (p.ej. el bug pre-existente
# `top25_importance` en figure5_panel_a). Las 12 metricas se calculan mucho antes
# (chunks run_nested_cv y cv_results_table), y knitr las deja en render_env aunque un
# chunk mas abajo reviente. Por eso envolvemos en tryCatch y SEGUIMOS a la captura:
# la comparacion del Chat 2 es sobre las 12 metricas, no sobre el HTML completo.
render_ok <- tryCatch({
  rmarkdown::render(
    input       = rmd,
    output_dir  = normalizePath(outdir),
    output_file = paste0(tag, ".html"),
    envir       = render_env
  )
  TRUE
}, error = function(e) {
  message(sprintf(">> AVISO: render() no completo el HTML (fallo un chunk posterior): %s",
                  conditionMessage(e)))
  message(">> Continuo: las 12 metricas se computan antes del fallo; las capturo desde render_env.")
  FALSE
})

# -- Persistir las 12 metricas en forma NUMERICA (no el "mean +/- sd" de texto del
#    chunk cv_results_table) para poder diffear exacto tras cada fase. No altera nada
#    del computo: solo lee cv_results_all del entorno de render y lo vuelca a disco. ----
if (exists("cv_results_all", envir = render_env)) {
  cv_results_all <- get("cv_results_all", envir = render_env)
  metrics <- do.call(rbind, lapply(cv_results_all, function(res) {
    if (is.null(res)) return(NULL)
    s <- res$summary
    data.frame(
      Model = res$model_name, Approach = res$approach, Microbiome = res$microbiome,
      AUROC_mean = s$AUROC_mean, AUROC_sd = s$AUROC_sd,
      PRAUC_mean = s$PRAUC_mean, PRAUC_sd = s$PRAUC_sd,
      Sensitivity_mean = s$Sensitivity_mean, Sensitivity_sd = s$Sensitivity_sd,
      Specificity_mean = s$Specificity_mean, Specificity_sd = s$Specificity_sd,
      Accuracy_mean = s$Accuracy_mean, Accuracy_sd = s$Accuracy_sd,
      Balanced_Accuracy_mean = s$Balanced_Accuracy_mean, Balanced_Accuracy_sd = s$Balanced_Accuracy_sd,
      threshold_mean = s$threshold_mean, threshold_sd = s$threshold_sd,
      Youden_mean = s$Youden_mean, Youden_sd = s$Youden_sd,
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(metrics) || nrow(metrics) == 0) {
    stop("cv_results_all existe pero no tiene combinaciones completas: el nested CV fallo ",
         "antes de producir resultados (revisa el error del chunk run_nested_cv arriba).")
  }
  metrics <- metrics[order(-metrics$AUROC_mean, -metrics$PRAUC_mean), ]
  csv_path <- file.path(outdir, paste0(tag, "_metrics.csv"))
  utils::write.csv(metrics, csv_path, row.names = FALSE)
  saveRDS(lapply(cv_results_all, function(r) if (is.null(r)) NULL else r$summary),
          file.path(outdir, paste0(tag, "_cv_summary.rds")))
  message(sprintf(">> Metricas de %d combinaciones guardadas en %s", nrow(metrics), csv_path))
} else {
  warning("cv_results_all no existe en el entorno de render; no se guardo la tabla de metricas.")
}

message(sprintf(">> Listo en %.1f min. HTML %s: %s",
                as.numeric(difftime(Sys.time(), t0, units = "mins")),
                if (isTRUE(render_ok)) "completo" else "INCOMPLETO (fallo un chunk de figura)",
                file.path(outdir, paste0(tag, ".html"))))
