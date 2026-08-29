#!/usr/bin/env Rscript
# =============================================================================
# Permutation test for model discrimination — DRIVER (Fase 0/1/2)
# -----------------------------------------------------------------------------
# Reemplaza a analysis/permutation_test_nested_cv.Rmd, que es un .Rmd sin YAML,
# no corre standalone y usa la firma vieja del motor (le faltan clr_zero_levels /
# genera_clean / abs_data, obligatorios desde Chat 2/3).
#
# MODELO OBJETIVO (Decision C2): glmnet_base | Approach3_DataDriven | ANCOM_Taxa
#   Es el mejor modelo del baseline actual (AUROC 0.760) — el resultado que se
#   reporta y que el revisor pidio validar. El script viejo apuntaba a
#   rf_base|Approach3|Full_Microbiome (hoy 0.680, 6.º).
#
# ESTADISTICO (Fase 2, Decision A1 + B2):
#   * AUROC_auto  = el que reporta el motor (pROC direction="auto").
#     "auto" elige por fold la orientacion que da AUC>=0.5 -> bajo etiquetas
#     permutadas infla la nula a ~0.73 (medido) en vez de 0.5. NO sirve de nula.
#   * AUROC_fixed = recalculado desde las predicciones out-of-fold que ahora
#     exporta el motor (result$test_predictions), con direction="<" FIJA
#     (controles<casos, la orientacion semantica de P(preterm=1)).
#     Mismo metodo EXACTO para el observado y para la nula.
#
# POLITICA DE FALLOS (Decision (a)): se promedian solo los folds computables.
#   Los folds se estratificaron con las etiquetas VERDADERAS; al permutar, un
#   fold de test puede quedar con una sola clase y su AUROC no existe. El motor
#   ahora saltea ese fold (guard en R/nested_cv.R) en vez de abortar la corrida.
#   Se reporta n_folds_valid por permutacion.
#
# ⚠ NO-FUGA — las TRES columnas de etiquetas se permutan:
#   1. clinical_data_all$preterm   (outcome del modelo)
#   2. subject_labels  (global)    (estratificacion del inner split 70/30)
#   3. abs_data$meta$preterm       (seleccion de taxa por ANCOM-BC2)
#   Si (3) no se permuta, ANCOM-BC2 elige taxa con las etiquetas verdaderas y la
#   nula queda contaminada. Referencia: tests/testthat/test-leakage-permutation.R
#   cv_folds NO se toca: el motor solo le extrae $id (R/nested_cv.R:130) -> la
#   particion queda FIJA y solo cambian las etiquetas.
#
# Uso:
#   PTB_PERM_MODE=verify Rscript scripts/permutation_test.R
#   PTB_PERM_MODE=run PTB_PERM_N=20  PTB_PERM_CORES=8 Rscript scripts/permutation_test.R
#   PTB_PERM_MODE=run PTB_PERM_N=999 PTB_PERM_CORES=8 Rscript scripts/permutation_test.R
# =============================================================================

suppressWarnings(suppressMessages({ library(here); library(parallel) }))

MODE    <- Sys.getenv("PTB_PERM_MODE", "verify")
N_PERM  <- as.integer(Sys.getenv("PTB_PERM_N", "20"))
N_CORES <- as.integer(Sys.getenv("PTB_PERM_CORES", "1"))
SEED0   <- 10000L
OUT_DIR <- here::here("analysis", "_output")
RMD     <- here::here("analysis", "integrated_preterm_prediction_workflow.Rmd")
TAG     <- Sys.getenv("PTB_PERM_TAG", sprintf("n%d", N_PERM))

# Baseline congelado (Fase 0), real data post-CLR-fix y post-Escherichia-Shigella.
# Fuente: analysis/_output/A_refactor_real_metrics.csv
TARGET <- list(model = "glmnet_base", approach = "Approach3_DataDriven",
               microbiome = "ANCOM_Taxa", frozen_auroc = 0.76)

say <- function(...) { cat(sprintf(...), "\n", sep = ""); flush(stdout()) }

# -----------------------------------------------------------------------------
# Scaffolding: ejecutar los chunks del .Rmd principal, POR NOMBRE
# -----------------------------------------------------------------------------
# Se reejecuta el setup real en vez de reimplementarlo: cualquier divergencia
# (una semilla, un filtro) invalidaria la comparacion observado-vs-nula. Llega
# hasta `nested_cv_helper_functions` porque ese chunk es el que hace source() de
# R/clr.R, R/threshold.R, R/feature_selection.R y R/nested_cv.R; el loop de 12
# combos vive en el chunk siguiente (`run_nested_cv`) y NO se ejecuta.
# Ademas trae los 6 helpers de screening de Approach 3, que son GLOBAL DEPs del
# motor y viven SOLO en el .Rmd, no en R/.
read_rmd_chunks <- function(path) {
  lines  <- readLines(path, warn = FALSE, encoding = "UTF-8")
  starts <- grep("^```\\{r", lines); closes <- grep("^```\\s*$", lines)
  lapply(starts, function(s) {
    e <- closes[closes > s][1]; hdr <- lines[s]
    list(name = trimws(sub("[,}].*$", "", sub("^```\\{r\\s*", "", hdr))), header = hdr,
         code = if (e > s + 1) lines[(s + 1):(e - 1)] else character(0))
  })
}

source_rmd_upto <- function(path, upto, envir = globalenv()) {
  chunks <- read_rmd_chunks(path); nms <- vapply(chunks, `[[`, "", "name")
  stopifnot(upto %in% nms)
  for (i in seq_len(which(nms == upto)[1])) {
    ch <- chunks[[i]]
    if (grepl("eval\\s*=\\s*FALSE", ch$header) || !length(ch$code)) next
    eval(parse(text = paste(ch$code, collapse = "\n")), envir = envir)
  }
  invisible(TRUE)
}

build_scaffolding <- function(rmd) {
  g <- globalenv()
  source_rmd_upto(rmd, "nested_cv_helper_functions", envir = g)
  # clr_zero_levels y abs_data se calculan al inicio del chunk `run_nested_cv`,
  # ANTES del loop. Replicados verbatim (el motor los exige como argumentos).
  assign("clr_taxa_all",
         setdiff(names(get("micro_genus_full", g)), c("index", "id", "shannon_diversity")), g)
  assign("clr_zero_levels",
         get("fit_clr_zerorepl", g)(get("micro_genus_full", g), get("clr_taxa_all", g)), g)
  assign("abs_data", get("load_abs_matrix", g)(get("cfg", g)), g)
  assign("glmnet_spec", get("models_spec", g)$glmnet_base, g)
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Estadistico: AUROC por fold con DIRECCION FIJA, desde las predicciones
# -----------------------------------------------------------------------------
# direction = "<" -> controles(0) por debajo de casos(1), que es la orientacion
# correcta para una probabilidad de la clase "1". A diferencia de "auto", NO mira
# las etiquetas para elegir orientacion, asi que bajo permutacion la nula puede
# caer por debajo de 0.5 y se centra en 0.5.
auroc_fixed_per_fold <- function(preds) {
  if (is.null(preds) || !nrow(preds)) return(numeric(0))
  vapply(split(preds, preds$fold), function(d) {
    tc <- as.character(d$true_class)
    if (length(unique(tc)) < 2) return(NA_real_)   # politica (a): fold no computable
    as.numeric(pROC::auc(pROC::roc(factor(tc, levels = c("0", "1")), d$pred_prob,
                                   levels = c("0", "1"), direction = "<",
                                   quiet = TRUE)))
  }, numeric(1))
}

# Las tres columnas `preterm` tienen TIPOS DISTINTOS (character / numeric /
# integer). Propagar la etiqueta barajada sin respetar el tipo destino rompe el
# motor ("argument is of length zero" en el screening de Approach 3).
cast_like <- function(x, template) {
  if (is.factor(template))       factor(as.character(x), levels = levels(template))
  else if (is.integer(template)) as.integer(x)
  else if (is.numeric(template)) as.numeric(x)
  else                           as.character(x)
}

# -----------------------------------------------------------------------------
# Una corrida del combo objetivo (permutada o no)
# -----------------------------------------------------------------------------
run_target <- function(perm_seed = NULL, keep_full = FALSE) {
  g     <- globalenv()
  clin  <- get("approach3_clinical_all", g)
  abs_d <- get("abs_data", g)
  sl    <- get("subject_labels", g)
  clin0 <- clin; abs0 <- abs_d

  if (!is.null(perm_seed)) {
    set.seed(perm_seed)
    sl$preterm <- sample(sl$preterm)                       # baraja a nivel SUJETO
    lab <- stats::setNames(sl$preterm, sl$id)
    clin$preterm       <- cast_like(unname(lab[as.character(clin$id)]), clin0$preterm)
    abs_d$meta$preterm <- cast_like(unname(lab[as.character(abs_d$meta$id)]),
                                    abs0$meta$preterm)     # ⚠ ANCOM-BC2
    stopifnot(!anyNA(clin$preterm), !anyNA(abs_d$meta$preterm))
  }

  old_sl <- get("subject_labels", g)                        # GLOBAL DEP del motor
  assign("subject_labels", sl, envir = g)
  on.exit(assign("subject_labels", old_sl, envir = g), add = TRUE)

  res <- tryCatch(
    get("train_with_nested_cv", g)(
      model_name = TARGET$model, model_spec = get("glmnet_spec", g),
      clinical_data_all = clin, microbiome_data_all = get("micro_genus_full", g),
      approach_name = TARGET$approach, microbiome_option = TARGET$microbiome,
      cv_folds = get("cv_folds", g), clr_zero_levels = get("clr_zero_levels", g),
      genera_clean = get("genera_clean", g), abs_data = abs_d),
    error = function(e) { message("perm failed: ", conditionMessage(e)); NULL })

  if (is.null(res) || is.null(res$summary))
    return(list(auto = NA_real_, fixed = NA_real_, n_folds_valid = 0L, err = TRUE))

  per_fold <- auroc_fixed_per_fold(res$test_predictions)
  out <- list(auto  = as.numeric(res$summary$AUROC_mean),   # direction="auto"
              fixed = if (all(is.na(per_fold))) NA_real_ else mean(per_fold, na.rm = TRUE),
              n_folds_valid = sum(!is.na(per_fold)),        # politica (a)
              err = FALSE)
  # keep_full solo para el observado: en la nula multiplicaria por 999 la memoria
  if (keep_full) {
    out$fold_results     <- res$fold_results
    out$auroc_fixed_fold <- per_fold
    out$test_predictions <- res$test_predictions
  }
  out
}

# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------
say("\n=== Scaffolding (chunks del .Rmd hasta nested_cv_helper_functions) ===")
t0 <- Sys.time()
invisible(capture.output(suppressWarnings(suppressMessages(build_scaffolding(RMD)))))
say("listo en %.1f s | subjects=%d PTB=%d folds=%d",
    as.numeric(difftime(Sys.time(), t0, units = "secs")),
    nrow(subject_labels), sum(subject_labels$preterm == "1"), nrow(cv_folds))

if (identical(MODE, "verify")) {
  say("\n=== FASE 0/2 — observado (SIN permutar) ===")
  say("Target: %s | %s | %s", TARGET$model, TARGET$approach, TARGET$microbiome)
  tt <- Sys.time()
  invisible(capture.output(obs <- run_target(NULL, keep_full = TRUE)))
  say("\n  AUROC_auto  (motor, direction=auto) : %.4f", obs$auto)
  say("  AUROC_fixed (direction='<' fija)    : %.4f", obs$fixed)
  say("  baseline congelado                  : %.4f", TARGET$frozen_auroc)
  say("  folds computables                   : %d", obs$n_folds_valid)
  say("  |auto - congelado| = %.2e   |fixed - auto| = %.2e",
      abs(obs$auto - TARGET$frozen_auroc), abs(obs$fixed - obs$auto))
  say("  runtime: %.1f s", as.numeric(difftime(Sys.time(), tt, units = "secs")))
  saveRDS(obs, file.path(OUT_DIR, "perm_observed.rds"))

  # TABLA SUPLEMENTARIA POR FOLD — segundo pedido del revisor ("per-fold metrics
  # should appear in a supplementary table"). Portado del chunk
  # `supplementary_fold_metrics` del .Rmd retirado: ahora sale del MISMO objeto
  # que produce el observado, en vez de depender de un `cv_results_all` que
  # tuviera que estar en memoria de otra corrida.
  ft <- obs$fold_results
  ft$AUROC_fixed <- as.numeric(obs$auroc_fixed_fold)[match(ft$fold, as.integer(names(obs$auroc_fixed_fold)))]
  ft <- ft[, c("fold", "n_test", "n_ptb", "threshold", "AUROC", "AUROC_fixed",
               "PRAUC", "Sensitivity", "Specificity", "Balanced_Accuracy", "Youden")]
  num <- setdiff(names(ft), c("fold", "n_test", "n_ptb"))
  ft[num] <- lapply(ft[num], function(x) round(x, 4))
  write.csv(ft, file.path(OUT_DIR, "supplementary_per_fold_metrics.csv"), row.names = FALSE)
  say("\n  Tabla suplementaria por fold -> analysis/_output/supplementary_per_fold_metrics.csv")
  print(ft)

} else if (identical(MODE, "report")) {
  # Reconstruye los productos finales desde el .rds ya calculado (no re-corre las
  # permutaciones). Portado de los chunks `permutation_results` y
  # `permutation_plot` del .Rmd retirado.
  nullf <- file.path(OUT_DIR, sprintf("perm_null_%s.rds", TAG))
  stopifnot(file.exists(nullf))
  r   <- readRDS(nullf)
  obs <- r$observed$fixed
  f   <- r$fixed[!is.na(r$fixed)]
  a   <- r$auto[!is.na(r$auto)]
  pv  <- (sum(f >= obs) + 1) / (length(f) + 1)          # Phipson & Smyth (2010)
  pv_auto <- (sum(a >= obs) + 1) / (length(a) + 1)

  say("\n=== REPORTE (%s) ===", TAG)
  say("  observado           : %.4f", obs)
  say("  nula fija           : %.4f +/- %.4f  (z vs 0.5 = %.2f)",
      mean(f), sd(f), (mean(f) - 0.5) / (sd(f) / sqrt(length(f))))
  say("  fraccion nula < 0.5 : %.3f  (con direction=auto: %.3f)", mean(f < 0.5), mean(a < 0.5))
  say("  p-value (fija)      : %.4f   [%d/%d >= observado]", pv, sum(f >= obs), length(f))
  say("  p-value (auto,previo): %.4f   [%d/%d]  <- lo que habria reportado el estadistico anterior",
      pv_auto, sum(a >= obs), length(a))
  say("  validas %d/%d | folds computables %.3f/5 | perms con fold degenerado %d",
      length(f), r$n_perm, mean(r$n_folds_valid), sum(r$n_folds_valid < nrow(cv_folds)))

  write.csv(data.frame(permutation = seq_along(r$fixed), seed = r$seeds,
                       auroc_fixed = r$fixed, auroc_auto = r$auto,
                       n_folds_valid = r$n_folds_valid),
            file.path(OUT_DIR, sprintf("permutation_null_%s.csv", TAG)), row.names = FALSE)
  write.csv(data.frame(model = TARGET$model, approach = TARGET$approach,
                       microbiome = TARGET$microbiome, observed_auroc = obs,
                       null_mean = mean(f), null_sd = sd(f), null_median = median(f),
                       null_min = min(f), null_max = max(f),
                       n_permutations = r$n_perm, n_valid = length(f),
                       n_exceeding = sum(f >= obs), p_value = pv,
                       p_value_prior_statistic = pv_auto),
            file.path(OUT_DIR, sprintf("permutation_summary_%s.csv", TAG)), row.names = FALSE)

  suppressWarnings(suppressMessages(library(ggplot2)))
  p <- ggplot(data.frame(null_auroc = f), aes(x = null_auroc)) +
    geom_histogram(bins = 40, fill = "grey72", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0.5, colour = "grey35", linetype = "dashed", linewidth = 0.7) +
    geom_vline(xintercept = obs, colour = "#C62828", linewidth = 1.1) +
    annotate("text", x = obs - 0.012, y = Inf, vjust = 1.8, hjust = 1,
             label = sprintf("Observed\nAUROC = %.3f", obs),
             colour = "#C62828", size = 3.5, fontface = "bold") +
    annotate("text", x = 0.5, y = Inf, vjust = 1.4, hjust = 0.5,
             label = "Chance\n(0.5)", colour = "grey35", size = 3, fontface = "italic") +
    annotate("label", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4,
             label = sprintf("p = %.3f\n(%d/%d >= observed)", pv, sum(f >= obs), length(f)),
             fill = "#C8E6C9", size = 3.4, fontface = "bold") +
    labs(title = "Permutation test: model discrimination vs. chance",
         subtitle = sprintf("%s | %s | %s - %d subject-level label permutations (direction-fixed AUROC)",
                            TARGET$model, TARGET$approach, TARGET$microbiome, length(f)),
         x = "Mean AUROC (null distribution)", y = "Count",
         caption = sprintf("Null: %.3f +/- %.3f | Observed: %.3f | one-sided permutation p = %.3f",
                           mean(f), sd(f), obs, pv)) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          plot.caption = element_text(size = 8, colour = "grey40"),
          panel.grid.minor = element_blank())
  ggsave(file.path(OUT_DIR, sprintf("permutation_null_%s.png", TAG)), p,
         width = 8, height = 5, dpi = 300)
  ggsave(file.path(OUT_DIR, sprintf("permutation_null_%s.pdf", TAG)), p, width = 8, height = 5)
  say("  CSVs + figura (.png/.pdf) escritos en analysis/_output/")

} else {
  seeds <- SEED0 + seq_len(N_PERM)
  say("\n=== FASE 2 — %d permutaciones | %d worker(s) ===", N_PERM, N_CORES)
  tt <- Sys.time()

  if (N_CORES > 1) {
    cl <- makeCluster(N_CORES)
    clusterExport(cl, c("read_rmd_chunks", "source_rmd_upto", "build_scaffolding",
                        "auroc_fixed_per_fold", "cast_like", "run_target",
                        "RMD", "TARGET"), envir = environment())
    invisible(clusterEvalQ(cl, {
      # Un worker por core: si ademas cada uno abre threads propios (OpenMP/BLAS/
      # data.table) se sobresuscriben los cores y el throughput cae. Se fija a 1
      # ANTES de cargar los paquetes (build_scaffolding los attachea).
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
                 MKL_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1")
      suppressWarnings(suppressMessages({ library(here); library(parallel) }))
      try(data.table::setDTthreads(1), silent = TRUE)
      invisible(capture.output(suppressWarnings(suppressMessages(build_scaffolding(RMD)))))
      TRUE
    }))
    say("workers listos (%.1f s de setup)", as.numeric(difftime(Sys.time(), tt, units = "secs")))
    out <- parLapply(cl, seeds, function(s) {
      invisible(utils::capture.output(                     # el motor es muy verboso
        r <- suppressWarnings(suppressMessages(run_target(s)))))
      r
    })
    stopCluster(cl)          # explicito: on.exit al final del script tiraba ruido
  } else {
    out <- lapply(seeds, function(s) { invisible(capture.output(r <- run_target(s))); r })
  }

  secs <- as.numeric(difftime(Sys.time(), tt, units = "secs"))
  auto  <- vapply(out, `[[`, numeric(1), "auto")
  fixed <- vapply(out, `[[`, numeric(1), "fixed")
  nfv   <- vapply(out, `[[`, integer(1), "n_folds_valid")

  obs_path <- file.path(OUT_DIR, "perm_observed.rds")
  observed <- if (file.exists(obs_path)) readRDS(obs_path) else NULL

  say("\n--- RUNTIME ---")
  say("  total %.1f min | %.1f s/perm efectivos (%d workers)",
      secs / 60, secs / N_PERM, N_CORES)
  say("  proyeccion 999 perms: %.2f h", secs / N_PERM * 999 / 3600)

  say("")
  say("--- NULA: AUROC_auto (el estadistico ANTERIOR, direction=auto) ---")
  say("  media %.4f | sd %.4f | rango [%.4f, %.4f] | validas %d/%d",
      mean(auto, na.rm = TRUE), sd(auto, na.rm = TRUE),
      min(auto, na.rm = TRUE), max(auto, na.rm = TRUE), sum(!is.na(auto)), N_PERM)

  say("\n--- NULA: AUROC_fixed (direction fija — el correcto) ---")
  say("  media %.4f | sd %.4f | rango [%.4f, %.4f] | validas %d/%d",
      mean(fixed, na.rm = TRUE), sd(fixed, na.rm = TRUE),
      min(fixed, na.rm = TRUE), max(fixed, na.rm = TRUE), sum(!is.na(fixed)), N_PERM)
  say("  >>> chequeo: la nula fija deberia centrarse cerca de 0.5 (la auto, ~0.73)")
  say("  folds computables por perm: media %.2f / %d (politica (a))",
      mean(nfv), nrow(cv_folds))
  say("  perms con algun fold degenerado: %d", sum(nfv < nrow(cv_folds)))

  if (!is.null(observed)) {
    v <- fixed[!is.na(fixed)]
    p <- (sum(v >= observed$fixed) + 1) / (length(v) + 1)   # Phipson & Smyth
    say("\n--- p-value (direction fija) ---")
    say("  observado %.4f | nulas >= observado: %d/%d | p = %.4f",
        observed$fixed, sum(v >= observed$fixed), length(v), p)
  }

  saveRDS(list(auto = auto, fixed = fixed, n_folds_valid = nfv, seeds = seeds,
               n_perm = N_PERM, secs = secs, target = TARGET, observed = observed),
          file.path(OUT_DIR, sprintf("perm_null_%s.rds", TAG)))
}

say("\nDONE (%s)", MODE)
