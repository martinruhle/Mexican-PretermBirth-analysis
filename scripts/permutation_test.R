#!/usr/bin/env Rscript
# =============================================================================
# permutation_test.R — test de permutacion con nula del MAXIMO de 18 combinaciones
# -----------------------------------------------------------------------------
# PREGUNTA: ¿hay discriminacion por encima del azar aun habiendo elegido la mejor
#   de las 18 combinaciones (12 con microbioma + 6 solo clinicas)? No mide el aporte
#   del microbioma (eso es docs/CLINICAL_ONLY_ARM.md). Los criterios de lectura se
#   escribieron ANTES de correr: results/permutation_test/README.md.
#
# ESTADISTICO: por combinacion, el AUROC medio de los folds externos que reporta el
#   motor (roc_ptb(): direccion fija "<", R/nested_cv.R). En cada permutacion se
#   re-corren las 18 (screening del A3, ANCOM-BC2 por fold, split interno, umbral,
#   ajuste) y se guarda el maximo. p ajustado de una combinacion =
#   (#{max_perm >= observado} + 1) / (m + 1)  (Phipson & Smyth 2010; max-T de
#   Westfall-Young en un paso). Se reporta para el maximo observado y para el modelo
#   principal (MAIN). El p sin ajustar de cada combinacion (contra su propia nula)
#   va solo como referencia.
#
# NO-FUGA — se permuta a nivel PARTICIPANTE y la MISMA permutacion llega a todas las
#   copias de la etiqueta: el $clinical de cada enfoque de approaches_config
#   (outcome), subject_labels (GLOBAL DEP del motor: estratifica el split 70/30) y
#   abs_data$meta$preterm (ANCOM-BC2). cv_folds NO se toca: particion fija; un fold
#   de test con una sola clase se saltea (guard del motor), igual en las 18.
#
# CACHE DE ANCOM-BC2: su resultado depende solo de (sujetos de outer-train,
#   etiquetas de esas muestras, OTU filtrada), y es el mismo para las 6
#   combinaciones con ANCOM de una permutacion. Se reemplaza run_ancombc_on_fold()
#   en el entorno global (el motor lo resuelve ahi) por una version con cache cuya
#   llave incluye esas tres cosas, sin tocar R/. El motor re-siembra despues de
#   ANCOM (set.seed(123 + fold + 1000 / + 2000)), asi que la cache no cambia el RNG
#   de lo que sigue. El modo verify lo comprueba (observado y una permutacion).
#
# ANDAMIAJE: chunks del .Rmd principal POR NOMBRE hasta nested_cv_helper_functions,
#   mas el comienzo del chunk run_nested_cv hasta ANTES del loop (ahi se crean
#   models_to_run, approaches_config, clr_zero_levels y abs_data). El brazo solo
#   clinico se arma como en scripts/clinical_only_arm.R: microbioma = solo llaves
#   (index, id), microbiome_option = "Full_Microbiome", genera_clean/abs_data NULL.
#
# MODOS (desde la raiz del repo, con los datos de la cohorte en data/raw/):
#   PTB_PERM_MODE=verify Rscript scripts/permutation_test.R
#     observado de las 18 == corrida de referencia (PTB_REFERENCE_RUN_TAG, por
#     defecto clinonly_2026-10-16) campo a campo y == el CSV versionado del brazo
#     solo clinico; cache == sin cache (observado y semilla SEED0 + 1); tiempos.
#     Guarda analysis/_output/perm_observed_<TAG>.rds.
#   PTB_PERM_MODE=run PTB_PERM_N=199 PTB_PERM_CORES=4 Rscript scripts/permutation_test.R
#     una .rds por semilla en analysis/_output/perm_<TAG>/ (reanudable: las semillas
#     ya guardadas se saltean) + progress.log.
#   PTB_PERM_MODE=report PTB_PERM_N=199 Rscript scripts/permutation_test.R
#     p-values, tablas y figura desde las .rds (no carga el pipeline).
# =============================================================================

suppressWarnings(suppressMessages({ library(here); library(parallel) }))

MODE     <- Sys.getenv("PTB_PERM_MODE", "verify")
N_PERM   <- as.integer(Sys.getenv("PTB_PERM_N", "199"))
N_CORES  <- as.integer(Sys.getenv("PTB_PERM_CORES", "4"))
SEED0    <- 20000L                     # semillas SEED0 + 1..N (la prueba previa: 10001-10999)
TAG      <- Sys.getenv("PTB_PERM_TAG", "max18")
REF_TAG  <- Sys.getenv("PTB_REFERENCE_RUN_TAG", "clinonly_2026-10-16")
OUT_DIR  <- here::here("analysis", "_output")
PERM_DIR <- file.path(OUT_DIR, paste0("perm_", TAG))
OBS_FILE <- file.path(OUT_DIR, paste0("perm_observed_", TAG, ".rds"))
RMD      <- here::here("analysis", "integrated_preterm_prediction_workflow.Rmd")
REF_CSV  <- here::here("results", "clinical_only_arm", "metrics_by_combination.csv")
MAIN     <- "glmnet_base|Approach3_DataDriven|ANCOM_Taxa"   # modelo principal reportado

# Nombres que el motor busca con exists(): si alguno existiera en el entorno global,
# el motor lo usaria (p.ej. microbiome_data_fold entraria como microbioma).
ENGINE_LOCALS <- c("microbiome_data_fold", "microbiome_to_use", "selected_taxa_by_fold",
                   "selected_variables_by_fold", "roc_curves_storage", "pr_curves_storage",
                   "final_model_fit")
# Campos que se comparan entre corridas (los mismos que guarda clinical_only_arm.R)
CMP_FIELDS <- c("model_name", "approach", "microbiome", "fold_results", "summary",
                "test_predictions", "selected_taxa", "selected_variables", "availability")

say <- function(...) { cat(sprintf(...), "\n", sep = ""); flush(stdout()) }
perm_file <- function(seed) file.path(PERM_DIR, sprintf("seed_%05d.rds", seed))

# -----------------------------------------------------------------------------
# Andamiaje: chunks del .Rmd principal, POR NOMBRE
# -----------------------------------------------------------------------------
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
  # Comienzo del chunk run_nested_cv, hasta ANTES del loop de las 12: define
  # models_to_run, approaches_config, microbiome_config, clr_zero_levels, abs_data.
  chunks <- read_rmd_chunks(rmd)
  code   <- chunks[[which(vapply(chunks, `[[`, "", "name") == "run_nested_cv")]]$code
  stop_at <- grep("^cv_results_all <- list\\(\\)", code)
  stopifnot("run_nested_cv: no se encontro el comienzo del loop" = length(stop_at) == 1)
  eval(parse(text = paste(code[seq_len(stop_at - 1)], collapse = "\n")), envir = g)
  # El motor resuelve run_ancombc_on_fold() en su entorno: tiene que ser el global
  # para que la version con cache lo reemplace.
  stopifnot(identical(environment(get("train_with_nested_cv", g)), g))
  assign("micro_none", get("micro_genus_full", g)[, c("index", "id")], g)
  assign("run_ancombc_uncached", get("run_ancombc_on_fold", g), g)
  invisible(TRUE)
}

# Las columnas `preterm` tienen TIPOS DISTINTOS (character / numeric / integer).
# Propagar la etiqueta barajada sin respetar el tipo destino rompe el motor
# ("argument is of length zero" en el screening de Approach 3).
cast_like <- function(x, template) {
  if (is.factor(template))       factor(as.character(x), levels = levels(template))
  else if (is.integer(template)) as.integer(x)
  else if (is.numeric(template)) as.numeric(x)
  else                           as.character(x)
}

# -----------------------------------------------------------------------------
# Cache de ANCOM-BC2
# -----------------------------------------------------------------------------
ancom_cached <- function(otu_data, meta_data, outer_train_subjects) {
  cache <- get(".ancom_cache", envir = globalenv())
  tr  <- meta_data$id %in% outer_train_subjects
  key <- paste(c(sort(unique(as.character(outer_train_subjects))), "|",
                 paste(meta_data$id[tr], meta_data$preterm[tr], sep = ":"), "|",
                 colnames(otu_data), "|", dim(otu_data)), collapse = ",")
  if (exists(key, envir = cache, inherits = FALSE)) {
    cache$.hits <- cache$.hits + 1L
    val <- get(key, envir = cache)
  } else {
    cache$.misses <- cache$.misses + 1L
    val <- tryCatch(get("run_ancombc_uncached", globalenv())(otu_data, meta_data,
                                                             outer_train_subjects),
                    error = function(e) e)
    assign(key, val, envir = cache)
  }
  if (inherits(val, "error")) stop(val)   # mismo error que la llamada sin cache
  val
}

use_ancom_cache <- function(on) {
  g <- globalenv()
  cache <- new.env(hash = TRUE); cache$.hits <- 0L; cache$.misses <- 0L
  assign(".ancom_cache", cache, g)
  assign("run_ancombc_on_fold", if (on) ancom_cached else get("run_ancombc_uncached", g), g)
  invisible(cache)
}

# -----------------------------------------------------------------------------
# Las 18 combinaciones (permutadas o no)
# -----------------------------------------------------------------------------
micro_options <- function() c(get("cfg", globalenv())$microbiome_options, "None")

keys_18 <- function() {
  g <- globalenv()
  unlist(lapply(get("models_to_run", g), function(m)
    lapply(get("approaches_config", g), function(a)
      paste(m$name, a$name, micro_options(), sep = "|"))))
}

permute_inputs <- function(perm_seed) {
  g  <- globalenv()
  sl <- get("subject_labels", g); ap <- get("approaches_config", g); ab <- get("abs_data", g)
  if (!is.null(perm_seed)) {
    set.seed(perm_seed)
    sl$preterm <- sample(sl$preterm)                       # baraja a nivel SUJETO
  }
  # Con o sin permutar, cada copia de la etiqueta sale de la MISMA tabla por sujeto
  # (sin permutar, deja los datos iguales: lo comprueba verify).
  lab <- stats::setNames(as.character(sl$preterm), as.character(sl$id))
  for (i in seq_along(ap)) {
    cl <- ap[[i]]$clinical
    cl$preterm <- cast_like(unname(lab[as.character(cl$id)]), cl$preterm)
    stopifnot(!anyNA(cl$preterm))
    ap[[i]]$clinical <- cl
  }
  ab$meta$preterm <- cast_like(unname(lab[as.character(ab$meta$id)]), ab$meta$preterm)  # ⚠ ANCOM-BC2
  stopifnot(!anyNA(ab$meta$preterm))
  list(subject_labels = sl, approaches = ap, abs_data = ab)
}

run_18 <- function(perm_seed = NULL, cache = TRUE) {
  g   <- globalenv()
  inp <- permute_inputs(perm_seed)
  old_sl <- get("subject_labels", g)                        # GLOBAL DEP del motor
  assign("subject_labels", inp$subject_labels, envir = g)
  on.exit(assign("subject_labels", old_sl, envir = g), add = TRUE)
  ac <- use_ancom_cache(cache)

  res <- list()
  for (model in get("models_to_run", g)) for (approach in inp$approaches)
    for (mi in micro_options()) {
      hit <- ENGINE_LOCALS[vapply(ENGINE_LOCALS, exists, logical(1), envir = g)]
      if (length(hit)) stop("Objects the engine looks up with exists() are defined: ",
                            paste(hit, collapse = ", "))
      none <- identical(mi, "None")
      k <- paste(model$name, approach$name, mi, sep = "|")
      r <- tryCatch(get("train_with_nested_cv", g)(
        model_name = model$name, model_spec = model$spec,
        clinical_data_all = approach$clinical,
        microbiome_data_all = get(if (none) "micro_none" else "micro_genus_full", g),
        approach_name = approach$name,
        microbiome_option = if (none) "Full_Microbiome" else mi,
        cv_folds = get("cv_folds", g), clr_zero_levels = get("clr_zero_levels", g),
        dict = get("dict", g), allowed_availability = get("availability_allowed", g),
        genera_clean = if (none) NULL else get("genera_clean", g),
        abs_data = if (none) NULL else inp$abs_data, use_pca = model$pca),
        error = function(e) { message(k, " failed: ", conditionMessage(e)); NULL })
      if (!is.null(r)) r$microbiome <- mi
      res[k] <- list(r)
    }
  attr(res, "ancom_calls") <- c(hits = ac$.hits, misses = ac$.misses)
  res
}

# AUROC medio por combinacion (el del motor) + chequeo: recalculado desde las
# predicciones out-of-fold con roc_ptb() tiene que dar lo mismo.
summarise_18 <- function(res, keys) {
  auc <- vapply(keys, function(k) {
    r <- res[[k]]; if (is.null(r$summary)) NA_real_ else as.numeric(r$summary$AUROC_mean)
  }, numeric(1))
  nf <- vapply(keys, function(k) {
    r <- res[[k]]; if (is.null(r$summary)) 0L else as.integer(r$summary$n_folds)
  }, integer(1))
  rd <- vapply(keys, function(k) {
    tp <- res[[k]]$test_predictions
    if (is.null(tp) || !nrow(tp)) return(NA_real_)
    pf <- vapply(split(tp, tp$fold), function(d) {
      if (length(unique(as.character(d$true_class))) < 2) return(NA_real_)
      as.numeric(pROC::auc(get("roc_ptb", globalenv())(d$true_class, d$pred_prob)))
    }, numeric(1))
    abs(mean(pf, na.rm = TRUE) - auc[[k]])
  }, numeric(1))
  list(auroc = auc, n_folds = nf,
       recompute_maxdiff = if (all(is.na(rd))) NA_real_ else max(rd, na.rm = TRUE),
       ancom_calls = attr(res, "ancom_calls"))
}

same_results <- function(a, b, keys) {
  bad <- unlist(lapply(keys, function(k) {
    f <- intersect(CMP_FIELDS, union(names(a[[k]]), names(b[[k]])))
    d <- f[!vapply(f, function(x) identical(a[[k]][[x]], b[[k]][[x]]), logical(1))]
    if (length(d)) paste0(k, " [", paste(d, collapse = ", "), "]")
  }))
  if (length(bad)) paste(bad, collapse = "; ") else TRUE
}

timed <- function(expr) {
  t <- Sys.time()
  invisible(utils::capture.output(v <- suppressWarnings(suppressMessages(expr))))
  list(value = v, mins = as.numeric(difftime(Sys.time(), t, units = "mins")))
}

# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------
if (identical(MODE, "verify")) {
  say("=== Scaffolding (chunks del .Rmd hasta nested_cv_helper_functions + inicio de run_nested_cv) ===")
  t0 <- Sys.time()
  invisible(capture.output(suppressWarnings(suppressMessages(build_scaffolding(RMD)))))
  KEYS <- keys_18()
  say("listo en %.1f s | sujetos %d (PTB %d) | folds %d | combinaciones %d | clases: %s",
      as.numeric(difftime(Sys.time(), t0, units = "secs")), nrow(subject_labels),
      sum(subject_labels$preterm == "1"), nrow(cv_folds), length(KEYS),
      paste(availability_allowed, collapse = ", "))
  stopifnot(length(KEYS) == 18, MAIN %in% KEYS)

  # Sin permutar, las etiquetas reconstruidas desde subject_labels == las de los datos
  inp0 <- permute_inputs(NULL)
  stopifnot(identical(inp0$approaches, approaches_config), identical(inp0$abs_data, abs_data))
  say("Etiquetas por sujeto consistentes en los 3 enfoques y en la metadata de ANCOM-BC2: TRUE")

  say("\n=== 1. Observado, sin cache ===")
  o1 <- timed(run_18(NULL, cache = FALSE))
  say("  %.1f min | llamadas a ANCOM-BC2: %s", o1$mins,
      paste(names(attr(o1$value, "ancom_calls")), attr(o1$value, "ancom_calls"), collapse = " "))
  ref_file <- file.path(OUT_DIR, paste0(REF_TAG, "_cv_folds.rds"))
  stopifnot("falta la corrida de referencia" = file.exists(ref_file))
  ref <- readRDS(ref_file)
  names(ref) <- vapply(ref, function(r) paste(r$model_name, r$approach, r$microbiome, sep = "|"), "")
  stopifnot(setequal(names(ref), KEYS))
  chk <- same_results(o1$value, ref, KEYS)
  if (!isTRUE(chk)) stop("El observado difiere de la corrida ", REF_TAG, ": ", chk)
  say("  18 combinaciones == corrida %s (%s): TRUE", REF_TAG, paste(CMP_FIELDS[-(1:3)], collapse = ", "))
  rc <- utils::read.csv(REF_CSV, stringsAsFactors = FALSE)
  rck <- paste(rc$Model, rc$Approach, rc$Microbiome, sep = "|")
  obs <- summarise_18(o1$value, KEYS)
  stopifnot(setequal(rck, KEYS),
            isTRUE(all.equal(unname(obs$auroc[rck]), rc$AUROC_mean, tolerance = 1e-12)))
  say("  AUROC == results/clinical_only_arm/metrics_by_combination.csv: TRUE")
  say("  AUROC del motor == recalculado desde las predicciones (max |diff| %.1e)", obs$recompute_maxdiff)
  stopifnot(obs$recompute_maxdiff < 1e-12)

  say("\n=== 2. Observado, con cache ===")
  o2 <- timed(run_18(NULL, cache = TRUE))
  ac <- attr(o2$value, "ancom_calls")
  chk <- same_results(o2$value, o1$value, KEYS)
  if (!isTRUE(chk)) stop("Con cache el observado difiere: ", chk)
  say("  %.1f min | ANCOM-BC2: %d calculadas, %d desde la cache | == sin cache: TRUE",
      o2$mins, ac[["misses"]], ac[["hits"]])
  stopifnot(ac[["hits"]] > 0)

  s1 <- SEED0 + 1L
  say("\n=== 3. Permutacion %d, sin cache y con cache ===", s1)
  p1 <- timed(run_18(s1, cache = FALSE))
  p2 <- timed(run_18(s1, cache = TRUE))
  chk <- same_results(p2$value, p1$value, KEYS)
  if (!isTRUE(chk)) stop("Con cache la permutacion difiere: ", chk)
  ac <- attr(p2$value, "ancom_calls")
  say("  sin cache %.1f min | con cache %.1f min (ANCOM-BC2: %d calculadas, %d desde la cache) | identicas: TRUE",
      p1$mins, p2$mins, ac[["misses"]], ac[["hits"]])

  obs$keys <- KEYS; obs$main <- MAIN; obs$reference_run <- REF_TAG
  obs$fold_results <- lapply(o1$value[KEYS], `[[`, "fold_results")
  obs$minutes <- c(observed_uncached = o1$mins, observed_cached = o2$mins,
                   perm_uncached = p1$mins, perm_cached = p2$mins)
  saveRDS(obs, OBS_FILE)
  say("\n  maximo observado %.4f (%s) | principal %.4f (%s)",
      max(obs$auroc), KEYS[which.max(obs$auroc)], obs$auroc[[MAIN]], MAIN)
  say("  -> %s", OBS_FILE)
  say("  proyeccion de %d permutaciones con cache: %.1f h de nucleo", N_PERM, p2$mins * N_PERM / 60)

} else if (identical(MODE, "run")) {
  seeds <- SEED0 + seq_len(N_PERM)
  dir.create(PERM_DIR, showWarnings = FALSE, recursive = TRUE)
  todo <- seeds[!file.exists(perm_file(seeds))]
  prog <- file.path(PERM_DIR, "progress.log")
  say("=== %d permutaciones (semillas %d-%d) | ya guardadas %d | faltan %d | %d worker(s) ===",
      N_PERM, min(seeds), max(seeds), N_PERM - length(todo), length(todo), N_CORES)
  cat(sprintf("%s start: %d to do, %d workers\n", format(Sys.time()), length(todo), N_CORES),
      file = prog, append = TRUE)

  one_seed <- function(s) {
    t <- Sys.time()
    out <- tryCatch({
      invisible(utils::capture.output(                     # el motor es muy verboso
        res <- suppressWarnings(suppressMessages(run_18(s, cache = TRUE)))))
      sm <- summarise_18(res, keys_18())
      sm$seed <- s; sm$minutes <- as.numeric(difftime(Sys.time(), t, units = "mins"))
      saveRDS(sm, perm_file(s))
      sprintf("seed %d ok %.1f min | max %.4f | folds %d | ANCOM %d/%d", s, sm$minutes,
              max(sm$auroc, na.rm = TRUE), min(sm$n_folds), sm$ancom_calls[["misses"]],
              sum(sm$ancom_calls))
    }, error = function(e) sprintf("seed %d FAILED: %s", s, conditionMessage(e)))
    cat(sprintf("%s pid %d %s\n", format(Sys.time()), Sys.getpid(), out), file = prog, append = TRUE)
    out
  }

  tt <- Sys.time()
  if (length(todo) && N_CORES > 1) {
    cl <- makeCluster(N_CORES, outfile = file.path(PERM_DIR, "workers.log"))
    clusterExport(cl, c("read_rmd_chunks", "source_rmd_upto", "build_scaffolding", "cast_like",
                        "ancom_cached", "use_ancom_cache", "micro_options", "keys_18",
                        "permute_inputs", "run_18", "summarise_18", "perm_file", "one_seed",
                        "RMD", "PERM_DIR", "ENGINE_LOCALS", "prog"), envir = environment())
    invisible(clusterEvalQ(cl, {
      # Un worker por core: sin threads propios (OpenMP/BLAS/data.table) para no
      # sobresuscribir. Se fija ANTES de cargar los paquetes.
      Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
                 MKL_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1")
      suppressWarnings(suppressMessages({ library(here); library(parallel) }))
      try(data.table::setDTthreads(1), silent = TRUE)
      invisible(capture.output(suppressWarnings(suppressMessages(build_scaffolding(RMD)))))
      TRUE
    }))
    say("workers listos (%.1f min de setup)", as.numeric(difftime(Sys.time(), tt, units = "mins")))
    out <- parLapplyLB(cl, todo, one_seed)
    stopCluster(cl)
  } else if (length(todo)) {
    invisible(capture.output(suppressWarnings(suppressMessages(build_scaffolding(RMD)))))
    out <- lapply(todo, one_seed)
  } else out <- list()
  done <- sum(file.exists(perm_file(seeds)))
  say("listo en %.1f h | guardadas %d/%d | fallidas en esta corrida %d",
      as.numeric(difftime(Sys.time(), tt, units = "hours")), done, N_PERM,
      sum(grepl("FAILED", unlist(out))))

} else if (identical(MODE, "report")) {
  stopifnot(file.exists(OBS_FILE))
  obs   <- readRDS(OBS_FILE); keys <- obs$keys
  seeds <- SEED0 + seq_len(N_PERM)
  have  <- file.exists(perm_file(seeds))
  P <- lapply(perm_file(seeds[have]), readRDS)
  stopifnot(length(P) > 0, all(vapply(P, function(p) identical(names(p$auroc), keys), logical(1))))
  A  <- do.call(rbind, lapply(P, `[[`, "auroc"))           # m x 18
  nf <- vapply(P, function(p) min(p$n_folds), integer(1))
  mx <- apply(A, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  ok <- !is.na(mx); m <- sum(ok)
  eps <- 1e-12                                             # empates exactos cuentan como >=
  p_adj <- function(x) (sum(mx[ok] >= x - eps) + 1) / (m + 1)
  mc_ci <- function(x) stats::binom.test(sum(mx[ok] >= x - eps), m)$conf.int
  best_key <- keys[which.max(obs$auroc)]

  summ <- do.call(rbind, lapply(keys, function(k) {
    a <- A[, k]; v <- a[!is.na(a)]; ci <- mc_ci(obs$auroc[[k]])
    data.frame(combination = k, observed_auroc = obs$auroc[[k]],
               null_mean = mean(v), null_sd = stats::sd(v), n_null = length(v),
               p_unadjusted = (sum(v >= obs$auroc[[k]] - eps) + 1) / (length(v) + 1),
               p_adjusted_max18 = p_adj(obs$auroc[[k]]),
               p_adjusted_mc_ci_low = ci[1], p_adjusted_mc_ci_high = ci[2],
               stringsAsFactors = FALSE)
  }))
  summ <- summ[order(-summ$observed_auroc), ]

  say("=== REPORTE (%s) | permutaciones con resultado %d/%d (guardadas %d) ===",
      TAG, m, N_PERM, sum(have))
  say("  nula del maximo de 18: media %.4f | sd %.4f | mediana %.4f | rango [%.4f, %.4f]",
      mean(mx[ok]), stats::sd(mx[ok]), stats::median(mx[ok]), min(mx[ok]), max(mx[ok]))
  for (k in unique(c(best_key, obs$main))) {
    ci <- mc_ci(obs$auroc[[k]])
    say("  %s: observado %.4f | %d/%d maximos >= observado | p ajustado = %.4f (IC 95%% Monte Carlo %.4f-%.4f)",
        k, obs$auroc[[k]], sum(mx[ok] >= obs$auroc[[k]] - eps), m, p_adj(obs$auroc[[k]]), ci[1], ci[2])
  }
  say("  folds computables por permutacion: media %.2f | permutaciones con algun fold salteado %d",
      mean(nf), sum(nf < 5))
  say("  combinaciones sin resultado: %d de %d", sum(is.na(A)), length(A))
  rdm <- max(vapply(P, function(p) p$recompute_maxdiff, numeric(1)), na.rm = TRUE)
  say("  AUROC del motor vs recalculado desde predicciones: max |diff| %.1e", rdm)
  say("  minutos por permutacion (un worker): mediana %.1f",
      stats::median(vapply(P, `[[`, numeric(1), "minutes")))

  nul <- data.frame(seed = seeds[have], n_folds_valid = nf, max_auroc = mx,
                    max_combination = keys[apply(A, 1, which.max)], A, check.names = FALSE)
  utils::write.csv(nul, file.path(OUT_DIR, sprintf("permutation_null_%s.csv", TAG)), row.names = FALSE)
  utils::write.csv(summ, file.path(OUT_DIR, sprintf("permutation_summary_%s.csv", TAG)), row.names = FALSE)

  suppressWarnings(suppressMessages(library(ggplot2)))
  lab <- function(k) sprintf("%s\n%.3f, p = %.3f", c(best = "Best of 18", main = "Main model")[[
    if (k == best_key) "best" else "main"]], obs$auroc[[k]], p_adj(obs$auroc[[k]]))
  vl <- data.frame(k = unique(c(best_key, obs$main)))
  vl$x <- unname(obs$auroc[vl$k]); vl$label <- vapply(vl$k, lab, "")
  p <- ggplot(data.frame(max_auroc = mx[ok]), aes(x = max_auroc)) +
    geom_histogram(bins = 30, fill = "grey72", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0.5, colour = "grey35", linetype = "dashed", linewidth = 0.6) +
    geom_vline(data = vl, aes(xintercept = x), colour = c("#C62828", "#1565C0")[seq_len(nrow(vl))],
               linewidth = 1) +
    geom_label(data = vl, aes(x = x, y = Inf, label = label), vjust = 1.2, size = 3,
               fill = "white", label.size = 0) +
    labs(title = "Permutation test: best of 18 combinations vs. chance",
         subtitle = sprintf("Maximum mean AUROC over the 18 combinations in each of %d subject-level label permutations",
                            m),
         x = "Maximum mean AUROC across the 18 combinations (null distribution)", y = "Count",
         caption = sprintf("Null maximum: %.3f ± %.3f | one-sided p = (b + 1)/(m + 1)",
                           mean(mx[ok]), stats::sd(mx[ok]))) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          plot.caption = element_text(size = 8, colour = "grey40"),
          panel.grid.minor = element_blank())
  ggsave(file.path(OUT_DIR, sprintf("permutation_null_%s.png", TAG)), p, width = 8, height = 5, dpi = 300)
  ggsave(file.path(OUT_DIR, sprintf("permutation_null_%s.pdf", TAG)), p, width = 8, height = 5)
  say("  CSVs + figura escritos en analysis/_output/")
} else {
  stop("PTB_PERM_MODE desconocido: ", MODE)
}

say("DONE (%s)", MODE)
