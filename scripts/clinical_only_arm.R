#!/usr/bin/env Rscript
# =============================================================================
# clinical_only_arm.R  —  Tablas de docs/CLINICAL_ONLY_ARM.md
# -----------------------------------------------------------------------------
# QUE ES: las 12 combinaciones del analisis principal (2 modelos x 3 enfoques x 2
# entradas de microbioma) incluyen siempre microbioma, asi que no separan su aporte.
# Este script agrega el brazo SOLO CLINICO: los mismos 2 modelos x 3 enfoques, con
# las mismas variables clinicas, los mismos folds, las mismas semillas y el mismo
# motor, sin ningun taxon ni la diversidad de Shannon. La comparacion fold a fold es
# pareada: el motor siembra cada fold con 123 + fold (split interno + 1000, ajuste
# + 2000) sin importar la entrada de microbioma, asi que los dos brazos ven los
# mismos sujetos en inner-train, inner-val y outer-test.
#
# COMO (sin tocar R/, el .Rmd ni config/):
#   1. reconstruye el entorno del pipeline corriendo los chunks del .Rmd POR NOMBRE
#      hasta run_nested_cv inclusive, con el perfil del config (clases de
#      disponibilidad at_visit por defecto). Ese chunk corre las 12 combinaciones
#      con microbioma tal como lo hace el pipeline.
#   2. llama a train_with_nested_cv() con cada modelo y enfoque definidos en ese
#      mismo chunk, pasandole como microbioma solo las llaves (index, id). El motor
#      toma como taxa toda columna que no sea llave, desenlace ni variable clinica,
#      asi que no encuentra ninguna y no aplica la CLR; con microbiome_option =
#      "Full_Microbiome" no corre ANCOM-BC2 (genera_clean y abs_data van en NULL:
#      si lo intentara, abortaria). Se etiqueta Microbiome = "None".
#   3. referencia: numero de muestras por participante (docs/SAMPLE_COUNT_REFERENCE.md),
#      score = -(numero de muestras), AUROC y PR-AUC con roc_ptb() / prauc_ptb() del
#      motor en los mismos folds.
#
# CHEQUEOS (abortan si fallan):
#   a. las 12 con microbioma reproducen una corrida guardada (PTB_REFERENCE_RUN_TAG,
#      por defecto la linea base sesfix_2026-10-16, si su _cv_folds.rds existe):
#      metricas por fold, predicciones, taxa de ANCOM-BC2 y variables por fold
#      identicas; y, con el perfil default, las filas "after" de
#      results/socioeconomic_indicators/metrics_by_combination.csv (versionado);
#   b. el brazo solo clinico no ve microbioma: antes de cada llamada no existe ningun
#      objeto con los nombres que el motor busca con exists() (exists() tambien mira
#      los entornos que lo envuelven, y un objeto suelto con esos nombres entraria al
#      modelo); el log de cada corrida dice 0 columnas de taxa en los 5 folds y en el
#      modelo final; los predictores del modelo final son exactamente las variables
#      clinicas del enfoque, sin taxa, shannon_diversity ni columnas CLR_;
#   c. cada combinacion solo clinica tiene resultado en los 5 folds y usa las mismas
#      variables clinicas que sus dos pares con microbioma (lista fija en los
#      enfoques 1 y 2; la misma seleccion por fold en el 3) y los mismos sujetos de
#      test por fold;
#   d. la referencia del conteo reproduce el AUROC por fold publicado en
#      results/sample_count_reference/auroc_by_fold.csv (variante "all").
#
# SALIDAS: tablas en markdown por consola; con PTB_RESULTS_DIR, los CSV; en
# analysis/_output/, <PTB_RUN_TAG>.log (la salida del pipeline y del motor) y
# <PTB_RUN_TAG>_cv_folds.rds (las 18 combinaciones, mismos campos que run_baseline.R).
#
# CORRER (necesita los datos de la cohorte en data/raw/; desde la raiz del repo):
#   PTB_RESULTS_DIR=results/clinical_only_arm Rscript scripts/clinical_only_arm.R
# =============================================================================

suppressWarnings(suppressMessages(library(here)))

RMD         <- here::here("analysis", "integrated_preterm_prediction_workflow.Rmd")
profile     <- Sys.getenv("PTB_PROFILE", "default")
run_tag     <- Sys.getenv("PTB_RUN_TAG", "clinonly")
ref_tag     <- Sys.getenv("PTB_REFERENCE_RUN_TAG", "sesfix_2026-10-16")
results_dir <- Sys.getenv("PTB_RESULTS_DIR", "")
ref_csv     <- here::here("results", "socioeconomic_indicators", "metrics_by_combination.csv")
count_csv   <- here::here("results", "sample_count_reference", "auroc_by_fold.csv")

say <- function(...) { cat(sprintf(...), "\n", sep = ""); flush(stdout()) }
f3  <- function(x) sprintf("%.3f", x)
pm  <- function(m, s) sprintf("%.3f ± %.3f", m, s)
sgn <- function(x) sprintf("%+.3f", x)

io <- new.env()
sys.source(here::here("R", "io.R"), envir = io)
outdir <- io$read_config(profile)$output$dir %||% "analysis/_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(outdir, paste0(run_tag, ".log"))
cat("", file = log_file)

# -----------------------------------------------------------------------------
# Chunks del .Rmd principal, POR NOMBRE (mismo lector que scripts/permutation_test.R)
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

source_rmd_upto <- function(path, upto, envir) {
  chunks <- read_rmd_chunks(path); nms <- vapply(chunks, `[[`, "", "name")
  stopifnot(upto %in% nms)
  for (i in seq_len(which(nms == upto)[1])) {
    ch <- chunks[[i]]
    if (grepl("eval\\s*=\\s*FALSE", ch$header) || !length(ch$code)) next
    eval(parse(text = paste(ch$code, collapse = "\n")), envir = envir)
  }
  invisible(TRUE)
}

# Entorno propio para el pipeline (como el envir de render() en run_baseline.R): los
# objetos de los chunks no se mezclan con los de este script.
S <- new.env(parent = globalenv())

say("=== Pipeline chunks up to run_nested_cv: the 12 combinations with microbiome ===")
t0 <- Sys.time()
local({
  con <- file(log_file, open = "a"); sink(con)
  on.exit({ sink(); close(con) })
  suppressWarnings(suppressMessages(source_rmd_upto(RMD, "run_nested_cv", S)))
})
say("done in %.1f min | subjects %d (preterm %d) | outer folds %d | availability classes: %s",
    as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(S$subject_labels),
    sum(S$subject_labels$preterm == "1"), nrow(S$cv_folds),
    paste(S$availability_allowed, collapse = ", "))

key <- function(r) paste(r$model_name, r$approach, r$microbiome, sep = "|")
micro_res <- S$cv_results_all
stopifnot("the pipeline did not produce the 12 combinations" =
            length(micro_res) == 12 && !any(vapply(micro_res, is.null, logical(1))))
names(micro_res) <- vapply(micro_res, key, "")

# -----------------------------------------------------------------------------
# Brazo solo clinico
# -----------------------------------------------------------------------------
# Nombres que el motor busca con exists(): si alguno existiera en S o mas arriba,
# el motor lo usaria (p.ej. microbiome_data_fold entraria como microbioma).
engine_locals <- c("microbiome_data_fold", "microbiome_to_use", "selected_taxa_by_fold",
                   "selected_variables_by_fold", "roc_curves_storage", "pr_curves_storage",
                   "final_model_fit")
taxa_all <- S$clr_taxa_all
micro_none <- S$micro_genus_full[, c("index", "id")]
stopifnot(nrow(micro_none) == nrow(S$micro_genus_full), !anyDuplicated(micro_none$index))

say("\n=== Clinical-only arm: 2 models x 3 approaches, microbiome input = keys only ===")
none_res <- list()
t1 <- Sys.time()
for (model in S$models_to_run) for (approach in S$approaches_config) {
  hit <- engine_locals[vapply(engine_locals, exists, logical(1), envir = S)]
  if (length(hit)) stop("Objects the engine looks up with exists() are defined: ",
                        paste(hit, collapse = ", "))
  k <- paste(model$name, approach$name, "None", sep = "|")
  out <- utils::capture.output(res <- suppressWarnings(suppressMessages(
    S$train_with_nested_cv(
      model_name = model$name, model_spec = model$spec,
      clinical_data_all = approach$clinical, microbiome_data_all = micro_none,
      approach_name = approach$name, microbiome_option = "Full_Microbiome",
      cv_folds = S$cv_folds, clr_zero_levels = S$clr_zero_levels,
      dict = S$dict, allowed_availability = S$availability_allowed,
      genera_clean = NULL, abs_data = NULL, use_pca = model$pca))))
  cat(sprintf("\n######## CLINICAL ONLY: %s ########\n", k), out, sep = "\n",
      file = log_file, append = TRUE)
  if (is.null(res)) stop(k, ": the engine returned no result")

  # b. sin microbioma: 0 columnas de taxa en cada fold y en el modelo final
  n_fold0 <- sum(grepl("Taxa columns identified for CLR: 0$", out))
  n_fold_any <- sum(grepl("Taxa columns identified for CLR: ", out))
  n_final0 <- sum(grepl("Microbiome taxa: 0$", out))
  if (n_fold0 != nrow(S$cv_folds) || n_fold_any != n_fold0 || n_final0 != 1)
    stop(k, ": taxa columns found (folds with 0 taxa: ", n_fold0, "/", n_fold_any,
         "; final model with 0 taxa: ", n_final0, ")")
  if (is.null(res$final_model)) stop(k, ": no final model")
  vi <- workflows::extract_recipe(res$final_model)$var_info
  pred_vars <- vi$variable[vi$role == "predictor"]
  if (!setequal(pred_vars, res$availability$clinical_variables) ||
      length(intersect(pred_vars, c(taxa_all, "shannon_diversity"))) ||
      any(startsWith(pred_vars, "CLR_")))
    stop(k, ": the final model's predictors are not the approach's clinical variables")

  res$microbiome <- "None"
  none_res[[k]] <- res
}
say("done in %.1f min", as.numeric(difftime(Sys.time(), t1, units = "mins")))

all_res <- c(none_res, micro_res)
fold_ids <- function(r) {
  tp <- r$test_predictions
  paste(tp$fold, tp$id)[order(tp$fold, tp$id)]
}

cat("\n== Checks ==\n")
cat("Clinical-only arm: 0 taxa columns in every fold and in the final model; final-model",
    "predictors = the approach's clinical variables: TRUE\n")

# c. mismas variables clinicas y mismos sujetos de test que los pares con microbioma
for (k in names(none_res)) {
  r <- none_res[[k]]
  if (nrow(r$fold_results) != nrow(S$cv_folds)) stop(k, ": fewer than 5 folds with results")
  for (m in c("ANCOM_Taxa", "Full_Microbiome")) {
    p <- micro_res[[sub("\\|None$", paste0("|", m), k)]]
    if (!identical(r$availability, p$availability))
      stop(k, ": clinical variables differ from ", m)
    if (!identical(r$selected_variables, p$selected_variables))
      stop(k, ": Approach 3 per-fold selection differs from ", m)
    if (!identical(fold_ids(r), fold_ids(p)))
      stop(k, ": outer-test subjects differ from ", m)
  }
}
cat("Clinical-only combinations: 5/5 folds; same clinical variables (fixed list or per-fold",
    "selection) and same outer-test subjects as their two microbiome counterparts: TRUE\n")

# a. las 12 con microbioma == corrida guardada y CSV versionado
ref_file <- file.path(outdir, paste0(ref_tag, "_cv_folds.rds"))
if (file.exists(ref_file)) {
  ref <- readRDS(ref_file)
  names(ref) <- vapply(ref, key, "")
  stopifnot(setequal(names(ref), names(micro_res)))
  bad <- unlist(lapply(names(ref), function(k) {
    f <- names(ref[[k]])
    d <- f[!vapply(f, function(x) identical(micro_res[[k]][[x]], ref[[k]][[x]]), logical(1))]
    if (length(d)) paste0(k, " [", paste(d, collapse = ", "), "]")
  }))
  if (length(bad)) stop("The 12 combinations differ from run ", ref_tag, ": ",
                        paste(bad, collapse = "; "))
  cat(sprintf("The 12 combinations with microbiome are identical to run %s (%s): TRUE\n",
              ref_tag, paste(names(ref[[1]])[-(1:3)], collapse = ", ")))
} else {
  cat(sprintf("Run %s not found in %s: check against the saved run SKIPPED\n", ref_tag, outdir))
}
if (identical(profile, "default") && file.exists(ref_csv)) {
  rc <- utils::read.csv(ref_csv, stringsAsFactors = FALSE)
  rc <- rc[rc$run == "after", ]
  cols <- grep("_(mean|sd)$", names(rc), value = TRUE)
  ok <- vapply(seq_len(nrow(rc)), function(i) {
    s <- micro_res[[paste(rc$Model[i], rc$Approach[i], rc$Microbiome[i], sep = "|")]]$summary
    isTRUE(all.equal(unlist(rc[i, cols]), unlist(s[1, cols]), tolerance = 1e-10,
                     check.attributes = FALSE))
  }, logical(1))
  if (nrow(rc) != 12 || !all(ok)) stop("The 12 combinations differ from ", ref_csv)
  cat("The 12 combinations match the 'after' rows of",
      "results/socioeconomic_indicators/metrics_by_combination.csv: TRUE\n")
}

# -----------------------------------------------------------------------------
# Referencia: numero de muestras por participante, mismos folds
# -----------------------------------------------------------------------------
folds <- do.call(rbind, lapply(seq_len(nrow(S$cv_folds)), function(i)
  data.frame(id = as.character(rsample::assessment(S$cv_folds$splits[[i]])$id), fold = i,
             stringsAsFactors = FALSE)))
n_samp <- table(as.character(S$approach1_clinical_all$id))
stopifnot(sum(n_samp) == nrow(S$micro_genus_full), setequal(names(n_samp), folds$id))
cnt <- merge(folds, data.frame(id = as.character(S$subject_labels$id),
                               preterm = as.character(S$subject_labels$preterm),
                               stringsAsFactors = FALSE), by = "id")
cnt$n_samples <- as.integer(n_samp[cnt$id])
cnt$score <- -cnt$n_samples
count_fold <- do.call(rbind, lapply(split(cnt, cnt$fold), function(d) data.frame(
  fold = d$fold[1], n_test = nrow(d), n_ptb = sum(d$preterm == "1"),
  AUROC = as.numeric(pROC::auc(S$roc_ptb(d$preterm, d$score))),
  PRAUC = S$prauc_ptb(d$preterm, d$score))))
if (file.exists(count_csv)) {
  pub <- utils::read.csv(count_csv, stringsAsFactors = FALSE)
  pub <- pub[pub$variant == "all", ]
  if (!isTRUE(all.equal(pub$AUROC[order(pub$fold)], count_fold$AUROC, tolerance = 1e-12)))
    stop("The sample-count reference does not reproduce ", count_csv)
  cat("Sample-count reference reproduces the published per-fold AUROC",
      "(results/sample_count_reference/auroc_by_fold.csv, variant 'all'): TRUE\n")
}

# -----------------------------------------------------------------------------
# Tablas
# -----------------------------------------------------------------------------
model_lab <- c(glmnet_base = "Elastic net (glmnet)", rf_base = "Random Forest")
appr_lab  <- c(Approach1_DREAM = "1 — DREAM (minimal)", Approach2_Literature = "2 — literature",
               Approach3_DataDriven = "3 — data-driven")
micro_lab <- c(None = "Clinical only", ANCOM_Taxa = "+ ANCOM-BC2 taxa",
               Full_Microbiome = "+ full microbiome")
micro_ord <- names(micro_lab)
sm <- c("AUROC", "PRAUC", "Sensitivity", "Specificity", "Balanced_Accuracy")

metrics <- do.call(rbind, lapply(all_res, function(r) {
  s <- r$summary
  cv <- if (r$approach == "Approach3_DataDriven") "selected per fold (screening of the pool)" else
    paste(r$availability$clinical_variables, collapse = " + ")
  data.frame(Model = r$model_name, Approach = r$approach, Microbiome = r$microbiome,
             clinical_variables = cv, n_folds = s$n_folds,
             s[, paste0(rep(sm, each = 2), c("_mean", "_sd"))], stringsAsFactors = FALSE)
}))
metrics <- metrics[order(metrics$Approach, metrics$Model,
                         match(metrics$Microbiome, micro_ord)), ]
rownames(metrics) <- NULL

by_fold <- do.call(rbind, lapply(all_res, function(r) {
  f <- r$fold_results
  data.frame(Model = r$model_name, Approach = r$approach, Microbiome = r$microbiome,
             fold = f$fold, n_test = f$n_test, n_ptb = f$n_ptb, threshold = f$threshold,
             f[, sm], stringsAsFactors = FALSE)
}))
by_fold <- by_fold[order(by_fold$Approach, by_fold$Model,
                         match(by_fold$Microbiome, micro_ord), by_fold$fold), ]
rownames(by_fold) <- NULL

paired <- do.call(rbind, lapply(names(none_res), function(k) {
  b <- by_fold[paste(by_fold$Model, by_fold$Approach, by_fold$Microbiome, sep = "|") == k, ]
  do.call(rbind, lapply(c("ANCOM_Taxa", "Full_Microbiome"), function(m) {
    w <- by_fold[by_fold$Model == b$Model[1] & by_fold$Approach == b$Approach[1] &
                   by_fold$Microbiome == m, ]
    w <- w[match(b$fold, w$fold), ]
    da <- w$AUROC - b$AUROC; dp <- w$PRAUC - b$PRAUC
    data.frame(Model = b$Model[1], Approach = b$Approach[1], Microbiome = m,
               AUROC_clinical_only = mean(b$AUROC), AUROC_with_microbiome = mean(w$AUROC),
               delta_AUROC_mean = mean(da), delta_AUROC_sd = stats::sd(da),
               delta_AUROC_min = min(da), delta_AUROC_max = max(da),
               folds_higher = sum(da > 0), folds_lower = sum(da < 0), folds_equal = sum(da == 0),
               PRAUC_clinical_only = mean(b$PRAUC), PRAUC_with_microbiome = mean(w$PRAUC),
               delta_PRAUC_mean = mean(dp), delta_PRAUC_sd = stats::sd(dp),
               stringsAsFactors = FALSE)
  }))
}))
paired <- paired[order(paired$Approach, paired$Model, match(paired$Microbiome, micro_ord)), ]
rownames(paired) <- NULL

get_m <- function(mo, ap, mi) metrics[metrics$Model == mo & metrics$Approach == ap &
                                        metrics$Microbiome == mi, ]
models <- names(model_lab); approaches <- names(appr_lab)
cnt_m <- c(mean(count_fold$AUROC), stats::sd(count_fold$AUROC),
           mean(count_fold$PRAUC), stats::sd(count_fold$PRAUC))

cat("\n== Approach 1 (gestational age at the visit + maternal age): mean ± SD across folds ==\n\n")
cat("| Model | Input | AUROC | PR-AUC |\n|---|---|---|---|\n")
for (mo in models) for (mi in micro_ord) {
  r <- get_m(mo, "Approach1_DREAM", mi)
  cat(sprintf("| %s | %s | %s | %s |\n", model_lab[[mo]], micro_lab[[mi]],
              pm(r$AUROC_mean, r$AUROC_sd), pm(r$PRAUC_mean, r$PRAUC_sd)))
}
cat(sprintf("| Reference | Number of samples per participant | %s | %s |\n",
            pm(cnt_m[1], cnt_m[2]), pm(cnt_m[3], cnt_m[4])))

cat("\n== Approach 1, AUROC by outer fold ==\n\n")
cat("| Fold | Subjects (preterm) |", paste(sprintf("%s: %s", rep(c("glmnet", "RF"), each = 3),
                                                   rep(c("clinical", "+ANCOM", "+full"), 2)),
                                           collapse = " | "), "| Sample count |\n")
cat("|---|---|", paste(rep("---", 7), collapse = "|"), "|\n", sep = "")
for (fo in sort(unique(by_fold$fold))) {
  v <- unlist(lapply(models, function(mo) lapply(micro_ord, function(mi)
    by_fold$AUROC[by_fold$Model == mo & by_fold$Approach == "Approach1_DREAM" &
                    by_fold$Microbiome == mi & by_fold$fold == fo])))
  cf <- count_fold[count_fold$fold == fo, ]
  cat(sprintf("| %d | %d (%d) | %s | %s |\n", fo, cf$n_test, cf$n_ptb,
              paste(f3(v), collapse = " | "), f3(cf$AUROC)))
}

for (met in c("AUROC", "PRAUC")) {
  cat(sprintf("\n== %s, mean ± SD across folds: clinical only vs with microbiome ==\n\n", met))
  cat("| Approach | Model |", paste(micro_lab[micro_ord], collapse = " | "), "|\n")
  cat("|---|---|---|---|---|\n")
  for (ap in approaches) for (mo in models) {
    v <- vapply(micro_ord, function(mi) {
      r <- get_m(mo, ap, mi); pm(r[[paste0(met, "_mean")]], r[[paste0(met, "_sd")]])
    }, "")
    cat(sprintf("| %s | %s | %s |\n", appr_lab[[ap]], model_lab[[mo]], paste(v, collapse = " | ")))
  }
  cat(sprintf("| Reference | Number of samples per participant | %s | | |\n",
              if (met == "AUROC") pm(cnt_m[1], cnt_m[2]) else pm(cnt_m[3], cnt_m[4])))
}

cat("\n== Paired difference by fold (with microbiome minus clinical only) ==\n\n")
cat("| Approach | Model | Microbiome input | ΔAUROC mean (min to max) | Folds higher / lower / equal | ΔPR-AUC mean |\n")
cat("|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(paired))) with(paired[i, ], cat(sprintf(
  "| %s | %s | %s | %s (%s to %s) | %d / %d / %d | %s |\n", appr_lab[[Approach]], model_lab[[Model]],
  micro_lab[[Microbiome]], sgn(delta_AUROC_mean), sgn(delta_AUROC_min), sgn(delta_AUROC_max),
  folds_higher, folds_lower, folds_equal, sgn(delta_PRAUC_mean))))

cat("\n== Sensitivity and specificity (threshold chosen on the inner validation set) ==\n\n")
cat("| Approach | Model | Input | Sensitivity | Specificity | Balanced accuracy |\n")
cat("|---|---|---|---|---|---|\n")
for (ap in approaches) for (mo in models) for (mi in micro_ord) {
  r <- get_m(mo, ap, mi)
  cat(sprintf("| %s | %s | %s | %s | %s | %s |\n", appr_lab[[ap]], model_lab[[mo]], micro_lab[[mi]],
              pm(r$Sensitivity_mean, r$Sensitivity_sd), pm(r$Specificity_mean, r$Specificity_sd),
              pm(r$Balanced_Accuracy_mean, r$Balanced_Accuracy_sd)))
}

cat("\n== Summary ==\n")
for (mi in micro_ord) {
  x <- metrics[metrics$Microbiome == mi, ]
  cat(sprintf("%-18s mean AUROC of the 6: %s | best: %s | %s (%s) | combinations above 0.5: %d\n",
              micro_lab[[mi]], f3(mean(x$AUROC_mean)), x$Model[which.max(x$AUROC_mean)],
              x$Approach[which.max(x$AUROC_mean)], f3(max(x$AUROC_mean)), sum(x$AUROC_mean > 0.5)))
}
cat(sprintf("Combinations (of 18) with mean AUROC above the sample-count reference (%s): %d\n",
            f3(cnt_m[1]), sum(metrics$AUROC_mean > cnt_m[1])))
cat(sprintf("PR-AUC of a classifier without discrimination = preterm prevalence: %s overall (14/43), %s mean across folds\n",
            f3(mean(S$subject_labels$preterm == "1")), f3(mean(count_fold$n_ptb / count_fold$n_test))))

# -----------------------------------------------------------------------------
# Archivos
# -----------------------------------------------------------------------------
keep <- c("model_name", "approach", "microbiome", "fold_results", "summary",
          "test_predictions", "selected_taxa", "selected_variables", "availability")
saveRDS(lapply(all_res, function(r) r[intersect(keep, names(r))]),
        file.path(outdir, paste0(run_tag, "_cv_folds.rds")))
cat(sprintf("\nSaved the 18 combinations to %s and the pipeline/engine output to %s\n",
            file.path(outdir, paste0(run_tag, "_cv_folds.rds")), log_file))

if (nzchar(results_dir)) {
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(metrics, file.path(results_dir, "metrics_by_combination.csv"), row.names = FALSE)
  utils::write.csv(by_fold, file.path(results_dir, "metrics_by_fold.csv"), row.names = FALSE)
  utils::write.csv(paired, file.path(results_dir, "microbiome_vs_clinical_only.csv"), row.names = FALSE)
  utils::write.csv(count_fold, file.path(results_dir, "sample_count_reference_by_fold.csv"),
                   row.names = FALSE)
  cat("Wrote metrics_by_combination.csv, metrics_by_fold.csv, microbiome_vs_clinical_only.csv and",
      "sample_count_reference_by_fold.csv to", results_dir, "\n")
}
say("Total %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
