#!/usr/bin/env Rscript
# =============================================================================
# sample_count_reference.R  —  Tablas de docs/SAMPLE_COUNT_REFERENCE.md
# -----------------------------------------------------------------------------
# QUE ES: una referencia para el AUROC de los modelos. Una participante con parto
# prematuro deja de dar muestras antes (el embarazo termina antes), asi que el NUMERO
# de muestras por participante ya discrimina el desenlace sin mirar el microbioma ni
# la clinica. Se evalua como un score fijo, sin nada ajustado:
#   score = -(numero de muestras)   -> menos muestras = mas probable prematuro,
# orientacion fijada de antemano y AUROC con roc_ptb() del motor (direction "<"),
# en los MISMOS 5 folds externos que los modelos, fold a fold (media ± SD, como el
# motor). Variantes:
#   all            todas las muestras, todas las participantes
#   all_reached_X  todas las muestras, solo las participantes cuyo parto fue en la
#                  semana X o despues (separa el efecto de sacar participantes)
#   before_X       solo las muestras tomadas antes de la semana X, entre esas mismas
#                  participantes: todas seguian embarazadas en toda la ventana, asi
#                  que la truncacion por el parto no entra en el conteo
# X = PTB_COUNT_WEEKS (por defecto 24, 28, 32). sdg_visita / sdg_parto estan en
# notacion semanas.dias (24.4 = 24 semanas y 4 dias): con X entero, "< X" y ">= X"
# valen igual en esa notacion que en semanas decimales.
#
# CHEQUEOS (abortan si fallan):
#   1. los folds se reconstruyen con el mismo codigo que el chunk nested_cv_setup del
#      .Rmd (semilla y v del config) y tienen que coincidir sujeto por sujeto con los
#      de una corrida guardada (PTB_FOLDS_RUN_TAG, por defecto la linea base
#      atvisit_2026-10-07) si su _cv_folds.rds existe;
#   2. el mismo calculo por fold, aplicado a las predicciones guardadas de esa corrida,
#      reproduce el AUROC por fold que reporto el motor en las 12 combinaciones.
#
# CORRER (necesita los datos de la cohorte en data/raw/):
#   PTB_RESULTS_DIR=results/sample_count_reference Rscript scripts/sample_count_reference.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(rsample)
})
source(file.path("R", "io.R"))
source(file.path("R", "nested_cv.R"))   # roc_ptb(): la misma metrica que el motor

cfg <- read_config(Sys.getenv("PTB_PROFILE", "default"))
outdir      <- cfg$output$dir %||% "analysis/_output"
results_dir <- Sys.getenv("PTB_RESULTS_DIR", "")
run_tag     <- Sys.getenv("PTB_FOLDS_RUN_TAG", "atvisit_2026-10-07")
model_csv   <- Sys.getenv("PTB_MODEL_METRICS",
                          file.path("results", "clinical_availability", "metrics_by_combination.csv"))
weeks <- as.numeric(strsplit(Sys.getenv("PTB_COUNT_WEEKS", "24,28,32"), ",")[[1]])
stopifnot(length(weeks) > 0, all(!is.na(weeks)), all(weeks == round(weeks)))

sid   <- cfg$columns$subject_id
y     <- cfg$columns$outcome
pos   <- cfg$columns$outcome_positive
ga_d  <- cfg$columns$gestational_age   # edad gestacional al parto
ga_v  <- cfg$columns$visit_ga          # edad gestacional en la toma de muestra

samples <- as.data.frame(load_dataset(cfg)$matrix)
samples <- data.frame(id = as.character(samples[[sid]]),
                      preterm = ifelse(as.character(samples[[y]]) == pos, "1", "0"),
                      ga_delivery = samples[[ga_d]], ga_visit = samples[[ga_v]],
                      stringsAsFactors = FALSE)
stopifnot(!anyNA(samples))
# una edad al parto y un desenlace por participante
stopifnot(all(tapply(samples$ga_delivery, samples$id, function(v) length(unique(v))) == 1),
          all(tapply(samples$preterm, samples$id, function(v) length(unique(v))) == 1))

# --- 1. outer folds, as in the chunk nested_cv_setup of the .Rmd ---------------
set.seed(cfg$cv$seed)
subject_labels <- samples %>%
  group_by(id) %>%
  summarize(preterm = ifelse(any(preterm == "1"), "1", "0")) %>%
  ungroup()
cv_folds <- vfold_cv(subject_labels, v = cfg$cv$outer_folds, strata = preterm)
folds <- bind_rows(lapply(seq_len(nrow(cv_folds)), function(i)
  data.frame(id = assessment(cv_folds$splits[[i]])$id, fold = i, stringsAsFactors = FALSE)))
stopifnot(!anyDuplicated(folds$id), setequal(folds$id, subject_labels$id))

cat("== Checks ==\n")
run_file <- file.path(outdir, paste0(run_tag, "_cv_folds.rds"))
run <- NULL
if (file.exists(run_file)) {
  run <- readRDS(run_file)
  key <- function(r) paste(r$model_name, r$approach, r$microbiome, sep = "|")
  names(run) <- vapply(run, key, "")
  fold_ids <- function(tp) sort(paste(tp$fold, tp$id))
  same <- vapply(run, function(r) identical(fold_ids(r$test_predictions),
                                            fold_ids(folds)), logical(1))
  if (!all(same)) stop("Rebuilt folds differ from the folds of run ", run_tag, ": ",
                       paste(names(run)[!same], collapse = ", "))
  cat(sprintf("Rebuilt outer folds identical to the %d combinations of run %s: TRUE\n",
              length(run), run_tag))
} else {
  cat(sprintf("Run %s not found in %s: fold check against the models SKIPPED\n",
              run_tag, outdir))
}

# AUROC per outer fold of a subject-level score (higher = preterm), as the engine does
fold_auroc <- function(df) {
  df %>%
    group_by(fold) %>%
    summarise(n_participants = n(), n_preterm = sum(preterm == "1"),
              n_term = sum(preterm == "0"), n_samples = sum(n_samples),
              AUROC = if (n_preterm > 0 && n_term > 0)
                as.numeric(pROC::auc(roc_ptb(preterm, score))) else NA_real_,
              .groups = "drop") %>%
    arrange(fold)
}

if (!is.null(run)) {
  engine_ok <- vapply(run, function(r) {
    tp <- r$test_predictions
    mine <- fold_auroc(data.frame(fold = tp$fold, preterm = as.character(tp$true_class),
                                  score = tp$pred_prob, n_samples = 0))
    eng <- r$fold_results[order(r$fold_results$fold), ]
    identical(mine$fold, as.integer(eng$fold)) &&
      identical(mine$AUROC, as.numeric(eng$AUROC))
  }, logical(1))
  if (!all(engine_ok)) stop("The per-fold AUROC does not reproduce the engine in: ",
                            paste(names(run)[!engine_ok], collapse = ", "))
  cat(sprintf("Per-fold AUROC of the saved predictions reproduces the engine in %d/%d combinations: TRUE\n",
              sum(engine_ok), length(run)))
}

# --- 2. variants ---------------------------------------------------------------
participants <- samples %>%
  group_by(id) %>%
  summarise(preterm = first(preterm), ga_delivery = first(ga_delivery), .groups = "drop") %>%
  left_join(folds, by = "id")

count_variant <- function(variant, week, keep_ids, window_week) {
  s <- samples[samples$id %in% keep_ids, ]
  inwin <- if (is.na(window_week)) rep(TRUE, nrow(s)) else s$ga_visit < window_week
  n <- tapply(inwin, factor(s$id, levels = keep_ids), sum)
  participants[match(keep_ids, participants$id), c("id", "preterm", "fold")] %>%
    mutate(variant = variant, week = week, n_samples = as.integer(n[id]), score = -n_samples)
}

variants <- list(count_variant("all", NA, participants$id, NA))
for (X in weeks) {
  reached <- participants$id[participants$ga_delivery >= X]
  variants <- c(variants, list(count_variant(paste0("all_reached_", X), X, reached, NA),
                               count_variant(paste0("before_", X), X, reached, X)))
}
scores <- bind_rows(variants)
labels <- c(all = "All samples, all participants")
short  <- c(all = "All")
for (X in weeks) {
  labels[paste0("all_reached_", X)] <- sprintf("All samples, participants who reached week %d", X)
  labels[paste0("before_", X)]      <- sprintf("Samples before week %d, participants who reached week %d", X, X)
  short[paste0("all_reached_", X)]  <- sprintf("All, reached %d", X)
  short[paste0("before_", X)]       <- sprintf("Before %d", X)
}
scores$variant <- factor(scores$variant, levels = names(labels))

by_fold <- scores %>%
  group_by(variant, week) %>%
  group_modify(~ fold_auroc(.x)) %>%
  ungroup()

summary_tab <- scores %>%
  group_by(variant, week) %>%
  group_modify(function(d, k) {
    r  <- roc_ptb(d$preterm, d$score)
    ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
    f  <- fold_auroc(d)
    data.frame(n_participants = nrow(d), n_preterm = sum(d$preterm == "1"),
               n_term = sum(d$preterm == "0"), n_samples = sum(d$n_samples),
               mean_samples_preterm = mean(d$n_samples[d$preterm == "1"]),
               mean_samples_term = mean(d$n_samples[d$preterm == "0"]),
               participants_no_sample = sum(d$n_samples == 0),
               AUROC_mean = mean(f$AUROC), AUROC_sd = stats::sd(f$AUROC),
               n_folds = sum(!is.na(f$AUROC)),
               AUROC_pooled = as.numeric(pROC::auc(r)),
               AUROC_pooled_ci_low = ci[1], AUROC_pooled_ci_high = ci[3])
  }) %>%
  ungroup()

distribution <- scores %>%
  count(variant, week, n_samples, preterm) %>%
  tidyr::pivot_wider(names_from = preterm, values_from = n, values_fill = 0) %>%
  rename(term = `0`, preterm = `1`) %>%
  arrange(variant, n_samples)

# --- 3. tables -----------------------------------------------------------------
f3 <- function(x) sprintf("%.3f", x)
cat("\n== Summary ==\n\n")
cat("| Variant | Participants (preterm / term) | Samples | Mean samples per participant (preterm / term) | Participants with no sample in the window | AUROC, mean ± SD across folds | Pooled AUROC (95% CI) |\n")
cat("|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(summary_tab))) with(summary_tab[i, ], cat(sprintf(
  "| %s | %d (%d / %d) | %d | %.2f / %.2f | %d | %s ± %s | %s (%.2f–%.2f) |\n",
  labels[[as.character(variant)]], n_participants, n_preterm, n_term, n_samples,
  mean_samples_preterm, mean_samples_term, participants_no_sample,
  f3(AUROC_mean), f3(AUROC_sd), f3(AUROC_pooled), AUROC_pooled_ci_low, AUROC_pooled_ci_high)))

cat("\n== AUROC by outer fold ==\n\n")
vlev <- levels(scores$variant)
cat("| Fold | Participants (preterm) |", paste(short[vlev], collapse = " | "), "|\n")
cat("|---|---|", paste(rep("---", length(vlev)), collapse = "|"), "|\n", sep = "")
for (k in sort(unique(by_fold$fold))) {
  bf <- by_fold[by_fold$fold == k, ]
  np <- paste(unique(sprintf("%d (%d)", bf$n_participants, bf$n_preterm)), collapse = " / ")
  cat(sprintf("| %d | %s | %s |\n", k, np,
              paste(f3(bf$AUROC[match(vlev, bf$variant)]), collapse = " | ")))
}

cat("\n== Samples per participant ==\n\n")
print(as.data.frame(distribution), row.names = FALSE)

# Margin above 0.5 kept when the count is restricted to the window before week X,
# among the same participants (before_X vs all_reached_X)
cat("\n== Margin above 0.5 kept within the window ==\n")
for (X in weeks) {
  a <- summary_tab[summary_tab$variant == paste0("all_reached_", X), ]
  b <- summary_tab[summary_tab$variant == paste0("before_", X), ]
  cat(sprintf("Week %d: per-fold mean %s -> %s (%.0f%% of the margin); pooled %s -> %s (%.0f%%)\n",
              X, f3(a$AUROC_mean), f3(b$AUROC_mean), 100 * (b$AUROC_mean - 0.5) / (a$AUROC_mean - 0.5),
              f3(a$AUROC_pooled), f3(b$AUROC_pooled), 100 * (b$AUROC_pooled - 0.5) / (a$AUROC_pooled - 0.5)))
}

# Descriptive: gestational age at the first sample. Weeks.days notation, so only
# order statistics are meaningful: quantile type 1 returns an observed value.
cat("\n== Gestational age at the first sample (weeks.days), median by outcome ==\n")
first_ga <- samples %>% group_by(id) %>%
  summarise(preterm = first(preterm), ga_delivery = first(ga_delivery),
            first_visit = min(ga_visit), .groups = "drop")
for (X in c(NA, weeks)) {
  fg <- if (is.na(X)) first_ga else first_ga[first_ga$ga_delivery >= X, ]
  med <- tapply(fg$first_visit, fg$preterm, stats::quantile, probs = 0.5, type = 1)
  cat(sprintf("  %-28s preterm %.1f (n = %d), term %.1f (n = %d)\n",
              if (is.na(X)) "all participants" else sprintf("participants who reached %d", X),
              med[["1"]], sum(fg$preterm == "1"), med[["0"]], sum(fg$preterm == "0")))
}

# --- 4. comparison with the models ---------------------------------------------
if (file.exists(model_csv)) {
  m <- utils::read.csv(model_csv, stringsAsFactors = FALSE)
  if ("analysis" %in% names(m)) m <- m[m$analysis == "main", ]
  m <- m[order(-m$AUROC_mean), ]
  cat(sprintf("\n== Models (%s, %d combinations) ==\n", model_csv, nrow(m)))
  cat(sprintf("Best: %s | %s | %s  AUROC %s ± %s; mean of the %d combinations %s\n",
              m$Model[1], m$Approach[1], m$Microbiome[1], f3(m$AUROC_mean[1]),
              f3(m$AUROC_sd[1]), nrow(m), f3(mean(m$AUROC_mean))))
  for (i in seq_len(nrow(summary_tab)))
    cat(sprintf("  combinations with mean AUROC above %-24s (%s): %d\n",
                as.character(summary_tab$variant[i]), f3(summary_tab$AUROC_mean[i]),
                sum(m$AUROC_mean > summary_tab$AUROC_mean[i])))
  if (!is.null(run)) {
    best <- run[[paste(m$Model[1], m$Approach[1], m$Microbiome[1], sep = "|")]]
    if (!is.null(best)) {
      fr <- best$fold_results[order(best$fold_results$fold), ]
      ref <- by_fold[by_fold$variant == "all", ]
      cat("  per fold, best model vs all samples:",
          paste(sprintf("%d: %s vs %s", fr$fold, f3(fr$AUROC), f3(ref$AUROC)), collapse = "; "), "\n")
    }
  }
} else {
  cat("\nModel metrics", model_csv, "not found: comparison with the models skipped\n")
}

if (nzchar(results_dir)) {
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(summary_tab, file.path(results_dir, "auroc_summary.csv"), row.names = FALSE)
  utils::write.csv(by_fold, file.path(results_dir, "auroc_by_fold.csv"), row.names = FALSE)
  utils::write.csv(distribution, file.path(results_dir, "samples_per_participant.csv"), row.names = FALSE)
  cat("\nWrote auroc_summary.csv, auroc_by_fold.csv and samples_per_participant.csv to",
      results_dir, "\n")
}
