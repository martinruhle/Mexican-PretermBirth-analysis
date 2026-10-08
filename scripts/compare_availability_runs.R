#!/usr/bin/env Rscript
# =============================================================================
# compare_availability_runs.R  —  Tablas de docs/CLINICAL_AVAILABILITY_RESULTS.md
# -----------------------------------------------------------------------------
# QUE ES: compara tres corridas de scripts/run_baseline.R que difieren SOLO en las
# clases de disponibilidad admitidas como predictores clinicos (columna `availability`
# de config/data_dictionary.csv):
#   BEFORE  todas las clases (at_visit, after_visit, outcome_defined): reproduce el
#           analisis anterior a la clasificacion
#   MAIN    analisis principal, solo at_visit (config clinical_availability$allowed)
#   SENS    analisis de sensibilidad, at_visit + after_visit
# Verifica que las tres usan los mismos folds y la misma seleccion de ANCOM-BC2, que el
# Approach 1 (dos variables at_visit) da lo mismo en las tres, y que ninguna variable de
# una clase no admitida entro a un modelo. Imprime las tablas en markdown y, si se da
# PTB_RESULTS_DIR, escribe los CSV.
#
# CORRER (con las tres corridas hechas, ver la seccion "Reproducing" del documento):
#   PTB_BEFORE_TAG=availctrl_2026-10-07 PTB_MAIN_TAG=atvisit_2026-10-07 \
#   PTB_SENS_TAG=aftervisit_2026-10-07 PTB_RESULTS_DIR=results/clinical_availability \
#     Rscript scripts/compare_availability_runs.R
# =============================================================================

source(file.path("R", "io.R"))
source(file.path("R", "availability.R"))
cfg  <- read_config(Sys.getenv("PTB_PROFILE", "default"))
dict <- as.data.frame(read_data_dictionary(cfg))
avail <- stats::setNames(dict$availability, dict$variable)
outdir <- cfg$output$dir %||% "analysis/_output"

tags <- c(before      = Sys.getenv("PTB_BEFORE_TAG", "availctrl_2026-10-07"),
          main        = Sys.getenv("PTB_MAIN_TAG",   "atvisit_2026-10-07"),
          sensitivity = Sys.getenv("PTB_SENS_TAG",   "aftervisit_2026-10-07"))
results_dir <- Sys.getenv("PTB_RESULTS_DIR", "")
# Optional: a run made before the engine read the availability column. BEFORE must
# reproduce it exactly (same code path with every class admitted).
reference_tag <- Sys.getenv("PTB_REFERENCE_TAG", "")

model_lab <- c(rf_base = "Random Forest", glmnet_base = "Elastic net (glmnet)")
appr_lab  <- c(Approach1_DREAM = "1 — DREAM (minimal)", Approach2_Literature = "2 — literature",
               Approach3_DataDriven = "3 — data-driven")
micro_lab <- c(ANCOM_Taxa = "ANCOM-BC2 taxa", Full_Microbiome = "Full microbiome")
key <- function(r) paste(r$model_name, r$approach, r$microbiome, sep = "|")

runs <- lapply(tags, function(t) {
  r <- readRDS(file.path(outdir, paste0(t, "_cv_folds.rds")))
  stats::setNames(r, vapply(r, key, ""))
})
combos <- names(runs$before)
stopifnot(all(vapply(runs, function(r) setequal(names(r), combos), logical(1))))
# fixed order: approach, microbiome, model
ord <- do.call(order, lapply(c("approach", "microbiome", "model_name"),
                             function(f) vapply(runs$before[combos], `[[`, "", f)))
combos <- combos[ord]

# --- 1. checks ---------------------------------------------------------------
cat("== Checks ==\n")
if (nzchar(reference_tag)) {
  ref <- readRDS(file.path(outdir, paste0(reference_tag, "_cv_folds.rds")))
  ref <- stats::setNames(ref, vapply(ref, key, ""))
  parts <- c("fold_results", "summary", "test_predictions", "selected_taxa", "selected_variables")
  same_ref <- setequal(names(ref), combos) && all(vapply(combos, function(k)
    identical(ref[[k]][parts], runs$before[[k]][parts]), logical(1)))
  cat(sprintf("BEFORE (%s) identical to the reference run %s: %s\n",
              tags[["before"]], reference_tag, same_ref))
}
fold_ids <- function(r) {
  tp <- r$test_predictions
  paste(tp$fold, tp$id)[order(tp$fold, tp$id)]
}
same_folds <- all(vapply(combos, function(k)
  identical(fold_ids(runs$before[[k]]), fold_ids(runs$main[[k]])) &&
  identical(fold_ids(runs$before[[k]]), fold_ids(runs$sensitivity[[k]])), logical(1)))
cat("Same outer folds and test subjects in the three runs:", same_folds, "\n")

ancom <- function(r) if (is.null(r$selected_taxa)) NULL else as.data.frame(r$selected_taxa)
same_taxa <- all(vapply(combos, function(k)
  identical(ancom(runs$before[[k]]), ancom(runs$main[[k]])) &&
  identical(ancom(runs$before[[k]]), ancom(runs$sensitivity[[k]])), logical(1)))
cat("Same ANCOM-BC2 taxa per fold in the three runs:", same_taxa, "\n")

a1 <- combos[grepl("Approach1", combos)]
same_a1 <- all(vapply(a1, function(k)
  identical(runs$before[[k]]$fold_results, runs$main[[k]]$fold_results) &&
  identical(runs$before[[k]]$fold_results, runs$sensitivity[[k]]$fold_results), logical(1)))
cat("Approach 1 fold results identical in the three runs:", same_a1, "\n")

# variables that reached a model: Approach 1/2 fixed list; Approach 3 per-fold selection
used_vars <- function(r) {
  if (!is.null(r$selected_variables)) {
    data.frame(fold = as.character(r$selected_variables$fold),
               variable = r$selected_variables$variable, stringsAsFactors = FALSE)
  } else {
    data.frame(fold = "all", variable = r$availability$clinical_variables,
               stringsAsFactors = FALSE)
  }
}
for (nm in names(runs)) {
  allowed <- runs[[nm]][[1]]$availability$allowed
  bad <- unlist(lapply(runs[[nm]], function(r) {
    v <- used_vars(r)$variable
    v[!avail[v] %in% allowed]
  }))
  cat(sprintf("%-11s admitted: %-40s variables of another class in a model: %d\n",
              nm, paste(allowed, collapse = ", "), length(bad)))
}

# --- 2. metrics --------------------------------------------------------------
metric_row <- function(r, analysis) {
  s <- r$summary
  data.frame(analysis = analysis, Model = r$model_name, Approach = r$approach,
             Microbiome = r$microbiome,
             AUROC_mean = s$AUROC_mean, AUROC_sd = s$AUROC_sd,
             PRAUC_mean = s$PRAUC_mean, PRAUC_sd = s$PRAUC_sd,
             Sensitivity_mean = s$Sensitivity_mean, Sensitivity_sd = s$Sensitivity_sd,
             Specificity_mean = s$Specificity_mean, Specificity_sd = s$Specificity_sd,
             Balanced_Accuracy_mean = s$Balanced_Accuracy_mean,
             Balanced_Accuracy_sd = s$Balanced_Accuracy_sd, stringsAsFactors = FALSE)
}
metrics <- do.call(rbind, lapply(names(runs), function(nm)
  do.call(rbind, lapply(combos, function(k) metric_row(runs[[nm]][[k]], nm)))))

f3 <- function(x) sprintf("%.3f", x)
pm <- function(m, s) sprintf("%.3f ± %.3f", m, s)
labs <- function(k) {
  r <- runs$before[[k]]
  c(model_lab[[r$model_name]], appr_lab[[r$approach]], micro_lab[[r$microbiome]])
}
md_row <- function(cells) cat("|", paste(cells, collapse = " | "), "|\n")

headline <- function(nm) {
  m <- metrics[metrics$analysis == nm, ]
  m <- m[order(-m$AUROC_mean, -m$PRAUC_mean), ]
  cat(sprintf("\n== Results: %s (mean ± SD across the 5 outer folds) ==\n\n", nm))
  md_row(c("Model", "Clinical variable set", "Microbiome input", "AUROC", "PRAUC",
           "Sensitivity", "Specificity"))
  md_row(rep("---", 7))
  for (i in seq_len(nrow(m))) {
    md_row(c(model_lab[[m$Model[i]]], appr_lab[[m$Approach[i]]], micro_lab[[m$Microbiome[i]]],
             pm(m$AUROC_mean[i], m$AUROC_sd[i]), pm(m$PRAUC_mean[i], m$PRAUC_sd[i]),
             pm(m$Sensitivity_mean[i], m$Sensitivity_sd[i]),
             pm(m$Specificity_mean[i], m$Specificity_sd[i])))
  }
}
headline("main")
headline("sensitivity")

compare <- function(a, b, title) {
  cat(sprintf("\n== %s: %s → %s (means) ==\n\n", title, a, b))
  md_row(c("Model", "Clinical variable set", "Microbiome input", "AUROC", "ΔAUROC", "PRAUC",
           "Sensitivity", "Specificity"))
  md_row(rep("---", 8))
  for (k in combos) {
    x <- metrics[metrics$analysis == a & paste(metrics$Model, metrics$Approach,
                                                metrics$Microbiome, sep = "|") == k, ]
    y <- metrics[metrics$analysis == b & paste(metrics$Model, metrics$Approach,
                                                metrics$Microbiome, sep = "|") == k, ]
    arrow <- function(col) sprintf("%s → %s", f3(x[[col]]), f3(y[[col]]))
    md_row(c(labs(k), arrow("AUROC_mean"), sprintf("%+.3f", y$AUROC_mean - x$AUROC_mean),
             arrow("PRAUC_mean"), arrow("Sensitivity_mean"), arrow("Specificity_mean")))
  }
  d <- merge(metrics[metrics$analysis == a, ], metrics[metrics$analysis == b, ],
             by = c("Model", "Approach", "Microbiome"))
  cat(sprintf("\nMean over the 12 combinations: AUROC %.3f → %.3f; PRAUC %.3f → %.3f; ",
              mean(d$AUROC_mean.x), mean(d$AUROC_mean.y),
              mean(d$PRAUC_mean.x), mean(d$PRAUC_mean.y)))
  cat(sprintf("combinations with AUROC > 0.5: %d → %d\n",
              sum(d$AUROC_mean.x > 0.5), sum(d$AUROC_mean.y > 0.5)))
}
compare("before", "main", "Before and after")
compare("main", "sensitivity", "Sensitivity analysis")

# --- 3. clinical variables ---------------------------------------------------
vars_long <- do.call(rbind, lapply(names(runs), function(nm) {
  do.call(rbind, lapply(combos, function(k) {
    r <- runs[[nm]][[k]]
    u <- used_vars(r)
    if (!nrow(u)) return(NULL)
    data.frame(analysis = nm, Model = r$model_name, Approach = r$approach,
               Microbiome = r$microbiome, u, availability = unname(avail[u$variable]),
               stringsAsFactors = FALSE)
  }))
}))

# Approach 2: one fixed list per run (the same in its four combinations)
cat("\n== Approach 2: final variable list ==\n\n")
a2 <- unique(vars_long[vars_long$Approach == "Approach2_Literature",
                       c("analysis", "variable", "availability")])
a2_vars <- unique(a2$variable)
md_row(c("Variable", "Class", "Before", "Main", "Sensitivity"))
md_row(rep("---", 5))
for (v in a2_vars) {
  inn <- vapply(names(runs), function(nm) if (any(a2$analysis == nm & a2$variable == v)) "✓" else "—", "")
  md_row(c(sprintf("`%s`", v), avail[[v]], inn))
}
cat("\nList length:", paste(sprintf("%s %d", names(runs),
    vapply(names(runs), function(nm) sum(a2$analysis == nm), 0L)), collapse = ", "), "\n")

# Approach 3: per-fold selection (identical across its four combinations?)
cat("\n== Approach 3: folds (of 5) in which each variable was selected ==\n\n")
a3 <- vars_long[vars_long$Approach == "Approach3_DataDriven", ]
for (nm in names(runs)) {
  s <- split(a3[a3$analysis == nm, c("fold", "variable")],
             paste(a3$Model, a3$Microbiome)[a3$analysis == nm])
  sel  <- lapply(s, function(x) paste(x$fold, x$variable))
  same <- all(vapply(sel, identical, logical(1), sel[[1]]))
  cat(sprintf("%s: selection identical in the four Approach 3 combinations: %s\n", nm, same))
}
a3_one <- a3[a3$Model == "glmnet_base" & a3$Microbiome == "ANCOM_Taxa", ]
a3_vars <- unique(a3_one$variable)
cnt <- function(nm, v) sum(a3_one$analysis == nm & a3_one$variable == v)
a3_vars <- a3_vars[order(match(avail[a3_vars], availability_levels()),
                         -vapply(a3_vars, function(v) cnt("main", v), 0L),
                         -vapply(a3_vars, function(v) cnt("before", v), 0L))]
cat("\n")
md_row(c("Variable", "Class", "Before", "Main", "Sensitivity"))
md_row(rep("---", 5))
for (v in a3_vars) {
  md_row(c(sprintf("`%s`", v), avail[[v]], vapply(names(runs), function(nm) as.character(cnt(nm, v)), "")))
}
for (nm in names(runs)) {
  per_fold <- tapply(a3_one$variable[a3_one$analysis == nm], a3_one$fold[a3_one$analysis == nm], length)
  cat(sprintf("%s: variables per fold %s\n", nm, paste(per_fold, collapse = ", ")))
}

# Approach 3 pool exclusions per run
cat("\n== Approach 3 pool: variables excluded by class ==\n")
for (nm in names(runs)) {
  ex <- runs[[nm]][[combos[grepl("Approach3", combos)][1]]]$availability$excluded
  cat(sprintf("%s (%d): %s\n", nm, nrow(ex),
              if (nrow(ex)) paste(sprintf("%s [%s]", ex$variable, ex$availability), collapse = ", ") else "none"))
}

# --- 4. files ----------------------------------------------------------------
if (nzchar(results_dir)) {
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(metrics, file.path(results_dir, "metrics_by_combination.csv"), row.names = FALSE)
  utils::write.csv(vars_long, file.path(results_dir, "clinical_variables_by_fold.csv"), row.names = FALSE)
  cat("\nWritten:", file.path(results_dir, c("metrics_by_combination.csv",
                                              "clinical_variables_by_fold.csv")), sep = "\n  ")
}
