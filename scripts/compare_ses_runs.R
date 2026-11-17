#!/usr/bin/env Rscript
# =============================================================================
# compare_ses_runs.R  —  Tablas de docs/SOCIOECONOMIC_INDICATORS.md
# -----------------------------------------------------------------------------
# QUE ES: compara dos corridas de scripts/run_baseline.R del analisis principal
# (clinical_availability = at_visit) que difieren SOLO en como se deriva
# low_education (y por lo tanto ses_risk_score) en el chunk feature_engineering:
#   BEFORE  nivel_academico comparado como texto: marca solo primaria
#   AFTER   R/clinical_features.R lee el codigo: marca primaria y secundaria
# Verifica que las dos usan los mismos folds y sujetos de test, la misma seleccion
# de ANCOM-BC2 y que el Approach 1 (que no usa estas variables) da lo mismo.
# Imprime las tablas en markdown y, si se da PTB_RESULTS_DIR, escribe los CSV.
#
# CORRER (con las dos corridas hechas, ver la seccion "Reproducing" del documento):
#   PTB_BEFORE_TAG=atvisit_2026-10-07 PTB_AFTER_TAG=sesfix_2026-10-16 \
#   PTB_RESULTS_DIR=results/socioeconomic_indicators Rscript scripts/compare_ses_runs.R
# =============================================================================

source(file.path("R", "io.R"))
source(file.path("R", "availability.R"))
cfg  <- read_config(Sys.getenv("PTB_PROFILE", "default"))
dict <- as.data.frame(read_data_dictionary(cfg))
avail <- stats::setNames(dict$availability, dict$variable)
outdir <- cfg$output$dir %||% "analysis/_output"

tags <- c(before = Sys.getenv("PTB_BEFORE_TAG", "atvisit_2026-10-07"),
          after  = Sys.getenv("PTB_AFTER_TAG",  "sesfix_2026-10-16"))
results_dir <- Sys.getenv("PTB_RESULTS_DIR", "")
ses_vars <- c("low_education", "ses_risk_score")

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
stopifnot(length(combos) == 12, setequal(names(runs$after), combos))
ord <- do.call(order, lapply(c("approach", "microbiome", "model_name"),
                             function(f) vapply(runs$before[combos], `[[`, "", f)))
combos <- combos[ord]

# --- 1. checks ---------------------------------------------------------------
cat("== Checks ==\n")
fold_ids <- function(r) {
  tp <- r$test_predictions
  paste(tp$fold, tp$id)[order(tp$fold, tp$id)]
}
same_folds <- all(vapply(combos, function(k)
  identical(fold_ids(runs$before[[k]]), fold_ids(runs$after[[k]])), logical(1)))
cat("Same outer folds and test subjects in both runs:", same_folds, "\n")

ancom <- function(r) if (is.null(r$selected_taxa)) NULL else as.data.frame(r$selected_taxa)
same_taxa <- all(vapply(combos, function(k)
  identical(ancom(runs$before[[k]]), ancom(runs$after[[k]])), logical(1)))
cat("Same ANCOM-BC2 taxa per fold in both runs:", same_taxa, "\n")

a1 <- combos[grepl("Approach1", combos)]
same_a1 <- all(vapply(a1, function(k)
  identical(runs$before[[k]]$fold_results, runs$after[[k]]$fold_results), logical(1)))
cat("Approach 1 fold results identical in both runs:", same_a1, "\n")

same_allowed <- identical(runs$before[[1]]$availability$allowed, runs$after[[1]]$availability$allowed)
cat("Same admitted availability classes:", same_allowed, "(",
    paste(runs$after[[1]]$availability$allowed, collapse = ", "), ")\n")

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
  cat(sprintf("%-7s variables of a class not admitted in a model: %d\n", nm, length(bad)))
}

# --- 2. metrics --------------------------------------------------------------
metric_row <- function(r, run) {
  s <- r$summary
  uses <- intersect(ses_vars, used_vars(r)$variable)
  data.frame(run = run, Model = r$model_name, Approach = r$approach,
             Microbiome = r$microbiome,
             uses_ses_variables = if (length(uses)) paste(uses, collapse = ";") else "",
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
row_of <- function(run, k) metrics[metrics$run == run &
  paste(metrics$Model, metrics$Approach, metrics$Microbiome, sep = "|") == k, ]

cat("\n== Before and after (mean ± SD across the 5 outer folds) ==\n\n")
md_row(c("Model", "Clinical variable set", "Microbiome input", "Uses the indicators",
         "AUROC before", "AUROC after", "ΔAUROC", "PRAUC before", "PRAUC after",
         "Sensitivity before", "Sensitivity after", "Specificity before", "Specificity after"))
md_row(rep("---", 13))
for (k in combos) {
  x <- row_of("before", k); y <- row_of("after", k)
  uses <- unique(c(strsplit(x$uses_ses_variables, ";")[[1]], strsplit(y$uses_ses_variables, ";")[[1]]))
  md_row(c(labs(k), if (length(uses)) paste(sprintf("`%s`", uses), collapse = ", ") else "—",
           pm(x$AUROC_mean, x$AUROC_sd), pm(y$AUROC_mean, y$AUROC_sd),
           sprintf("%+.3f", y$AUROC_mean - x$AUROC_mean),
           pm(x$PRAUC_mean, x$PRAUC_sd), pm(y$PRAUC_mean, y$PRAUC_sd),
           pm(x$Sensitivity_mean, x$Sensitivity_sd), pm(y$Sensitivity_mean, y$Sensitivity_sd),
           pm(x$Specificity_mean, x$Specificity_sd), pm(y$Specificity_mean, y$Specificity_sd)))
}
d <- merge(metrics[metrics$run == "before", ], metrics[metrics$run == "after", ],
           by = c("Model", "Approach", "Microbiome"))
cat(sprintf("\nMean over the 12 combinations: AUROC %.3f → %.3f; PRAUC %.3f → %.3f; ",
            mean(d$AUROC_mean.x), mean(d$AUROC_mean.y), mean(d$PRAUC_mean.x), mean(d$PRAUC_mean.y)))
cat(sprintf("combinations with AUROC > 0.5: %d → %d\n",
            sum(d$AUROC_mean.x > 0.5), sum(d$AUROC_mean.y > 0.5)))
for (a in names(appr_lab)) {
  s <- d[d$Approach == a, ]
  cat(sprintf("%s: mean ΔAUROC %+.3f (range %+.3f to %+.3f); combinations changed: %d of %d\n",
              a, mean(s$AUROC_mean.y - s$AUROC_mean.x), min(s$AUROC_mean.y - s$AUROC_mean.x),
              max(s$AUROC_mean.y - s$AUROC_mean.x),
              sum(s$AUROC_mean.y != s$AUROC_mean.x | s$Sensitivity_mean.y != s$Sensitivity_mean.x |
                    s$Specificity_mean.y != s$Specificity_mean.x), nrow(s)))
}
best <- function(run) {
  m <- metrics[metrics$run == run, ]
  m <- m[order(-m$AUROC_mean, -m$PRAUC_mean), ][1, ]
  sprintf("%s | %s | %s %s", m$Model, m$Approach, m$Microbiome, pm(m$AUROC_mean, m$AUROC_sd))
}
cat("Best before:", best("before"), "\nBest after: ", best("after"), "\n")

# per-fold AUROC of the combinations whose results changed
fold_auc <- function(r) {
  fr <- as.data.frame(r$fold_results)
  stats::setNames(fr$AUROC, fr$fold)
}
cat("\n== Per-fold AUROC of the combinations that changed ==\n\n")
for (k in combos) {
  b <- fold_auc(runs$before[[k]]); a <- fold_auc(runs$after[[k]])
  if (identical(b, a)) next
  cat(sprintf("%-55s %s  →  %s\n", k, paste(f3(b), collapse = " "), paste(f3(a), collapse = " ")))
}

# --- 3. clinical variables ---------------------------------------------------
vars_long <- do.call(rbind, lapply(names(runs), function(nm) {
  do.call(rbind, lapply(combos, function(k) {
    r <- runs[[nm]][[k]]
    u <- used_vars(r)
    if (!nrow(u)) return(NULL)
    data.frame(run = nm, Model = r$model_name, Approach = r$approach,
               Microbiome = r$microbiome, u, availability = unname(avail[u$variable]),
               stringsAsFactors = FALSE)
  }))
}))

cat("\n== Approach 2: final variable list ==\n\n")
a2 <- unique(vars_long[vars_long$Approach == "Approach2_Literature", c("run", "variable")])
md_row(c("Variable", "Before", "After"))
md_row(rep("---", 3))
for (v in unique(a2$variable)) {
  md_row(c(sprintf("`%s`", v), vapply(names(runs), function(nm)
    if (any(a2$run == nm & a2$variable == v)) "✓" else "—", "")))
}
cat("\nList length:", paste(sprintf("%s %d", names(runs),
    vapply(names(runs), function(nm) sum(a2$run == nm), 0L)), collapse = ", "), "\n")

cat("\n== Approach 3: folds (of 5) in which each variable was selected ==\n\n")
a3 <- vars_long[vars_long$Approach == "Approach3_DataDriven", ]
for (nm in names(runs)) {
  s <- split(a3[a3$run == nm, c("fold", "variable")], paste(a3$Model, a3$Microbiome)[a3$run == nm])
  sel  <- lapply(s, function(x) paste(x$fold, x$variable))
  cat(sprintf("%s: selection identical in the four Approach 3 combinations: %s\n", nm,
              all(vapply(sel, identical, logical(1), sel[[1]]))))
}
a3_one <- a3[a3$Model == "glmnet_base" & a3$Microbiome == "ANCOM_Taxa", ]
cnt <- function(nm, v) sum(a3_one$run == nm & a3_one$variable == v)
a3_vars <- unique(a3_one$variable)
a3_vars <- a3_vars[order(-vapply(a3_vars, function(v) cnt("after", v), 0L),
                         -vapply(a3_vars, function(v) cnt("before", v), 0L))]
cat("\n")
md_row(c("Variable", "Before", "After"))
md_row(rep("---", 3))
for (v in a3_vars) md_row(c(sprintf("`%s`", v), cnt("before", v), cnt("after", v)))
for (nm in names(runs)) {
  x <- a3_one[a3_one$run == nm, ]
  cat(sprintf("%s: variables per fold %s; folds with an indicator: %s\n", nm,
              paste(tapply(x$variable, x$fold, length), collapse = ", "),
              paste(sort(unique(x$fold[x$variable %in% ses_vars])), collapse = ", ")))
}
cat("\nPer-fold changes in the Approach 3 selection:\n")
for (f in sort(unique(a3_one$fold))) {
  b <- a3_one$variable[a3_one$run == "before" & a3_one$fold == f]
  a <- a3_one$variable[a3_one$run == "after"  & a3_one$fold == f]
  if (setequal(a, b)) next
  cat(sprintf("  fold %s: out %s; in %s\n", f,
              paste(setdiff(b, a), collapse = ", "), paste(setdiff(a, b), collapse = ", ")))
}

# --- 4. files ----------------------------------------------------------------
if (nzchar(results_dir)) {
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(metrics, file.path(results_dir, "metrics_by_combination.csv"), row.names = FALSE)
  utils::write.csv(vars_long, file.path(results_dir, "clinical_variables_by_fold.csv"), row.names = FALSE)
  cat("\nWritten:", file.path(results_dir, c("metrics_by_combination.csv",
                                              "clinical_variables_by_fold.csv")), sep = "\n  ")
}
