#!/usr/bin/env Rscript
# =============================================================================
# check_variable_availability.R  —  Evidencia de la columna `availability`
# -----------------------------------------------------------------------------
# QUE ES: recalcula desde los datos la evidencia que cita config/data_dictionary.csv
# (columnas availability / availability_source / derived_from) y cuenta en cuantos
# folds entro al modelo cada variable que NO se conoce en la visita. Es la fuente de
# los numeros de docs/CLINICAL_VARIABLE_AVAILABILITY.md. Solo lee; no escribe nada.
#
# QUE IMPRIME:
#   1. Variables clinicas: si varian entre visitas de una misma participante (medidas
#      en la visita) o son constantes (un dato por embarazo), junto a su clase.
#   2. Indicadores 0/1: en cuantas participantes a termino y prematuras valen 1.
#   3. Approach 2: candidatas del config con su clase.
#   4. Approach 3: seleccion clinica por fold de una corrida guardada, con su clase.
#
# REQUISITOS: data/raw/ (matriz de la cohorte, fuera de git) y, para el punto 4, el
# <tag>_cv_folds.rds que escribe scripts/run_baseline.R en analysis/_output/.
#
# CORRER:  PTB_RUN_TAG=thresholdfix_2026-09-30 Rscript scripts/check_variable_availability.R
# =============================================================================

source(file.path("R", "io.R"))
cfg  <- read_config(Sys.getenv("PTB_PROFILE", "default"))
tag  <- Sys.getenv("PTB_RUN_TAG", "thresholdfix_2026-09-30")
dict <- as.data.frame(read_data_dictionary(cfg))
m    <- as.data.frame(load_dataset(cfg)$matrix)

id <- cfg$columns$subject_id; out <- cfg$columns$outcome
avail <- stats::setNames(dict$availability, dict$variable)
subj_out <- tapply(m[[out]], m[[id]], `[`, 1)

# --- 1. variacion entre visitas ------------------------------------------------
clin <- dict$variable[dict$role == "clinical"]
ev <- do.call(rbind, lapply(clin, function(v) {
  n_obs <- tapply(m[[v]], m[[id]], function(z) sum(!is.na(z)))
  n_val <- tapply(m[[v]], m[[id]], function(z) length(unique(z[!is.na(z)])))
  data.frame(variable = v, availability = avail[[v]],
             multi_sample = sum(n_obs >= 2), varying = sum(n_obs >= 2 & n_val > 1))
}))
cat("\n== 1. Variacion entre visitas (participantes con 2 o mas muestras) ==\n")
print(ev[order(ev$availability, ev$varying == 0), ], row.names = FALSE)

# --- 2. indicadores 0/1 por desenlace ------------------------------------------
is01 <- vapply(clin, function(v) is.numeric(m[[v]]) && all(m[[v]] %in% c(0, 1, NA)), logical(1))
fl <- do.call(rbind, lapply(clin[is01], function(v) {
  s <- tapply(m[[v]], m[[id]], function(z) if (all(is.na(z))) NA else max(z, na.rm = TRUE))
  data.frame(variable = v, availability = avail[[v]],
             term_1 = sum(s == 1 & subj_out == 0, na.rm = TRUE), term_n = sum(subj_out == 0),
             preterm_1 = sum(s == 1 & subj_out == 1, na.rm = TRUE), preterm_n = sum(subj_out == 1))
}))
cat("\n== 2. Indicadores 0/1: participantes con valor 1 ==\n")
print(fl, row.names = FALSE)

# --- 3. Approach 2: candidatas -------------------------------------------------
a2 <- cfg$approaches$approach2_literature
a2_vars <- c(names(a2$core_vars), names(a2$candidates))
cat("\n== 3. Approach 2: candidatas del config ==\n")
print(data.frame(variable = a2_vars, availability = unname(avail[a2_vars])), row.names = FALSE)
cat("(La lista final del Approach 2 es fija en los 5 folds; la imprime el chunk approach2_setup.)\n")

# --- 4. Approach 3: seleccion por fold -----------------------------------------
rds <- file.path(cfg$output$dir %||% "analysis/_output", paste0(tag, "_cv_folds.rds"))
if (!file.exists(rds)) {
  cat("\n== 4. Sin", rds, "- correr scripts/run_baseline.R con PTB_RUN_TAG =", tag, "==\n")
} else {
  res <- readRDS(rds)
  a3  <- Filter(function(r) !is.null(r$selected_variables), res)
  sel <- lapply(a3, function(r) r$selected_variables[, c("fold", "variable")])
  same <- all(vapply(sel, function(s) identical(s, sel[[1]]), logical(1)))
  s <- sel[[1]]
  tab <- as.data.frame(table(variable = s$variable), stringsAsFactors = FALSE)
  names(tab)[2] <- "folds"
  tab$availability <- unname(avail[tab$variable])
  cat(sprintf("\n== 4. Approach 3 (%s): %d combinaciones, seleccion identica entre ellas: %s ==\n",
              tag, length(a3), same))
  print(tab[order(tab$availability, -tab$folds), ], row.names = FALSE)
  cat("Variables por fold:", paste(tapply(s$variable, s$fold, length), collapse = ", "), "\n")
  late <- tapply(avail[s$variable] != "at_visit", s$fold, sum)
  cat("Variables no disponibles en la visita, por fold:", paste(late, collapse = ", "), "\n")
}
