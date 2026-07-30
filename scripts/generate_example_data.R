#!/usr/bin/env Rscript
# =============================================================================
# generate_example_data.R  —  Datos SINTETICOS de ejemplo (Chat 5)
# -----------------------------------------------------------------------------
# QUE ES: genera un dataset 100% sintetico (NINGUN dato de paciente real) que
# permite correr el pipeline COMPLETO de punta a punta via el perfil de config
# `example:`, sin acceso a data/raw/. Criterio de exito = las 12 combinaciones y
# todas las figuras corren SIN ERROR. NO reproduce el baseline real (datos
# distintos); solo debe correr limpio.
#
# COMO SE DERIVA: la lista de columnas NO se hardcodea — se lee de
# config/data_dictionary.csv (fuente de verdad del contrato):
#   * 97 taxa  = role == "microbiome"  (nombres EXACTOS, incl. Escherichia-Shigella,
#                prefijos f__/o__/c__/d__ y el duplicado f__Bifidobacteriaceae.1)
#   * 70 no-taxa = el resto (id, index, preterm, sdg_parto + clinicas)
#
# ENTREGABLES (data/example/):
#   example_genus_rel.csv       matriz ancha: no-taxa + 97 taxa (abundancias rel, filas ~1)
#   example_genus_abs.csv       misma forma pero taxa = CONTEOS enteros (ANCOM-BC2)
#   example_metadata_long.csv   id, index, preterm, sdg_parto, sdg_visita (figuras poblacion)
#
# DISENO (respeta invariantes del contrato):
#   - n=42 sujetos, 2-3 muestras longitudinales c/u (~100 muestras), ~33% preterm.
#   - preterm (0/1) y sdg_parto son a nivel SUJETO (constantes en sus muestras);
#     sdg_visita varia por muestra y cubre los 3 trimestres.
#   - rel = counts / rowSums(counts): rel suma 1 por muestra y ambos archivos son
#     consistentes (mismas columnas no-taxa, mismos 97 taxa, mismos index/id).
#   - Microbioma vaginal realista: dominado por Lactobacillus, disperso (~77% ceros),
#     con senal debil inyectada en Gardnerella/Prevotella/Sneathia (dysbiosis mas
#     frecuente en preterm). Ningun taxon queda all-zero (rompe cmultRepl).
#   - Senal clinica debil edad/imc/complicaciones <-> preterm (approach2/3 no planos).
#   - Completitud 100% (>= el 80% que exige el contrato).
#
# CORRER:  Rscript scripts/generate_example_data.R
# =============================================================================

suppressWarnings(suppressMessages({
  library(readr)
}))

set.seed(20260730)   # semilla fija -> generacion reproducible

# --- Localizar el root del proyecto (sin setwd, sin rutas absolutas) ---------
root <- tryCatch(here::here(), error = function(e) getwd())
dict_path <- file.path(root, "config", "data_dictionary.csv")
out_dir   <- file.path(root, "data", "example")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 0. LEER EL CONTRATO: columnas taxa vs no-taxa (por NOMBRE, del diccionario)
# =============================================================================
dict <- utils::read.csv(dict_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("variable", "role") %in% names(dict)))

taxa    <- dict$variable[dict$role == "microbiome"]     # 97 taxa (orden del diccionario)
nontaxa <- dict$variable[dict$role != "microbiome"]     # 70 no-taxa
n_taxa  <- length(taxa)
stopifnot(n_taxa == 97, length(nontaxa) == 70)
stopifnot(anyDuplicated(taxa) == 0)                     # incl. el duplicado ...".1"
stopifnot("Escherichia-Shigella" %in% taxa,
          "f__Bifidobacteriaceae"   %in% taxa,
          "f__Bifidobacteriaceae.1" %in% taxa)

# =============================================================================
# 1. ESTRUCTURA DE LA COHORTE  (sujetos, muestras, outcome, trimestres)
# =============================================================================
n_subj <- 42L
# id (sujeto) e index (muestra) son ENTEROS, igual que en los datos reales: el .Rmd
# asume index numerico (chunk figure2_pca_calculation: as.integer(rownames(...))).
subj_id <- seq_len(n_subj)            # 1..42

# ~33% preterm a nivel sujeto (14/42). Ambas clases presentes; suficientes
# preterm para que los splits estratificados 5-fold + inner 70/30 tengan ambas
# clases en cada fold (requisito de roc()/optimize_threshold_cv, no tryCatch-eados).
n_preterm <- 14L
preterm_subj <- rep(0L, n_subj)
preterm_subj[sample(seq_len(n_subj), n_preterm)] <- 1L

# Edad gestacional al parto (semanas), nivel SUJETO:
#   term    (preterm=0): 37.0 - 41.0
#   preterm (preterm=1): 30.0 - 36.9 (una fraccion < 34 -> early_preterm)
sdg_parto_subj <- numeric(n_subj)
sdg_parto_subj[preterm_subj == 0] <- round(runif(sum(preterm_subj == 0), 37.0, 41.0), 1)
sdg_parto_subj[preterm_subj == 1] <- round(runif(sum(preterm_subj == 1), 30.0, 36.8), 1)
early_preterm_subj <- as.integer(sdg_parto_subj < 34)

# Muestras por sujeto (2 o 3), visitas en trimestres distintos para poblar la
# figura de trimestres (T1 <13, T2 13-26, T3 >=27 semanas).
n_visits_subj <- sample(2:3, n_subj, replace = TRUE, prob = c(0.35, 0.65))

samples <- do.call(rbind, lapply(seq_len(n_subj), function(k) {
  nv  <- n_visits_subj[k]
  # GA de visita por trimestre (garantiza cobertura de los 3 trimestres en la cohorte)
  tri_pool <- list(
    T1 = runif(1,  7, 12.5),   # primer trimestre
    T2 = runif(1, 15, 25),     # segundo trimestre
    T3 = runif(1, 28, 34)      # tercer trimestre
  )
  if (nv == 3) {
    gas <- c(tri_pool$T1, tri_pool$T2, tri_pool$T3)
  } else {
    # 2 visitas: elegir 2 trimestres distintos (mezcla para cubrir todos en la cohorte)
    picks <- sample(c("T1", "T2", "T3"), 2)
    gas <- sort(vapply(picks, function(p) tri_pool[[p]], numeric(1)))
  }
  # una visita nunca despues del parto
  gas <- pmin(gas, sdg_parto_subj[k] - 0.5)
  data.frame(
    id           = subj_id[k],
    visita       = seq_len(nv),
    sdg_visita   = round(gas, 1),
    stringsAsFactors = FALSE
  )
}))
rownames(samples) <- NULL
n_samples <- nrow(samples)
samples$index <- seq_len(n_samples)          # index entero unico por muestra (1..n)

# si = posicion del sujeto de cada muestra en los vectores subject-level (1..n_subj).
# Todas las expansiones sujeto->muestra usan samples-de-sujeto via si (sin ambiguedad
# nombre/posicion).
si <- match(samples$id, subj_id)

# Propagar atributos de sujeto a cada muestra
samples$preterm       <- preterm_subj[si]
samples$sdg_parto     <- sdg_parto_subj[si]
samples$early_preterm <- early_preterm_subj[si]

# =============================================================================
# 2. VARIABLES CLINICAS  (respetando type/unit/dominio del diccionario)
# =============================================================================
# --- Nivel SUJETO (constantes en las muestras de un mismo sujeto) ------------
# Edad: relacion en U con PTB (extremos <20 y >=35 mas frecuentes en preterm) -> senal debil.
edad_subj <- numeric(n_subj); names(edad_subj) <- subj_id
for (k in seq_len(n_subj)) {
  if (preterm_subj[k] == 1 && runif(1) < 0.55) {
    edad_subj[k] <- if (runif(1) < 0.5) round(runif(1, 16, 20)) else round(runif(1, 35, 43))
  } else {
    edad_subj[k] <- round(rnorm(1, 27, 5))
  }
}
edad_subj <- pmin(pmax(edad_subj, 15), 45)

# IMC pregestacional: preterm levemente mas alto (senal debil).
imc_preg_subj <- ifelse(preterm_subj == 1, rnorm(n_subj, 27.0, 4.2), rnorm(n_subj, 24.2, 3.4))
imc_preg_subj <- round(pmin(pmax(imc_preg_subj, 16.5), 41), 1)
talla_subj    <- round(rnorm(n_subj, 158, 6), 0)                       # cm
peso_preg_subj <- round(imc_preg_subj * (talla_subj / 100)^2, 1)       # kg (consistente con imc)

imc_categ_subj <- cut(imc_preg_subj, breaks = c(-Inf, 18.5, 25, 30, Inf),
                      labels = c("Bajo peso", "Normal", "Sobrepeso", "Obesidad"),
                      right = FALSE)
imc_categ_subj <- as.character(imc_categ_subj)

nivel_acad_subj <- sample(1:4, n_subj, replace = TRUE, prob = c(0.15, 0.35, 0.30, 0.20))
marital_subj    <- sample(c(1, 2, 5), n_subj, replace = TRUE, prob = c(0.55, 0.10, 0.35))
workout_subj    <- sample(c(1, 2), n_subj, replace = TRUE, prob = c(0.45, 0.55))
sex_baby_subj   <- sample(c("F", "M"), n_subj, replace = TRUE)

# Complicaciones obstetricas (0/1), nivel sujeto; algunas mas frecuentes en preterm.
bern <- function(p_term, p_pt) rbinom(n_subj, 1, ifelse(preterm_subj == 1, p_pt, p_term))
comp1tribleed_subj <- bern(0.12, 0.36)
compvaginf_subj    <- bern(0.15, 0.42)
compsexualinf_subj <- bern(0.06, 0.14)
comppreeclam_subj  <- bern(0.05, 0.22)
rpm_preterm_subj   <- as.integer(preterm_subj == 1 & rbinom(n_subj, 1, 0.30) == 1)
rpm_subj           <- as.integer(preterm_subj == 0 & rbinom(n_subj, 1, 0.10) == 1)
diabetes_gest_subj <- bern(0.07, 0.16)
oligohidramnios_subj <- bern(0.05, 0.15)
rciu_subj          <- bern(0.06, 0.20)
obito_subj         <- rbinom(n_subj, 1, 0.02)
bajo_peso_nac_subj <- as.integer(preterm_subj == 1 & rbinom(n_subj, 1, 0.5) == 1)

# Suplementacion (0/1) y frecuencia, nivel sujeto.
dietsuppl3mon_subj <- rbinom(n_subj, 1, 0.55)
vitaminsup_subj    <- rbinom(n_subj, 1, 0.60)
supintakefreq_subj <- round(ifelse(vitaminsup_subj == 1, runif(n_subj, 0.5, 2.0), 0), 1)

# Peso al nacer y desenlace (nivel sujeto)
birthweight_subj <- round(ifelse(preterm_subj == 1, rnorm(n_subj, 2450, 450),
                                 rnorm(n_subj, 3250, 380)))
birthweight_subj <- pmax(birthweight_subj, 900)
peso_nac_categ_subj <- ifelse(birthweight_subj < 2500, "Bajo",
                       ifelse(birthweight_subj > 4000, "Macrosomia", "Adecuado"))
desenlace_subj <- ifelse(preterm_subj == 1, "Pretermino", "Termino")

# --- Nivel VISITA (varian por muestra) ---------------------------------------
ga <- samples$sdg_visita
# Peso e IMC a la visita: crecen con la EG.
peso_visita <- round(peso_preg_subj[si] + (ga / 40) * rnorm(n_samples, 9, 2), 1)
imc_visita  <- round(peso_visita / (talla_subj[si] / 100)^2, 1)

# Hemoglobina (g/dL) y ajuste por altitud (CDMX ~ -1.1); anemia si Hb ajustada < 11.
hb <- round(rnorm(n_samples, 12.6, 1.2), 1)
hb <- pmin(pmax(hb, 8.5), 16)
hb_alti_adj <- round(hb - runif(n_samples, 0.9, 1.5), 1)          # ruido -> r<0.95 con hb
anemia_visita <- as.integer(hb_alti_adj < 11)

# Dieta / recordatorio 24 h (nivel visita), valores plausibles positivos.
pos <- function(mean, sd) pmax(round(rnorm(n_samples, mean, sd), 1), 0)
calories      <- pos(2050, 420)
calories_fat  <- pos(620, 160)
carbohydrates <- pos(270, 60)
protein       <- pos(78, 20)
fat           <- pos(69, 18)
vitamin_b1    <- pos(1.3, 0.4)
vitamin_b2    <- pos(1.5, 0.4)
vitamin_b6    <- pos(1.6, 0.5)
vitamin_b12   <- pos(3.4, 1.2)
folic_ac_correg <- pos(520, 160)
folate        <- pos(420, 120)
folate_dfe    <- pos(640, 180)
folate_food   <- pos(300, 90)
choline       <- pos(320, 80)
cystine       <- pos(1.1, 0.3)
glycine       <- pos(3.2, 0.9)
methionine    <- pos(1.7, 0.5)
serine        <- pos(3.0, 0.8)

# Medidas fetales (nivel visita), crecen con la EG.
pfetal <- round(pmax(exp(0.16 * ga) * runif(n_samples, 0.9, 1.1), 5))    # g (crecimiento exp)
fcf    <- round(rnorm(n_samples, 143, 8))                                # bpm
ccef   <- round(ga * 8.5 + rnorm(n_samples, 0, 6))                       # mm
dbip   <- round(ga * 2.4 + rnorm(n_samples, 0, 2))                       # mm
cabd   <- round(ga * 8.9 + rnorm(n_samples, 0, 7))                       # mm
longfe <- round(ga * 1.9 + rnorm(n_samples, 0, 2))                       # mm

# Features rezagadas / rolling / tiempo (NO usadas por ningun chunk; se emiten
# con valores plausibles y completos para contract-completeness).
ord <- order(samples$id, samples$visita)
lag_peso <- rep(NA_real_, n_samples); lag_hb <- rep(NA_real_, n_samples)
lag_imc  <- rep(NA_real_, n_samples); roll_peso <- rep(NA_real_, n_samples)
roll_imc <- rep(NA_real_, n_samples); tsf <- rep(NA_real_, n_samples)
for (id in subj_id) {
  idx <- which(samples$id == id)
  idx <- idx[order(samples$visita[idx])]
  p <- peso_visita[idx]; h <- hb[idx]; iv <- imc_visita[idx]; g <- ga[idx]
  lag_peso[idx] <- c(p[1], head(p, -1))          # primera visita = valor actual (sin NA)
  lag_hb[idx]   <- c(h[1], head(h, -1))
  lag_imc[idx]  <- c(iv[1], head(iv, -1))
  roll_peso[idx] <- cumsum(p) / seq_along(p)      # promedio movil acumulado
  roll_imc[idx]  <- cumsum(iv) / seq_along(iv)
  tsf[idx]       <- g - g[1]                       # semanas desde la primera visita
}

# =============================================================================
# 3. MICROBIOMA: conteos (multinomial) -> abundancias relativas
# =============================================================================
# 3.1 Prevalencia objetivo por taxon (fraccion de muestras donde esta presente).
contaminants <- c("Methylobacterium", "Methylorubrum", "Ralstonia", "Mesorhizobium",
                  "Microbacterium", "Bradyrhizobium", "Sphingomonas", "Pseudomonas",
                  "Acinetobacter", "o__Chloroplast")
signal_taxa  <- c("Gardnerella", "Prevotella", "Sneathia")   # elevados en dysbiosis/preterm

prev <- setNames(numeric(n_taxa), taxa)
prev["Lactobacillus"]  <- 1.0                                # dominante, siempre presente
prev[signal_taxa]      <- c(0.70, 0.62, 0.52)
prev[contaminants]     <- runif(length(contaminants), 0.08, 0.22)   # pasan prevalencia -> se quitan por decontaminante

assigned <- c("Lactobacillus", signal_taxa, contaminants)
rest <- setdiff(taxa, assigned)

# core intermitente (prev alta, pasan filtro): ~18 taxa
n_core <- 18L
core_more <- sample(rest, n_core); rest <- setdiff(rest, core_more)
prev[core_more] <- runif(n_core, 0.30, 0.70)

# intermitentes (prev media, pasan filtro >=0.05): ~25 taxa
n_interm <- 25L
interm <- sample(rest, n_interm); rest <- setdiff(rest, interm)
prev[interm] <- runif(n_interm, 0.15, 0.35)

# raros (presentes en 3-5 muestras -> NO pasan el filtro de prevalencia del 5%
# [necesita >= ceiling(0.05*n) muestras], pero con >=2 valores positivos por taxon,
# que es lo que exige cmultRepl(GBM) en la PCA de la figura 2 sobre los 97 taxa).
prev[rest] <- sample(3:5, length(rest), replace = TRUE) / n_samples

# 3.2 meanlog por tier (peso relativo dentro de la fraccion no-Lactobacillus).
tier_ml <- setNames(rep(-2.0, n_taxa), taxa)     # raros: peso bajo
tier_ml[c(signal_taxa, contaminants, core_more)] <- 0.0
tier_ml[interm] <- -1.0

# 3.3 Conjuntos de presencia por taxon (tamano = round(prev * n_samples), >=1).
present_sets <- lapply(taxa, function(t) {
  k <- max(1L, round(prev[t] * n_samples))
  if (t == "Lactobacillus") seq_len(n_samples) else sort(sample(seq_len(n_samples), min(k, n_samples)))
})
names(present_sets) <- taxa
# matriz logica muestras x taxa: TRUE si el taxon esta "presente" en la muestra
present_mat <- matrix(FALSE, n_samples, n_taxa, dimnames = list(NULL, taxa))
for (t in taxa) present_mat[present_sets[[t]], t] <- TRUE

# 3.4 Dysbiosis (mas frecuente en preterm) -> senal microbioma<->preterm (debil).
dys_prob   <- ifelse(samples$preterm == 1, 0.55, 0.25)
dysbiotic  <- rbinom(n_samples, 1, dys_prob)
dys_sev    <- ifelse(dysbiotic == 1, runif(n_samples, 0.5, 1.0), runif(n_samples, 0.0, 0.3))
lac_prop   <- ifelse(dysbiotic == 1, runif(n_samples, 0.10, 0.50), runif(n_samples, 0.70, 0.95))

# 3.5 Profundidad de secuenciacion por muestra (~1e4-5e4 reads).
depth <- round(runif(n_samples, 15000, 45000))

# 3.6 Generar conteos por muestra (multinomial sobre el vector de probabilidad).
counts <- matrix(0L, n_samples, n_taxa, dimnames = list(samples$index, taxa))
for (i in seq_len(n_samples)) {
  pres <- taxa[present_mat[i, ]]
  non_lac <- setdiff(pres, "Lactobacillus")
  raw <- rlnorm(length(non_lac), meanlog = tier_ml[non_lac], sdlog = 0.8)
  names(raw) <- non_lac
  # boost de senal en muestras dysbioticas
  boost <- non_lac %in% signal_taxa
  raw[boost] <- raw[boost] * (1 + 6 * dys_sev[i])
  p <- setNames(numeric(n_taxa), taxa)
  p["Lactobacillus"] <- lac_prop[i]
  if (length(non_lac) > 0) p[non_lac] <- raw / sum(raw) * (1 - lac_prop[i])
  p <- p / sum(p)
  counts[i, ] <- as.integer(stats::rmultinom(1, size = depth[i], prob = p))
}

# 3.7 Seguridad: cada taxon con >= min_pos muestras positivas. cmultRepl(GBM) —que
# corre sobre los 97 taxa en la PCA de la figura 2— falla en columnas con < 2 valores
# positivos ("not enough information to compute t hyper-parameter"). min_pos=3 da margen
# y queda por debajo del umbral de prevalencia (ceiling(0.05*n) muestras), asi que los
# taxa raros se siguen filtrando (no inflan genera_clean).
min_pos <- 3L
for (j in seq_len(n_taxa)) {
  pos_idx <- which(counts[, j] > 0)
  if (length(pos_idx) < min_pos) {
    need <- min_pos - length(pos_idx)
    cand <- setdiff(present_sets[[taxa[j]]], pos_idx)         # preferir su conjunto de presencia
    if (length(cand) < need) cand <- setdiff(seq_len(n_samples), pos_idx)
    add_idx <- cand[seq_len(need)]
    counts[add_idx, j] <- counts[add_idx, j] + (stats::rpois(need, 6) + 1L)
  }
}
stopifnot(all(colSums(counts > 0) >= min_pos), all(rowSums(counts) > 0))

# 3.8 Abundancias relativas = counts / rowSums(counts)  (cada fila suma 1).
rel <- counts / rowSums(counts)

# =============================================================================
# 4. ENSAMBLAR LAS TABLAS NO-TAXA (mismas columnas para rel y abs)
# =============================================================================
# data.frame con check.names=FALSE para preservar nombres exactos.
meta_df <- data.frame(
  id                      = samples$id,
  index                   = samples$index,
  preterm                 = samples$preterm,
  sdg_parto               = samples$sdg_parto,
  visita                  = samples$visita,
  dg_visita               = round(samples$sdg_visita * 7),
  sdg_visita              = samples$sdg_visita,
  dg_parto                = round(samples$sdg_parto * 7),
  desenlace_parto         = desenlace_subj[si],
  edad_cronologicamujer   = edad_subj[si],
  peso_pregestacional_kg  = peso_preg_subj[si],
  talla_mujer_cm          = talla_subj[si],
  imc_pregestacional      = imc_preg_subj[si],
  imc_pregest_categ       = imc_categ_subj[si],
  peso_kg                 = peso_visita,
  calories                = calories,
  calories_fat            = calories_fat,
  carbohydrates           = carbohydrates,
  protein                 = protein,
  fat                     = fat,
  vitamin_b1              = vitamin_b1,
  vitamin_b2              = vitamin_b2,
  vitamin_b6              = vitamin_b6,
  vitamin_b12             = vitamin_b12,
  folic_ac_correg         = folic_ac_correg,
  folate                  = folate,
  folate_dfe              = folate_dfe,
  folate_food             = folate_food,
  choline                 = choline,
  cystine                 = cystine,
  glycine                 = glycine,
  methionine              = methionine,
  serine                  = serine,
  dietsuppl3mon           = dietsuppl3mon_subj[si],
  vitaminsup              = vitaminsup_subj[si],
  supintakefreq           = supintakefreq_subj[si],
  comp1tribleed           = comp1tribleed_subj[si],
  compvaginf              = compvaginf_subj[si],
  compsexualinf           = compsexualinf_subj[si],
  comppreeclam            = comppreeclam_subj[si],
  rpm_preterm             = rpm_preterm_subj[si],
  rpm                     = rpm_subj[si],
  diabetes_gest           = diabetes_gest_subj[si],
  obito                   = obito_subj[si],
  oligohidramnios         = oligohidramnios_subj[si],
  rciu                    = rciu_subj[si],
  bajo_peso_nac           = bajo_peso_nac_subj[si],
  nivel_academico         = nivel_acad_subj[si],
  maritalstat             = marital_subj[si],
  workouthome             = workout_subj[si],
  sex_baby                = sex_baby_subj[si],
  birthweightgr           = birthweight_subj[si],
  peso_nacimiento         = peso_nac_categ_subj[si],
  pfetal                  = pfetal,
  fcf                     = fcf,
  ccef                    = ccef,
  dbip                    = dbip,
  cabd                    = cabd,
  longfe                  = longfe,
  hemoglobin_g_dl         = hb,
  hemoglobin_alti_adj     = hb_alti_adj,
  anemia_visita           = anemia_visita,
  early_preterm           = samples$early_preterm,
  imc_visita              = imc_visita,
  lag_peso_kg             = round(lag_peso, 1),
  lag_hemoglobin          = round(lag_hb, 1),
  lag_imc_visita          = round(lag_imc, 1),
  rolling_avg_peso        = round(roll_peso, 1),
  rolling_avg_imc         = round(roll_imc, 1),
  time_since_first        = round(tsf, 1),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Sanity: cubrir exactamente las 70 no-taxa del diccionario (orden libre).
stopifnot(setequal(names(meta_df), nontaxa))

# =============================================================================
# 5. ESCRIBIR LOS 3 CSV  (nombres de taxa EXACTOS; hyphen y .1 preservados)
# =============================================================================
assemble <- function(taxa_values_matrix) {
  taxa_df <- as.data.frame(taxa_values_matrix, check.names = FALSE, stringsAsFactors = FALSE)
  full <- cbind(meta_df, taxa_df)
  names(full) <- c(names(meta_df), taxa)   # forzar nombres exactos post-cbind
  full
}

rel_out <- assemble(rel)
abs_out <- assemble(counts)
meta_long_out <- data.frame(
  id         = samples$id,
  index      = samples$index,
  preterm    = samples$preterm,
  sdg_parto  = samples$sdg_parto,
  sdg_visita = samples$sdg_visita,
  check.names = FALSE, stringsAsFactors = FALSE
)

readr::write_csv(rel_out,        file.path(out_dir, "example_genus_rel.csv"))
readr::write_csv(abs_out,        file.path(out_dir, "example_genus_abs.csv"))
readr::write_csv(meta_long_out,  file.path(out_dir, "example_metadata_long.csv"))

# =============================================================================
# 6. DIAGNOSTICOS
# =============================================================================
pct_zero <- mean(counts == 0) * 100
prev_pass <- sum(colSums(counts > 0) >= ceiling(0.05 * n_samples))
tri <- cut(samples$sdg_visita, c(-Inf, 13, 27, Inf),
           labels = c("T1", "T2", "T3"), right = FALSE)

cat("\n=== generate_example_data.R — RESUMEN ===\n")
cat(sprintf("Sujetos: %d  (preterm=%d, %.1f%%; early_preterm=%d)\n",
            n_subj, n_preterm, 100 * n_preterm / n_subj, sum(early_preterm_subj)))
cat(sprintf("Muestras: %d  (muestras preterm=%d)\n", n_samples, sum(samples$preterm)))
cat(sprintf("Taxa: %d  | %% ceros en conteos: %.1f%%\n", n_taxa, pct_zero))
cat(sprintf("Taxa que pasan prevalencia (>=5%% muestras): %d  (esperado genera_clean ~%d tras quitar %d contaminantes)\n",
            prev_pass, prev_pass - sum(colSums(counts[, contaminants, drop = FALSE] > 0) >= ceiling(0.05 * n_samples)),
            length(contaminants)))
cat("Cobertura de trimestres (muestras):\n"); print(table(tri))
cat(sprintf("rowSums(rel) en [%.6f, %.6f] (debe ~1)\n", min(rowSums(rel)), max(rowSums(rel))))
cat(sprintf("Columnas rel: %d (70 no-taxa + 97 taxa = 167)\n", ncol(rel_out)))
cat(sprintf("NA en rel_out: %d | NA en abs_out: %d\n", sum(is.na(rel_out)), sum(is.na(abs_out))))
cat("Escritos en: ", normalizePath(out_dir), "\n", sep = "")
cat("  example_genus_rel.csv, example_genus_abs.csv, example_metadata_long.csv\n\n")
