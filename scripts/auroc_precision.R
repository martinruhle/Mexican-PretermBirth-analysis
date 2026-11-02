#!/usr/bin/env Rscript
# =============================================================================
# auroc_precision.R  —  Tablas de docs/AUROC_PRECISION.md
# -----------------------------------------------------------------------------
# QUE ES: cuanta precision tiene el AUROC con el tamano de esta cohorte, y cual es la
# menor diferencia de AUROC entre dos modelos evaluados en las mismas participantes
# (con y sin microbioma) que se detectaria con alfa 0.05 bilateral y potencia 80 %.
#
# PARTES (A-E no necesitan los datos de la cohorte; F si):
#   A. formulas de Hanley-McNeil (HM) y del aporte detectable pareado
#   B. reproduccion del calculo previo: delta = 2.8 * EE_HM(A), el mismo EE para los
#      dos modelos
#   C. contraste por simulacion: scores binormales pareados, test de DeLong pareado
#   D. el procedimiento de validacion cruzada: SD de la nula de permutacion versionada
#      (results/permutation_test/permutation_null_max18.csv: las 18 combinaciones bajo
#      la misma permutacion) contra la SD que daria solo el muestreo de los folds de test,
#      y SD entre folds observada (results/clinical_only_arm/metrics_by_fold.csv) contra
#      la esperada por ese muestreo
#   E. tabla por tamano de cohorte
#   F. conteos por ventana de semanas y correlacion r entre los AUROC de los pares con y
#      sin microbioma (DeLong sobre las predicciones guardadas de la corrida del brazo
#      solo clinico). Sin data/raw/ o sin esa corrida, F se saltea y la tabla usa r fijas.
#
# sdg_visita / sdg_parto estan en notacion semanas.dias (24.4 = 24 semanas y 4 dias):
# comparar con un entero ("< 32" = antes de 32+0, "< 33" = hasta 32+6) da lo mismo que
# en semanas decimales.
#
# CORRER (desde la raiz del repo):
#   PTB_RESULTS_DIR=results/auroc_precision Rscript scripts/auroc_precision.R
# Variables: PTB_NSIM (replicas por punto de la simulacion, 2000), PTB_PAIRS_RUN_TAG
# (corrida con las 18 combinaciones, clinonly_2026-10-16).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(rsample)
})
source(file.path("R", "io.R"))
source(file.path("R", "nested_cv.R"))   # roc_ptb(): la misma metrica que el motor

results_dir <- Sys.getenv("PTB_RESULTS_DIR", "")
nsim        <- as.integer(Sys.getenv("PTB_NSIM", "2000"))
pairs_tag   <- Sys.getenv("PTB_PAIRS_RUN_TAG", "clinonly_2026-10-16")
perm_csv    <- file.path("results", "permutation_test", "permutation_null_max18.csv")
fold_csv    <- file.path("results", "clinical_only_arm", "metrics_by_fold.csv")
stopifnot(file.exists(perm_csv), file.exists(fold_csv), nsim >= 200)

alpha <- 0.05; power <- 0.80
za <- qnorm(1 - alpha / 2); zb <- qnorm(power)
A_grid <- c(0.60, 0.70, 0.80)
f3 <- function(x) ifelse(is.na(x), "n/r", sprintf("%.3f", x))
out <- list()

# --- A. formulas ----------------------------------------------------------------
# Varianza del AUROC de Hanley & McNeil (1982), con Q1 y Q2 del modelo exponencial
hm_var <- function(A, n1, n0) {
  q1 <- A / (2 - A); q2 <- 2 * A^2 / (1 + A)
  (A * (1 - A) + (n1 - 1) * (q1 - A^2) + (n0 - 1) * (q2 - A^2)) / (n1 * n0)
}
# Media de K AUROC por fold con un modelo fijo: Var = sum(V_k) / K^2
hm_var_folds <- function(A, n1k, n0k) sum(hm_var(A, n1k, n0k)) / length(n1k)^2

# Metodo previo: delta = 2.8 * EE(A). Con r = 0.5 y el mismo EE para los dos modelos,
# sqrt(EE^2 + EE^2 - 2 * 0.5 * EE^2) = EE.
delta_prior <- function(A, n1, n0) 2.8 * sqrt(hm_var(A, n1, n0))

# Aporte detectable pareado (varianza de la diferencia de dos AUROC correlacionados de
# Hanley & McNeil 1983), con V(a) = var_fun(a):
#   V0 = 2 V(A) (1 - r)                              bajo H0 (los dos modelos en A)
#   V1 = V(A) + V(A + d) - 2 r sqrt(V(A) V(A + d))  bajo H1
# critical = "null": z_a sqrt(V0) + z_b sqrt(V1) = d (Obuchowski 1998: el valor critico
#   sale de la varianza bajo H0; es el calculo que se fijo en el borrador)
# critical = "alternative": (z_a + z_b) sqrt(V1) = d (test de Wald con la varianza estimada
#   de los datos, como el test de DeLong de pROC::roc.test; agregado despues de ver la
#   simulacion)
# NA = ningun delta con A + delta <= 1 alcanza la potencia.
delta_paired <- function(A, var_fun, r, critical = c("null", "alternative")) {
  critical <- match.arg(critical)
  V <- var_fun
  v1 <- function(d) max(V(A) + V(A + d) - 2 * r * sqrt(V(A) * V(A + d)), 0)
  g <- if (critical == "null") function(d) za * sqrt(2 * V(A) * (1 - r)) + zb * sqrt(v1(d)) - d else
    function(d) (za + zb) * sqrt(v1(d)) - d
  if (g(1 - A) > 0) return(NA_real_)
  stats::uniroot(g, c(0, 1 - A), tol = 1e-12)$root
}
pooled_var <- function(n1, n0) function(a) hm_var(a, n1, n0)

# --- B. reproduccion del calculo previo -----------------------------------------
cat("== B. Prior computation: delta = 2.8 * SE_HM(A), same SE for both models ==\n")
prior_rows <- data.frame(
  row = c("This cohort", "This cohort + DREAM validation Hispanic subgroup",
          "~130 (1/3 preterm)", "~200 (1/3 preterm)", "~300 (1/3 preterm)", "~500 (1/3 preterm)"),
  n1 = c(14, 29, 43, 67, 100, 167), n0 = c(29, 52, 87, 133, 200, 333))
stated <- data.frame(row = prior_rows$row[c(1, 1, 2:6)], A = c(0.70, 0.80, rep(0.70, 5)),
                     stated = c(0.251, 0.220, 0.176, 0.142, 0.114, 0.093, 0.072))
prior <- merge(expand.grid(row = prior_rows$row, A = c(0.70, 0.80), stringsAsFactors = FALSE),
               prior_rows, by = "row")
prior$SE <- sqrt(hm_var(prior$A, prior$n1, prior$n0))
prior$CI_halfwidth <- 1.96 * prior$SE
prior$delta_prior  <- delta_prior(prior$A, prior$n1, prior$n0)
prior$A_plus_delta <- prior$A + prior$delta_prior
prior <- merge(prior, stated, by = c("row", "A"), all.x = TRUE)
prior$reproduced <- ifelse(is.na(prior$stated), NA, round(prior$delta_prior, 3) == prior$stated)
prior <- prior[order(prior$A, prior$n1), ]
print(prior[, c("row", "A", "n1", "n0", "SE", "CI_halfwidth", "delta_prior", "A_plus_delta",
                "stated", "reproduced")], row.names = FALSE, digits = 4)
out$prior_reproduction <- prior

# El mismo metodo con la proporcion de la cohorte (14/43), n1 sin redondear y redondeado
prop_rows <- bind_rows(lapply(c(81, 130, 200, 300, 500), function(N) data.frame(
  N = N, n1_exact = N * 14 / 43, n1_rounded = round(N * 14 / 43),
  delta_prior_exact = delta_prior(0.70, N * 14 / 43, N * 29 / 43),
  delta_prior_rounded = delta_prior(0.70, round(N * 14 / 43), N - round(N * 14 / 43)))))
cat("\nPrior method at A = 0.70 with the cohort proportion 14/43:\n")
print(prop_rows, row.names = FALSE, digits = 4)
out$prior_cohort_proportion <- prop_rows

# --- C. contraste por simulacion ------------------------------------------------
# DeLong pareado (higher score = preterm, como roc_ptb()). Devuelve los dos AUROC y su
# matriz de covarianza.
delong2 <- function(s1, s2, y) {
  psi <- function(s) { a <- s[y]; b <- s[!y]; outer(a, b, ">") + 0.5 * outer(a, b, "==") }
  P1 <- psi(s1); P2 <- psi(s2)
  V10 <- cbind(rowMeans(P1), rowMeans(P2)); V01 <- cbind(colMeans(P1), colMeans(P2))
  list(auc = c(mean(P1), mean(P2)), S = stats::cov(V10) / sum(y) + stats::cov(V01) / sum(!y))
}

# Verificacion contra pROC en datos al azar
set.seed(2310)
chk <- replicate(20, {
  y <- sample(c(rep(TRUE, 14), rep(FALSE, 29))); s1 <- rnorm(43) + y; s2 <- 0.6 * s1 + rnorm(43)
  d <- delong2(s1, s2, y)
  r1 <- roc_ptb(as.integer(y), s1); r2 <- roc_ptb(as.integer(y), s2)
  z_mine <- (d$auc[1] - d$auc[2]) / sqrt(d$S[1, 1] + d$S[2, 2] - 2 * d$S[1, 2])
  z_proc <- unname(pROC::roc.test(r1, r2, method = "delong", paired = TRUE)$statistic)
  c(abs(d$S[1, 1] - pROC::var(r1, method = "delong")),
    abs(d$S[1, 2] - pROC::cov(r1, r2, method = "delong")), abs(z_mine - z_proc))
})
cat(sprintf("\n== C. Simulation ==\nDeLong implementation vs pROC (20 random data sets): max |diff| %.1e\n",
            max(chk)))
stopifnot(max(chk) < 1e-10)

# Composicion de los folds externos (la misma en las 18 combinaciones)
fb <- utils::read.csv(fold_csv, stringsAsFactors = FALSE)
comp <- unique(fb[, c("fold", "n_test", "n_ptb")])
stopifnot(nrow(comp) == 5, !anyDuplicated(comp$fold))
comp <- comp[order(comp$fold), ]
n1k <- comp$n_ptb; n0k <- comp$n_test - comp$n_ptb
n1 <- sum(n1k); n0 <- sum(n0k)
cat(sprintf("Outer folds (preterm / term): %s; total %d / %d\n",
            paste(sprintf("%d/%d", n1k, n0k), collapse = ", "), n1, n0))

# Una replica: AUROC agrupado de los dos modelos, z de DeLong y media de los AUROC por fold
sim_once <- function(mu1, mu2, rho, y, fold_y1, fold_y0) {
  e1 <- rnorm(length(y)); e2 <- rho * e1 + sqrt(1 - rho^2) * rnorm(length(y))
  s1 <- e1 + mu1 * y; s2 <- e2 + mu2 * y
  d <- delong2(s1, s2, y)
  vd <- d$S[1, 1] + d$S[2, 2] - 2 * d$S[1, 2]
  fm <- function(s) { P <- outer(s[y], s[!y], ">")
    mean(vapply(seq_along(n1k), function(k) mean(P[fold_y1 == k, fold_y0 == k, drop = FALSE]), 0)) }
  c(auc1 = d$auc[1], auc2 = d$auc[2], z = (d$auc[2] - d$auc[1]) / sqrt(vd),
    fold1 = fm(s1), fold2 = fm(s2))
}
y_sim <- c(rep(TRUE, n1), rep(FALSE, n0))
fold_y1 <- rep(seq_along(n1k), n1k); fold_y0 <- rep(seq_along(n0k), n0k)
run_sim <- function(A1, A2, rho, n = nsim)
  t(replicate(n, sim_once(sqrt(2) * qnorm(A1), sqrt(2) * qnorm(A2), rho, y_sim, fold_y1, fold_y0)))

# rho (correlacion de los scores dentro de cada clase) que produce la r deseada entre los
# AUROC estimados, bajo H0 (los dos en A)
rho_for_r <- function(A, r_target) {
  rhos <- seq(0.3, 0.8, by = 0.05)
  rs <- vapply(rhos, function(rho) { s <- run_sim(A, A, rho, n = 2 * nsim); stats::cor(s[, "auc1"], s[, "auc2"]) }, 0)
  stats::approx(rs, rhos, xout = r_target, ties = mean)$y
}

cross0 <- function(gap, grid) if (any(gap <= 0) && any(gap > 0))
  stats::approx(gap, grid, xout = 0, ties = mean)$y else NA_real_
set.seed(2311)
sim_rows <- list(); sim_fold <- list()
r_targets <- c(0.5)
for (A in A_grid) for (rt in r_targets) {
  rho <- rho_for_r(A, rt)
  d_an <- delta_paired(A, pooled_var(n1, n0), rt)
  s0 <- run_sim(A, A, rho)
  r_emp <- stats::cor(s0[, "auc1"], s0[, "auc2"])
  # grilla desde 0.10 por debajo del delta analitico hasta casi 1 - A
  hi <- 1 - A - 0.005
  lo <- max(0.01, (if (is.na(d_an)) hi - 0.15 else d_an - 0.10))
  grid <- seq(lo, hi, by = 0.01)
  gs <- lapply(grid, function(d) run_sim(A, A + d, rho))
  pw <- vapply(gs, function(x) mean(abs(x[, "z"]) > za), 0)
  d_sim <- if (any(pw >= power)) stats::approx(pw, grid, xout = power, ties = mean)$y else NA_real_
  # "oraculo": la misma ecuacion del metodo elegido con las SD verdaderas (simuladas) de la
  # diferencia bajo H0 y bajo H1, en vez de las de Hanley-McNeil
  sd0 <- stats::sd(s0[, "auc2"] - s0[, "auc1"])
  sd1 <- vapply(gs, function(x) stats::sd(x[, "auc2"] - x[, "auc1"]), 0)
  d_oracle <- cross0(za * sd0 + zb * sd1 - grid, grid)
  d_oracle_w <- cross0((za + zb) * sd1 - grid, grid)
  # lo mismo para la media de los AUROC por fold (la metrica del motor); se multiplica
  # por f en la parte D
  sim_fold[[length(sim_fold) + 1]] <- list(
    grid = grid, sd0 = stats::sd(s0[, "fold2"] - s0[, "fold1"]),
    sd1 = vapply(gs, function(x) stats::sd(x[, "fold2"] - x[, "fold1"]), 0))
  sim_rows[[length(sim_rows) + 1]] <- data.frame(
    A = A, r_target = rt, rho_scores = rho, r_empirical = r_emp,
    type1_error_delong = mean(abs(s0[, "z"]) > za),
    SE_pooled_sim = stats::sd(s0[, "auc1"]), SE_pooled_HM = sqrt(hm_var(A, n1, n0)),
    SE_foldmean_sim = stats::sd(s0[, "fold1"]), SE_foldmean_HM = sqrt(hm_var_folds(A, n1k, n0k)),
    SD_diff_pooled_sim = stats::sd(s0[, "auc2"] - s0[, "auc1"]),
    SD_diff_foldmean_sim = stats::sd(s0[, "fold2"] - s0[, "fold1"]),
    SD_diff_HM = sqrt(2 * hm_var(A, n1, n0) * (1 - rt)),
    delta_analytic = d_an, delta_oracle_binormal = d_oracle,
    delta_analytic_wald = delta_paired(A, pooled_var(n1, n0), rt, "alternative"),
    delta_oracle_binormal_wald = d_oracle_w, delta_simulated_delong = d_sim,
    power_max_on_grid = max(pw), delta_prior = delta_prior(A, n1, n0))
}
sim_tab <- bind_rows(sim_rows)
print(as.data.frame(sim_tab), row.names = FALSE, digits = 3)
out$simulation_check <- sim_tab

# --- D. el procedimiento de validacion cruzada ----------------------------------
cat("\n== D. Cross-validation procedure ==\n")
nul <- utils::read.csv(perm_csv, check.names = FALSE, stringsAsFactors = FALSE)
keys <- setdiff(names(nul), c("seed", "n_folds_valid", "max_auroc", "max_combination"))
stopifnot(length(keys) == 18)

# Nula solo por el muestreo de los folds de test: etiquetas permutadas a nivel participante
# sobre los mismos folds, scores continuos sin relacion con la etiqueta. Dentro de un fold
# con n1/n0, Var(AUROC) = (n1 + n0 + 1) / (12 n1 n0); los folds de una sola clase se
# saltean, como en el motor.
set.seed(2312)
fold_of <- rep(seq_along(n1k), n1k + n0k)
n_draw <- 200000
draws <- replicate(n_draw, {
  lab <- sample(c(rep(1L, n1), rep(0L, n0)))
  a <- tabulate(fold_of[lab == 1L], nbins = length(n1k)); b <- (n1k + n0k) - a
  ok <- a > 0 & b > 0
  c(var = sum((a[ok] + b[ok] + 1) / (12 * a[ok] * b[ok])) / sum(ok)^2, nvalid = sum(ok))
})
V_test_perm <- mean(draws["var", ])
p_skip_expected <- mean(draws["nvalid", ] < length(n1k))
V_test_strat <- sum((n1k + n0k + 1) / (12 * n1k * n0k)) / length(n1k)^2   # folds estratificados

sd_perm <- vapply(keys, function(k) stats::sd(nul[[k]], na.rm = TRUE), 0)
mean_perm <- vapply(keys, function(k) mean(nul[[k]], na.rm = TRUE), 0)
n_perm <- vapply(keys, function(k) sum(!is.na(nul[[k]])), 0)
infl <- function(rows) sqrt(mean(vapply(keys, function(k) stats::var(nul[[k]][rows], na.rm = TRUE), 0)) /
                            V_test_perm)
f_all <- infl(seq_len(nrow(nul)))
set.seed(2313)
f_boot <- replicate(2000, infl(sample(nrow(nul), replace = TRUE)))
f_ci <- stats::quantile(f_boot, c(0.025, 0.975))
cat(sprintf("Permutations: %d; with a skipped fold: %d (expected share under random permutation %.3f)\n",
            nrow(nul), sum(nul$n_folds_valid < 5), p_skip_expected))
cat(sprintf("Null SD from test-fold sampling alone: %.4f (stratified folds %.4f)\n",
            sqrt(V_test_perm), sqrt(V_test_strat)))
cat(sprintf("Null SD of the procedure (permutation null): median over the 18 combinations %.4f (range %.4f-%.4f)\n",
            stats::median(sd_perm), min(sd_perm), max(sd_perm)))
cat(sprintf("Inflation factor f = SD_procedure / SD_test_sampling (pooled over the 18): %.3f (bootstrap 95%% CI %.3f-%.3f)\n",
            f_all, f_ci[1], f_ci[2]))
perm_tab <- data.frame(combination = keys, n_perm = n_perm, null_mean = mean_perm, null_sd = sd_perm,
                       null_sd_test_sampling = sqrt(V_test_perm), inflation = sd_perm / sqrt(V_test_perm),
                       row.names = NULL)
out$permutation_null_sd <- perm_tab

# Delta detectable de la metrica del motor con varianzas binormales: SD simulada de la
# diferencia de las medias por fold, multiplicada por f
sim_tab$delta_engine_binormal <- vapply(sim_fold, function(z)
  cross0(za * f_all * z$sd0 + zb * f_all * z$sd1 - z$grid, z$grid), 0)
sim_tab$delta_engine_binormal_wald <- vapply(sim_fold, function(z)
  cross0((za + zb) * f_all * z$sd1 - z$grid, z$grid), 0)
sim_tab$delta_engine_HM <- vapply(seq_len(nrow(sim_tab)), function(i) delta_paired(
  sim_tab$A[i], function(a) f_all^2 * hm_var_folds(a, n1k, n0k), sim_tab$r_target[i]), 0)
sim_tab$grid_max <- vapply(sim_fold, function(z) max(z$grid), 0)
cat("Detectable delta of the engine metric (mean of 5 folds, x f), r = 0.5:\n")
print(as.data.frame(sim_tab[, c("A", "delta_engine_HM", "delta_engine_binormal", "delta_engine_binormal_wald",
                                "grid_max")]), row.names = FALSE, digits = 3)
out$simulation_check <- sim_tab

# Pares con y sin microbioma bajo la misma permutacion: SD de la diferencia y r bajo H0
pair_of <- function(k) sub("\\|(ANCOM_Taxa|Full_Microbiome)$", "|None", k)
micro_keys <- keys[!grepl("\\|None$", keys)]
pair_null <- bind_rows(lapply(micro_keys, function(k) {
  a <- nul[[k]]; b <- nul[[pair_of(k)]]; ok <- !is.na(a) & !is.na(b)
  data.frame(combination = k, clinical_only = pair_of(k), n_perm = sum(ok),
             r_null = stats::cor(a[ok], b[ok]), sd_diff_null = stats::sd(a[ok] - b[ok]))
}))
cat(sprintf("Pairs with vs without microbiome under the same permutation (12): r median %.3f (range %.3f-%.3f); SD of the difference median %.4f (range %.4f-%.4f)\n",
            stats::median(pair_null$r_null), min(pair_null$r_null), max(pair_null$r_null),
            stats::median(pair_null$sd_diff_null), min(pair_null$sd_diff_null), max(pair_null$sd_diff_null)))
out$pairs_permutation_null <- pair_null

# Seleccion: maximo de 18 bajo H0 y Bonferroni sobre los 12 pares
cat(sprintf("Null of one combination: mean %.3f (median over the 18); null of the maximum of 18: mean %.3f, SD %.3f\n",
            stats::median(mean_perm), mean(nul$max_auroc, na.rm = TRUE), stats::sd(nul$max_auroc, na.rm = TRUE)))
bonf <- (qnorm(1 - alpha / 2 / 12) + zb) / (za + zb)
cat(sprintf("Bonferroni over 12 pairs multiplies the detectable delta by about %.2f (z %.3f instead of %.3f)\n",
            bonf, qnorm(1 - alpha / 2 / 12), za))

# SD entre folds observada vs la esperada por el muestreo del fold de test (binormal con
# el AUROC medio de cada combinacion)
fb$key <- paste(fb$Model, fb$Approach, fb$Microbiome, sep = "|")
set.seed(2314)
fold_sim <- lapply(split(fb, fb$key), function(d) {
  d <- d[order(d$fold), ]
  A <- mean(d$AUROC)
  mu <- sqrt(2) * qnorm(min(max(A, 0.001), 0.999))
  s2_sim <- replicate(4000, {
    au <- vapply(seq_len(nrow(d)), function(k) {
      x <- rnorm(d$n_ptb[k], mu); w <- rnorm(d$n_test[k] - d$n_ptb[k])
      mean(outer(x, w, ">")) }, 0)
    stats::var(au) })
  list(row = data.frame(combination = d$key[1], AUROC_mean = A, sd_observed = stats::sd(d$AUROC),
                        var_observed = stats::var(d$AUROC), var_expected = mean(s2_sim)),
       s2 = s2_sim)
})
fold_sd <- bind_rows(lapply(fold_sim, `[[`, "row"))
ratio_obs <- sum(fold_sd$var_observed) / sum(fold_sd$var_expected)
ratio_null <- Reduce(`+`, lapply(fold_sim, `[[`, "s2")) / sum(fold_sd$var_expected)
cat(sprintf("Between-fold SD: observed median %.3f, expected from test-fold sampling median %.3f\n",
            stats::median(fold_sd$sd_observed), stats::median(sqrt(fold_sd$var_expected))))
cat(sprintf("Sum over the 18 of observed / expected between-fold variance: %.3f (95%% range if the folds only differed by test sampling, combinations taken as independent: %.3f-%.3f)\n",
            ratio_obs, stats::quantile(ratio_null, 0.025), stats::quantile(ratio_null, 0.975)))
out$between_fold_sd <- as.data.frame(fold_sd)

# --- F. datos de la cohorte: conteos por ventana y r de los pares ----------------
cat("\n== F. Cohort data ==\n")
cfg <- read_config(Sys.getenv("PTB_PROFILE", "default"))
outdir <- cfg$output$dir %||% "analysis/_output"
samples <- tryCatch(as.data.frame(load_dataset(cfg)$matrix), error = function(e) NULL)
r_delong <- NULL; windows <- NULL
if (is.null(samples)) {
  cat("Cohort data not found: window counts and the DeLong correlation SKIPPED\n")
} else {
  cc <- cfg$columns
  s <- data.frame(id = as.character(samples[[cc$subject_id]]),
                  preterm = as.character(samples[[cc$outcome]]) == cc$outcome_positive,
                  ga_delivery = samples[[cc$gestational_age]], ga_visit = samples[[cc$visit_ga]],
                  stringsAsFactors = FALSE)
  stopifnot(!anyNA(s),
            all(tapply(s$ga_delivery, s$id, function(v) length(unique(v))) == 1),
            all(tapply(s$preterm, s$id, function(v) length(unique(v))) == 1),
            all(s$ga_visit <= s$ga_delivery))
  p <- s %>% group_by(id) %>%
    summarise(preterm = first(preterm), ga_delivery = first(ga_delivery),
              first_visit = min(ga_visit), n_samples = n(), .groups = "drop")
  cat(sprintf("Participants %d (preterm %d, term %d); samples %d\n",
              nrow(p), sum(p$preterm), sum(!p$preterm), nrow(s)))
  cat(sprintf("preterm == (delivery before 37+0) in all participants: %s\n",
              all(p$preterm == (p$ga_delivery < 37))))

  # participantes con >= 1 especimen antes de `before` (semanas.dias < before; Inf = todas)
  has_sample <- function(before) p$id %in% unique(s$id[s$ga_visit < before])
  # base = participantes elegibles (todas, o las que llegaron a X); entran las que tienen
  # >= 1 especimen antes de `before`; case = desenlace que cuenta como caso
  win <- function(window, specimens, outcome, base, before, case) {
    keep <- base & has_sample(before)
    data.frame(window = window, specimens = specimens, outcome = outcome,
               eligible = sum(base), participants = sum(keep), cases = sum(keep & case),
               non_cases = sum(keep & !case), preterm = sum(keep & p$preterm), term = sum(keep & !p$preterm),
               eligible_without_sample = sum(base & !keep),
               eligible_without_sample_preterm = sum(base & !keep & p$preterm),
               samples_in_window = sum(s$id %in% p$id[keep] & s$ga_visit < before))
  }
  all_p <- rep(TRUE, nrow(p))
  windows <- bind_rows(
    win("All participants", "all", "preterm (< 37+0)", all_p, Inf, p$preterm),
    win("Deliveries before 32+0", "all", "delivery < 32+0", all_p, Inf, p$ga_delivery < 32),
    win("Deliveries before 34+0", "all", "delivery < 34+0", all_p, Inf, p$ga_delivery < 34),
    win("DREAM SC1 window", "before 32+0", "preterm (< 37+0)", all_p, 32, p$preterm),
    win("DREAM SC1 window", "through 32+6", "preterm (< 37+0)", all_p, 33, p$preterm),
    win("DREAM SC2 window", "before 28+0", "delivery < 32+0", all_p, 28, p$ga_delivery < 32),
    win("DREAM SC2 window", "through 28+6", "delivery < 32+0", all_p, 29, p$ga_delivery < 32),
    bind_rows(lapply(c(24, 28, 32), function(X) bind_rows(
      win(sprintf("Reached week %d", X), "all", "preterm (< 37+0)", p$ga_delivery >= X, Inf, p$preterm),
      win(sprintf("Landmark week %d", X), sprintf("before %d+0", X), "preterm (< 37+0)",
          p$ga_delivery >= X, X, p$preterm),
      win(sprintf("Landmark week %d", X), sprintf("through %d+6", X), "preterm (< 37+0)",
          p$ga_delivery >= X, X + 1, p$preterm)))))
  print(windows, row.names = FALSE)

  # folds reconstruidos como en el chunk nested_cv_setup y comparados con la corrida
  set.seed(cfg$cv$seed)
  subject_labels <- s %>% mutate(preterm = ifelse(preterm, "1", "0")) %>% group_by(id) %>%
    summarize(preterm = ifelse(any(preterm == "1"), "1", "0")) %>% ungroup()
  cv_folds <- vfold_cv(subject_labels, v = cfg$cv$outer_folds, strata = preterm)
  folds <- bind_rows(lapply(seq_len(nrow(cv_folds)), function(i)
    data.frame(id = assessment(cv_folds$splits[[i]])$id, fold = i, stringsAsFactors = FALSE)))
  fold_comp <- folds %>% left_join(p, by = "id") %>% group_by(fold) %>%
    summarise(n1 = sum(preterm), n0 = sum(!preterm), .groups = "drop")
  stopifnot(identical(as.integer(fold_comp$n1), as.integer(n1k)),
            identical(as.integer(fold_comp$n0), as.integer(n0k)))
  cat("Rebuilt outer folds have the composition of results/clinical_only_arm/metrics_by_fold.csv: TRUE\n")

  run_file <- file.path(outdir, paste0(pairs_tag, "_cv_folds.rds"))
  if (!file.exists(run_file)) {
    cat("Run", pairs_tag, "not found in", outdir, ": DeLong correlation SKIPPED\n")
  } else {
    run <- readRDS(run_file)
    names(run) <- vapply(run, function(r) paste(r$model_name, r$approach, r$microbiome, sep = "|"), "")
    stopifnot(setequal(names(run), keys))
    same <- vapply(run, function(r) identical(sort(paste(r$test_predictions$fold, r$test_predictions$id)),
                                              sort(paste(folds$fold, folds$id))), logical(1))
    if (!all(same)) stop("Rebuilt folds differ from run ", pairs_tag)
    cat(sprintf("Rebuilt outer folds identical to the %d combinations of run %s: TRUE\n", length(run), pairs_tag))
    # el AUROC por fold de las predicciones guardadas reproduce el del motor
    eng_ok <- vapply(run, function(r) {
      tp <- r$test_predictions
      mine <- vapply(sort(unique(tp$fold)), function(k) {
        t <- tp[tp$fold == k, ]; as.numeric(pROC::auc(roc_ptb(t$true_class, t$pred_prob))) }, 0)
      fr <- r$fold_results[order(r$fold_results$fold), ]
      identical(mine, as.numeric(fr$AUROC))
    }, logical(1))
    stopifnot(all(eng_ok))
    cat("Per-fold AUROC of the saved predictions reproduces the engine in 18/18 combinations: TRUE\n")

    r_delong <- bind_rows(lapply(micro_keys, function(k) {
      a <- run[[k]]$test_predictions; b <- run[[pair_of(k)]]$test_predictions
      b <- b[match(paste(a$fold, a$id), paste(b$fold, b$id)), ]
      stopifnot(identical(a$id, b$id), identical(as.character(a$true_class), as.character(b$true_class)))
      y <- as.character(a$true_class) == "1"
      d <- delong2(b$pred_prob, a$pred_prob, y)               # 1 = solo clinico, 2 = con microbioma
      # por fold, combinado como la media de los folds: r = sum C_k / sqrt(sum V1k sum V2k)
      Sk <- lapply(sort(unique(a$fold)), function(k) {
        w <- a$fold == k; delong2(b$pred_prob[w], a$pred_prob[w], y[w])$S })
      V1 <- sum(vapply(Sk, function(S) S[1, 1], 0)); V2 <- sum(vapply(Sk, function(S) S[2, 2], 0))
      C  <- sum(vapply(Sk, function(S) S[1, 2], 0))
      data.frame(combination = k, clinical_only = pair_of(k),
                 AUROC_pooled_clinical_only = d$auc[1], AUROC_pooled_with_microbiome = d$auc[2],
                 r_pooled = d$S[1, 2] / sqrt(d$S[1, 1] * d$S[2, 2]),
                 r_folds = C / sqrt(V1 * V2),
                 SE_diff_pooled = sqrt(d$S[1, 1] + d$S[2, 2] - 2 * d$S[1, 2]),
                 SE_diff_foldmean = sqrt(V1 + V2 - 2 * C) / length(Sk))
    }))
    # verificacion de la covarianza de DeLong contra pROC en un par
    k1 <- micro_keys[1]; a <- run[[k1]]$test_predictions; b <- run[[pair_of(k1)]]$test_predictions
    b <- b[match(paste(a$fold, a$id), paste(b$fold, b$id)), ]
    ra <- roc_ptb(a$true_class, a$pred_prob); rb <- roc_ptb(b$true_class, b$pred_prob)
    stopifnot(abs(r_delong$r_pooled[1] - pROC::cov(rb, ra, method = "delong") /
                    sqrt(pROC::var(rb, method = "delong") * pROC::var(ra, method = "delong"))) < 1e-10)
    cat(sprintf("DeLong r between the paired AUROC (12 pairs): pooled median %.3f (range %.3f-%.3f); by fold median %.3f (range %.3f-%.3f)\n",
                stats::median(r_delong$r_pooled), min(r_delong$r_pooled), max(r_delong$r_pooled),
                stats::median(r_delong$r_folds), min(r_delong$r_folds), max(r_delong$r_folds)))
    print(r_delong, row.names = FALSE, digits = 3)
    out$pairs_delong <- r_delong
  }
  # CSV versionado: con <= 1 caso de parto < 32+0, ese desenlace se informa como "<= 1" y la
  # subtarea 2 del DREAM como no evaluable, sin mas conteos (la consola de arriba los muestra)
  low <- windows$outcome == "delivery < 32+0" & windows$cases <= 1
  sc2 <- low & windows$window == "DREAM SC2 window"
  windows_csv <- windows
  windows_csv$cases <- ifelse(low, "<= 1", as.character(windows$cases))
  windows_csv$non_cases[low] <- NA
  windows_csv[sc2, c("participants", "preterm", "term", "eligible_without_sample",
                     "eligible_without_sample_preterm", "samples_in_window")] <- NA
  windows_csv$note <- ifelse(sc2, "not evaluable (<= 1 case)", "")
  out$window_counts <- windows_csv
}

# --- E. tabla por tamano de cohorte --------------------------------------------
cat("\n== E. Cohort-size table ==\n")
size_rows <- bind_rows(
  data.frame(row = "This cohort", n1 = n1, n0 = n0),
  if (!is.null(windows)) {
    w <- windows[windows$window %in% c("DREAM SC1 window", "Landmark week 32"), ]
    data.frame(row = paste(w$window, w$specimens, sep = ", "), n1 = w$cases, n0 = w$non_cases)
  },
  data.frame(row = "This cohort + DREAM validation Hispanic subgroup (prior row)", n1 = 29, n0 = 52),
  bind_rows(lapply(c(81, 130, 200, 300, 500), function(N) data.frame(
    row = c(sprintf("N = %d, cohort proportion 14/43", N), sprintf("N = %d, one third preterm", N)),
    n1 = c(round(N * 14 / 43), round(N / 3)), n0 = c(N - round(N * 14 / 43), N - round(N / 3))))))
# Varianzas por fila:
#   pooled = AUROC agrupado de las n1 + n0 participantes, modelo fijo (Hanley-McNeil)
#   engine = la metrica del motor: media de los AUROC de 5 folds externos estratificados
#            (la composicion real en "This cohort"; n1/5 y n0/5 por fold en el resto),
#            multiplicada por f^2 (reentrenamiento, parte D; supuesto: el mismo f en H1 y
#            en otros tamanos)
r_sets <- c("0.3" = 0.3, "0.5" = 0.5, "0.7" = 0.7, null_median = stats::median(pair_null$r_null))
if (!is.null(r_delong)) r_sets <- c(r_sets, delong_pooled_median = stats::median(r_delong$r_pooled),
                                    delong_folds_median = stats::median(r_delong$r_folds))
engine_var <- function(i) {
  n1i <- size_rows$n1[i]; n0i <- size_rows$n0[i]
  if (size_rows$row[i] == "This cohort") function(a) f_all^2 * hm_var_folds(a, n1k, n0k) else
    function(a) f_all^2 * hm_var(a, n1i / 5, n0i / 5) / 5
}
precision_tab <- bind_rows(lapply(seq_len(nrow(size_rows)), function(i) bind_rows(lapply(A_grid, function(A) {
  n1i <- size_rows$n1[i]; n0i <- size_rows$n0[i]
  data.frame(row = size_rows$row[i], N = n1i + n0i, n1 = n1i, n0 = n0i,
             preterm_share = n1i / (n1i + n0i), A = A,
             SE_pooled = sqrt(hm_var(A, n1i, n0i)), CI_halfwidth_pooled = za * sqrt(hm_var(A, n1i, n0i)),
             SE_engine = sqrt(engine_var(i)(A)), CI_halfwidth_engine = za * sqrt(engine_var(i)(A)),
             delta_prior = delta_prior(A, n1i, n0i))
}))))
delta_tab <- bind_rows(lapply(seq_len(nrow(size_rows)), function(i) bind_rows(lapply(A_grid, function(A) {
  n1i <- size_rows$n1[i]; n0i <- size_rows$n0[i]
  bind_rows(lapply(names(r_sets), function(nm) bind_rows(lapply(c("null", "alternative"), function(cr) data.frame(
    row = size_rows$row[i], N = n1i + n0i, n1 = n1i, n0 = n0i, preterm_share = n1i / (n1i + n0i), A = A,
    r_label = nm, r = r_sets[[nm]], critical = cr,
    delta_pooled = delta_paired(A, pooled_var(n1i, n0i), r_sets[[nm]], cr),
    delta_engine = delta_paired(A, engine_var(i), r_sets[[nm]], cr))))))
}))))
cat("Precision (A = 0.70):\n")
print(precision_tab[precision_tab$A == 0.70, ], row.names = FALSE, digits = 3)
r_main <- if ("delong_folds_median" %in% names(r_sets)) "delong_folds_median" else "0.5"
cat(sprintf("\nDetectable delta, r = %s (%.3f), A = 0.70 and 0.80:\n", r_main, r_sets[[r_main]]))
print(delta_tab[delta_tab$r_label == r_main & delta_tab$A >= 0.70,
                c("row", "N", "n1", "n0", "A", "critical", "delta_pooled", "delta_engine")],
      row.names = FALSE, digits = 3)
out$precision_by_size <- precision_tab
out$detectable_delta_by_size <- delta_tab
out$parameters <- data.frame(
  parameter = c("alpha", "power", "z_alpha", "z_beta", "nsim", "inflation_f", "inflation_f_ci_low",
                "inflation_f_ci_high", "null_sd_test_sampling_permuted", "null_sd_test_sampling_stratified",
                "expected_share_skipped_fold", paste0("r_", names(r_sets))),
  value = c(alpha, power, za, zb, nsim, f_all, f_ci[[1]], f_ci[[2]], sqrt(V_test_perm), sqrt(V_test_strat),
            p_skip_expected, unname(r_sets)))

if (nzchar(results_dir)) {
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(out)) utils::write.csv(out[[nm]], file.path(results_dir, paste0(nm, ".csv")), row.names = FALSE)
  cat("\nWrote", paste0(names(out), ".csv", collapse = ", "), "to", results_dir, "\n")
}
