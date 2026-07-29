# ============================================================================
# Data contract validation (Chat 3, playbook §4)
# ----------------------------------------------------------------------------
# Turns "breaks silently on someone else's data" into "stops with a message that
# says exactly what to fix". Replaces the positional-indexing assumptions of the
# old load_data chunk with explicit, name-based, config-driven checks:
#   * required columns present (subject_id / sample_id / outcome / gestational_age)
#   * outcome is binary 0/1
#   * every role = "microbiome" taxon in the dictionary exists in the matrix
#   * taxa columns live in the relative-abundance range [0, 1]  (this is the SAME
#     check that governs the CLR pseudocount decision — counts would blow the range)
#   * sample_id / subject_id agree between the wide matrix and the metadata table
#
# Pure function; reads nothing from the global environment.
# ============================================================================

#' Validate the input datasets against the config data contract
#'
#' Fails hard and clearly when the matrix or metadata do not satisfy the contract
#' declared in `cfg` + `config/data_dictionary.csv`. Intended to be called right
#' after [load_dataset()], before any domain split or modelling.
#'
#' @param matrix The wide input matrix (`load_dataset(cfg)$matrix`).
#' @param metadata The longitudinal metadata table
#'   (`load_dataset(cfg)$metadata`); pass `NULL` to skip the cross-table checks.
#' @param cfg Configuration list from [read_config()].
#' @param dict Optional pre-loaded dictionary (avoids re-reading).
#' @param tol Numeric tolerance for the taxa range check. Default `1e-6`.
#'
#' @return `matrix`, invisibly, if all checks pass; otherwise an informative error.
#' @export
validate_input_data <- function(matrix, metadata, cfg, dict = NULL, tol = 1e-6) {

  subject_id <- cfg$columns$subject_id
  sample_id  <- cfg$columns$sample_id
  outcome    <- cfg$columns$outcome
  gest_age   <- cfg$columns$gestational_age

  # -- 1. required columns present in the matrix ----------------------------
  required <- c(subject_id, sample_id, outcome, gest_age)
  missing  <- setdiff(required, names(matrix))
  if (length(missing) > 0) {
    stop("Columnas requeridas ausentes en la matriz: ",
         paste(missing, collapse = ", "),
         "\nRevisa config/config.yml -> columns.", call. = FALSE)
  }

  # -- 2. outcome is binary 0/1 ---------------------------------------------
  out_levels <- as.character(unique(stats::na.omit(matrix[[outcome]])))
  if (!all(out_levels %in% c("0", "1"))) {
    stop("El outcome '", outcome, "' debe ser binario 0/1; encontrado: ",
         paste(out_levels, collapse = ", "), call. = FALSE)
  }

  # -- 3. every dictionary taxon exists in the matrix -----------------------
  if (is.null(dict)) dict <- read_data_dictionary(cfg)
  micro_cols <- microbiome_columns(cfg, dict = dict)
  if (length(micro_cols) == 0) {
    stop("El diccionario no declara ninguna columna con role = 'microbiome'.",
         call. = FALSE)
  }
  missing_micro <- setdiff(micro_cols, names(matrix))
  if (length(missing_micro) > 0) {
    stop("Taxa del diccionario ausentes en la matriz: ",
         paste(utils::head(missing_micro, 10), collapse = ", "),
         if (length(missing_micro) > 10) " ..." else "",
         "\n(", length(missing_micro), " en total). ",
         "Revisa config/data_dictionary.csv (role = microbiome).", call. = FALSE)
  }

  # -- 4. taxa columns are relative abundances in [0, 1] --------------------
  taxa_mat <- as.matrix(matrix[, micro_cols, drop = FALSE])
  storage.mode(taxa_mat) <- "numeric"
  rng <- range(taxa_mat, na.rm = TRUE)
  if (rng[1] < -tol || rng[2] > 1 + tol) {
    stop(sprintf(paste0(
      "Las columnas de taxa deben ser abundancias relativas en [0, 1]; ",
      "rango observado = [%.4g, %.4g].\n",
      "Si son porcentajes (0-100) o conteos, conviértelas antes: este es el ",
      "mismo rango que fija el pseudocount del CLR."), rng[1], rng[2]),
      call. = FALSE)
  }

  # -- 5. cross-table key agreement (matrix <-> metadata) -------------------
  if (!is.null(metadata)) {
    meta_missing <- setdiff(c(sample_id, subject_id), names(metadata))
    if (length(meta_missing) > 0) {
      stop("La metadata no tiene las columnas clave: ",
           paste(meta_missing, collapse = ", "),
           "\nRevisa config/config.yml -> columns.", call. = FALSE)
    }
    m_samp <- unique(as.character(matrix[[sample_id]]))
    d_samp <- unique(as.character(metadata[[sample_id]]))
    only_matrix <- setdiff(m_samp, d_samp)
    only_meta   <- setdiff(d_samp, m_samp)
    if (length(only_matrix) > 0 || length(only_meta) > 0) {
      stop("Los sample_id ('", sample_id, "') no casan entre matriz y metadata.\n",
           if (length(only_matrix) > 0)
             paste0("  En la matriz pero no en la metadata (",
                    length(only_matrix), "): ",
                    paste(utils::head(only_matrix, 8), collapse = ", "), "\n"),
           if (length(only_meta) > 0)
             paste0("  En la metadata pero no en la matriz (",
                    length(only_meta), "): ",
                    paste(utils::head(only_meta, 8), collapse = ", "), "\n"),
           call. = FALSE)
    }
    m_subj <- unique(as.character(matrix[[subject_id]]))
    d_subj <- unique(as.character(metadata[[subject_id]]))
    subj_only_matrix <- setdiff(m_subj, d_subj)
    subj_only_meta   <- setdiff(d_subj, m_subj)
    if (length(subj_only_matrix) > 0 || length(subj_only_meta) > 0) {
      stop("Los subject_id ('", subject_id, "') no casan entre matriz y metadata.\n",
           if (length(subj_only_matrix) > 0)
             paste0("  En la matriz pero no en la metadata (",
                    length(subj_only_matrix), "): ",
                    paste(utils::head(subj_only_matrix, 8), collapse = ", "), "\n"),
           if (length(subj_only_meta) > 0)
             paste0("  En la metadata pero no en la matriz (",
                    length(subj_only_meta), "): ",
                    paste(utils::head(subj_only_meta, 8), collapse = ", "), "\n"),
           call. = FALSE)
    }
  }

  message(sprintf(
    "✓ Contrato de datos validado: %d taxa, %d muestras, %d sujetos%s.",
    length(micro_cols), nrow(matrix),
    length(unique(matrix[[subject_id]])),
    if (!is.null(metadata)) " (sample_id/subject_id casan con la metadata)" else ""))

  invisible(matrix)
}
