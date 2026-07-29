# ============================================================================
# Data input / output — config-driven, name-based (Chat 3)
# ----------------------------------------------------------------------------
# Replaces the reproducibility-breaking parts of the .Rmd `load_data` chunk:
#   * setwd("C:/Users/marti/Documents/Datos_mexicanos/") + implicit relative reads
#     -> paths come from config, resolved with here::here() (root-relative).
#   * positional domain split  micro <- data[, 1:99];  clin <- data[, 98:167]
#     (columns 98-99 = `index`/`id` were duplicated into BOTH blocks)
#     -> split_domains() selects by NAME using the role = "microbiome" rows of
#        config/data_dictionary.csv. No positions, no `^g__` regex (that prefix
#        does not exist in these data).
#   * the ANCOM branch's hardcoded absolute abs-matrix path + positional indexing
#     (read.csv(row.names = 98), [, 1:97], [, 98:166]) -> load_abs_matrix(), by name.
#
# Pure functions: everything enters by argument, nothing reads globals, no setwd,
# no absolute paths. Depends on: here, yaml (or config), readr, utils.
# ============================================================================

#' Read the ptbpredict configuration profile
#'
#' Prefers `config::get()` when the `config` package is installed (the target
#' state after renv is restored in Chat 8); otherwise falls back to parsing the
#' YAML directly with `yaml`, selecting the requested profile and, for non-default
#' profiles, shallow-merging it over `default` (mirroring `config`'s inheritance).
#'
#' @param profile Configuration profile name. Defaults to the `PTB_PROFILE`
#'   environment variable, or `"default"`.
#' @param file Path to the config YAML. Defaults to `config/config.yml` under the
#'   project root (resolved with [here::here()]).
#'
#' @return A named list with the resolved configuration.
#' @export
read_config <- function(profile = Sys.getenv("PTB_PROFILE", "default"),
                        file = here::here("config", "config.yml")) {
  if (requireNamespace("config", quietly = TRUE)) {
    return(config::get(config = profile, file = file))
  }
  yml <- yaml::yaml.load_file(file)
  if (!profile %in% names(yml)) {
    stop("Perfil '", profile, "' no existe en ", file, call. = FALSE)
  }
  cfg <- yml[[profile]]
  if (profile != "default" && "default" %in% names(yml)) {
    cfg <- utils::modifyList(yml[["default"]], cfg)
  }
  cfg
}

#' Read the data dictionary (contract) declared in the config
#'
#' @param cfg Configuration list from [read_config()].
#'
#' @return A data.frame with (at least) `variable` and `role` columns.
#' @export
read_data_dictionary <- function(cfg) {
  dict_path <- cfg$data$data_dictionary %||% "config/data_dictionary.csv"
  readr::read_csv(here::here(dict_path), show_col_types = FALSE,
                  progress = FALSE)
}

#' Names of the microbiome (taxa) columns, from the data dictionary
#'
#' The single source of truth for which columns are microbiome vs clinical is the
#' `role == "microbiome"` rows of `config/data_dictionary.csv` — never a `^g__`
#' regex (no such prefix here) nor positional indices.
#'
#' @param cfg Configuration list from [read_config()].
#' @param dict Optional pre-loaded dictionary (avoids re-reading).
#'
#' @return Character vector of taxa column names, in dictionary order.
#' @export
microbiome_columns <- function(cfg, dict = NULL) {
  if (is.null(dict)) dict <- read_data_dictionary(cfg)
  as.character(dict$variable[dict$role == "microbiome"])
}

#' Load the input datasets declared in the config
#'
#' Reads the wide matrix (`matrix_path`) and the longitudinal metadata
#' (`metadata_path`) with `readr::read_csv` (preserving exact column names — e.g.
#' the hyphen in `Escherichia-Shigella` — unlike base `read.csv`). The optional
#' participant table (`participant_path`) is read only when a path is configured;
#' it is `null` by default because it is loaded-but-unused in the current pipeline.
#'
#' @param cfg Configuration list from [read_config()].
#'
#' @return A named list with `matrix`, `metadata`, and `participant` (the last is
#'   `NULL` when `participant_path` is `null`).
#' @export
load_dataset <- function(cfg) {
  matrix_path <- here::here(cfg$data$matrix_path)
  meta_path   <- here::here(cfg$data$metadata_path)

  wide_matrix <- readr::read_csv(matrix_path, show_col_types = FALSE,
                                 progress = FALSE)
  metadata    <- readr::read_csv(meta_path, show_col_types = FALSE,
                                 progress = FALSE)

  participant <- NULL
  if (!is.null(cfg$data$participant_path)) {
    participant <- readr::read_csv(here::here(cfg$data$participant_path),
                                   show_col_types = FALSE, progress = FALSE)
  }

  list(matrix = wide_matrix, metadata = metadata, participant = participant)
}

#' Split a wide matrix into microbiome and clinical domains, by name
#'
#' Microbiome domain = the `role == "microbiome"` taxa columns plus the sample_id
#' and subject_id keys. Clinical domain = every remaining (non-taxa) column, which
#' naturally carries the keys, the outcome and gestational-age columns, and all
#' clinical variables. Column order within each domain follows the original matrix,
#' so this reproduces the legacy `data[, 1:99]` / `data[, 98:167]` objects exactly
#' — minus the bug where the shared keys `index`/`id` (matrix columns 98-99) landed
#' in BOTH blocks via the overlapping positional ranges.
#'
#' @param data The wide matrix (e.g. `load_dataset(cfg)$matrix`).
#' @param cfg Configuration list from [read_config()].
#' @param dict Optional pre-loaded dictionary (avoids re-reading).
#'
#' @return A named list with `microbiome` and `clinical` data.frames.
#' @export
split_domains <- function(data, cfg, dict = NULL) {
  sample_id  <- cfg$columns$sample_id
  subject_id <- cfg$columns$subject_id
  micro_cols <- microbiome_columns(cfg, dict = dict)

  missing_keys <- setdiff(c(sample_id, subject_id), names(data))
  if (length(missing_keys) > 0) {
    stop("Columnas clave ausentes en la matriz: ",
         paste(missing_keys, collapse = ", "), call. = FALSE)
  }
  # intersect() keeps names(data) order -> original matrix (positional) order.
  taxa_present <- intersect(names(data), micro_cols)

  microbiome <- data[, c(taxa_present, sample_id, subject_id), drop = FALSE]
  clinical   <- data[, setdiff(names(data), taxa_present), drop = FALSE]

  list(microbiome = microbiome, clinical = clinical)
}

#' Load the absolute-count matrix for ANCOM-BC2, split by name
#'
#' Replaces the ANCOM branch's hardcoded absolute path + positional indexing
#' (`read.csv(abs_file, row.names = 98)`, `[, 1:97]`, `[, 98:166]`). The sample_id
#' column becomes the row names (as `row.names = 98` did, but selected by NAME),
#' and taxa vs metadata are separated with the dictionary's `role == "microbiome"`
#' list. Read with `check.names = FALSE`, so taxa names are preserved exactly
#' (base `read.csv` would mangle `Escherichia-Shigella` -> `Escherichia.Shigella`,
#' which silently dropped that taxon from the legacy `intersect(genera_clean, ...)`).
#'
#' @param cfg Configuration list from [read_config()].
#' @param dict Optional pre-loaded dictionary (avoids re-reading).
#'
#' @return A named list with `otu` (samples x taxa, absolute counts; row names =
#'   sample_id) and `meta` (subject_id + clinical columns; same row order/names).
#' @export
load_abs_matrix <- function(cfg, dict = NULL) {
  abs_path  <- here::here(cfg$data$abs_matrix_path)
  sample_id <- cfg$columns$sample_id

  df <- utils::read.csv(abs_path, row.names = sample_id, check.names = FALSE,
                        stringsAsFactors = FALSE)

  micro_cols   <- microbiome_columns(cfg, dict = dict)
  taxa_present <- intersect(colnames(df), micro_cols)  # matrix order preserved
  otu  <- df[, taxa_present, drop = FALSE]
  meta <- df[, setdiff(colnames(df), taxa_present), drop = FALSE]

  list(otu = otu, meta = meta)
}

# Null-coalescing helper (base R has none). Local to R/ so functions above are
# self-contained without importing rlang.
`%||%` <- function(x, y) if (is.null(x)) y else x
