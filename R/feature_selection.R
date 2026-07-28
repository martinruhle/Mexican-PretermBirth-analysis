# ============================================================================
# ANCOM-BC2 differential-abundance feature selection, fitted WITHIN a fold
# ----------------------------------------------------------------------------
# Extracted verbatim from the chunk `nested_cv_helper_functions` of
# analysis/integrated_preterm_prediction_workflow.Rmd (Chat 2). Logic unchanged;
# only roxygen docs were added.
#
# Leakage invariant: this is fit ONLY on the outer-training subjects of the
# current fold (`outer_train_subjects`); the taxa it selects are then applied to
# the held-out data. It must never see the outer-test subjects.
#
# Depends on the analysis packages being attached (phyloseq::sample_data,
# phyloseq::otu_table, phyloseq::phyloseq, phyloseq::merge_samples,
# ANCOMBC::ancombc2, and dplyr), exactly as in the source .Rmd.
# ============================================================================

#' Run ANCOM-BC2 on a fold's training subjects and return significant taxa
#'
#' Filters the OTU/metadata tables to the current fold's outer-training subjects,
#' aggregates repeated (longitudinal) samples to one row per subject, and runs
#' ANCOM-BC2 adjusting for maternal age and pre-gestational BMI. Returns the taxa
#' differentially abundant by preterm status (p < 0.10), excluding the likely
#' artefact `o__Chloroplast`.
#'
#' @param otu_data Data.frame of absolute abundances, samples in rows and taxa in
#'   columns, row order aligned with `meta_data`.
#' @param meta_data Data.frame of sample metadata aligned to `otu_data`, including
#'   `id`, `preterm`, `edad_cronologicamujer` and `imc_pregestacional`.
#' @param outer_train_subjects Character/numeric vector of subject `id`s that make
#'   up the outer-training set for this fold.
#'
#' @return Character vector of selected taxa names (possibly empty), ordered by
#'   ascending p-value.
#'
#' @export
run_ancombc_on_fold <- function(otu_data, meta_data, outer_train_subjects) {

  # Suppress ANCOM-BC2 verbose output
  suppressMessages({
    suppressWarnings({

      # Filter to training subjects only
      training_indices <- which(meta_data$id %in% outer_train_subjects)
      otu_table_train <- otu_data[training_indices, ]
      meta_data_train <- meta_data[training_indices, ]

      # Format metadata
      meta_data_train$preterm <- as.factor(meta_data_train$preterm)
      meta_data_train$id <- as.factor(meta_data_train$id)
      meta_data_train$edad_cronologicamujer <- as.numeric(meta_data_train$edad_cronologicamujer)
      meta_data_train$imc_pregestacional <- as.numeric(meta_data_train$imc_pregestacional)

      # Impute missing values (median for numeric)
      if(any(is.na(meta_data_train$edad_cronologicamujer))) {
        median_age <- median(meta_data_train$edad_cronologicamujer, na.rm = TRUE)
        meta_data_train$edad_cronologicamujer[is.na(meta_data_train$edad_cronologicamujer)] <- median_age
      }

      if(any(is.na(meta_data_train$imc_pregestacional))) {
        median_bmi <- median(meta_data_train$imc_pregestacional, na.rm = TRUE)
        meta_data_train$imc_pregestacional[is.na(meta_data_train$imc_pregestacional)] <- median_bmi
      }

      # Create phyloseq object
      meta_data_ps <- sample_data(meta_data_train)
      otu_table_t <- t(otu_table_train)
      otu_table_ps <- otu_table(otu_table_t, taxa_are_rows = TRUE)

      pseq_train <- phyloseq(otu_table_ps, meta_data_ps)

      # Aggregate by subject (account for repeated measures)
      pseq_agregado <- merge_samples(pseq_train, group = "id")

      # Reconstruct metadata (one row per subject)
      meta_sujetos <- as(sample_data(pseq_train), "data.frame") %>%
        distinct(id, .keep_all = TRUE)
      rownames(meta_sujetos) <- meta_sujetos$id
      sample_data(pseq_agregado) <- meta_sujetos

      # Run ANCOM-BC2
      formula_ancombc <- "preterm + edad_cronologicamujer + imc_pregestacional"

      ancombc_output <- ancombc2(
        data = pseq_agregado,
        fix_formula = formula_ancombc,
        p_adj_method = "BH",
        lib_cut = 1000,
        pseudo = 0,
        pseudo_sens = TRUE,
        prv_cut = 0.05,
        group = "preterm",
        struc_zero = FALSE,
        neg_lb = FALSE,
        alpha = 0.10,
        verbose = FALSE
      )

      # Extract results
      ancombc_res <- ancombc_output$res

      # Select significant taxa (p < 0.10)
      # EXCLUDE o__Chloroplast (likely artifact)
      significant_taxa <- ancombc_res %>%
        filter(
          p_preterm1 < 0.10,
          taxon != "o__Chloroplast"
        ) %>%
        arrange(p_preterm1) %>%
        pull(taxon)

      return(as.character(significant_taxa))

    })
  })
}
