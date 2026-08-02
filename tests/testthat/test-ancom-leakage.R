# Component-level no-leakage (#7): run_ancombc_on_fold() selects taxa using ONLY
# the outer-training subjects. We prove it by perturbing the held-out (non-train)
# subjects' counts arbitrarily and checking the selected taxa are unchanged: if
# any held-out sample influenced the result, the two runs would differ.
#
# Uses the committed synthetic example (abs-count matrix). One ANCOM-BC2 run is
# ~13s here, so two runs ~25-30s; gated only on package availability.

test_that("run_ancombc_on_fold uses only the outer-training subjects (no held-out leakage)", {
  skip_if_not_installed("ANCOMBC")
  skip_if_not_installed("phyloseq")

  cfg <- read_config("example")
  abs <- load_abs_matrix(cfg)                       # $otu (samples x taxa), $meta (id + clinical)

  # subject-level labels straight from the abs metadata
  lab  <- tapply(abs$meta$preterm, abs$meta$id, function(x) as.character(x[1]))
  ptb  <- names(lab)[lab == "1"]
  term <- names(lab)[lab == "0"]

  set.seed(7)
  held_out       <- c(sample(ptb, 2), sample(term, 4))   # excluded from training
  train_subjects <- setdiff(names(lab), held_out)         # both classes present

  invisible(capture.output(
    taxa1 <- run_ancombc_on_fold(abs$otu, abs$meta, train_subjects)))

  # Perturb ONLY the held-out subjects' counts, then re-run on the same train set.
  otu_perturbed <- abs$otu
  rows_held <- which(abs$meta$id %in% held_out)
  otu_perturbed[rows_held, ] <- otu_perturbed[rows_held, ] * 1000L + 777L
  invisible(capture.output(
    taxa2 <- run_ancombc_on_fold(otu_perturbed, abs$meta, train_subjects)))

  # Selection is a function of the training subjects alone -> byte-identical.
  expect_identical(taxa1, taxa2)
})
