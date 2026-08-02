# Light integration (#6): the committed synthetic example must load, validate and
# split end-to-end. This is a SMOKE test — it proves the load/contract pipeline
# RUNS on real files (I/O included), not that it is leakage-free. No-leakage is
# covered by the component tests (test-clr / test-threshold / test-ancom-leakage)
# and end-to-end by test-leakage-permutation.

test_that("the example profile loads, validates and splits into 97 taxa + keys", {
  cfg <- read_config("example")

  ds <- load_dataset(cfg)
  expect_s3_class(ds$matrix, "data.frame")
  expect_s3_class(ds$metadata, "data.frame")

  # data contract passes on the committed example (reads the dictionary from disk)
  expect_message(validate_input_data(ds$matrix, ds$metadata, cfg),
                 "Contrato de datos validado")

  sp <- split_domains(ds$matrix, cfg)
  taxa <- setdiff(names(sp$microbiome), c("index", "id"))
  expect_length(taxa, 97)
  expect_true(all(c("index", "id") %in% names(sp$microbiome)))
})

test_that("load_abs_matrix preserves exact taxa names on the example (no read.csv mangling)", {
  cfg <- read_config("example")
  abs <- load_abs_matrix(cfg)

  # the two names that base read.csv would corrupt / that were silently dropped
  # by the legacy positional ANCOM branch
  expect_true("Escherichia-Shigella"     %in% colnames(abs$otu))
  expect_true("f__Bifidobacteriaceae.1"  %in% colnames(abs$otu))
  expect_length(intersect(microbiome_columns(cfg), colnames(abs$otu)), 97)
})
