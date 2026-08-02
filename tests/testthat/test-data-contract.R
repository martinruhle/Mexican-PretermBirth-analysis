# validate_input_data(): the 5 contract checks (playbook §4).
# Error substrings are copied verbatim from R/data_contract.R — they document
# the messages the function actually raises, not aspirational ones.

test_that("a valid fixture passes: emits the OK message and returns the matrix invisibly", {
  cfg <- make_cfg(); dict <- make_dict(); m <- make_matrix(); meta <- make_metadata()

  expect_message(validate_input_data(m, meta, cfg, dict = dict),
                 "Contrato de datos validado")

  suppressMessages(vis <- withVisible(validate_input_data(m, meta, cfg, dict = dict)))
  expect_false(vis$visible)          # returns invisibly
  expect_equal(vis$value, m)         # returns the matrix unchanged
})

test_that("check 1: missing required column is reported", {
  cfg <- make_cfg(); dict <- make_dict(); meta <- make_metadata()
  m <- make_matrix(); m$sdg_parto <- NULL          # drop gestational_age column
  expect_error(validate_input_data(m, meta, cfg, dict = dict),
               "Columnas requeridas ausentes en la matriz")
})

test_that("check 2: non-binary outcome is reported", {
  cfg <- make_cfg(); dict <- make_dict(); meta <- make_metadata()
  m <- make_matrix(); m$preterm[1] <- 2            # 0/1/2 is not binary
  expect_error(validate_input_data(m, meta, cfg, dict = dict),
               "debe ser binario 0/1")
})

test_that("check 3: dictionary taxon absent from the matrix is reported", {
  cfg <- make_cfg(); meta <- make_metadata(); m <- make_matrix()
  dict <- rbind(make_dict(),
                data.frame(variable = "Bogus_taxon", role = "microbiome",
                           stringsAsFactors = FALSE))
  expect_error(validate_input_data(m, meta, cfg, dict = dict),
               "Taxa del diccionario ausentes en la matriz")
})

test_that("check 4: taxa outside the relative-abundance range [0, 1] are reported", {
  cfg <- make_cfg(); dict <- make_dict(); meta <- make_metadata()
  m <- make_matrix(); m[["Lactobacillus"]][1] <- 5   # count-scale value, not a proportion
  expect_error(validate_input_data(m, meta, cfg, dict = dict),
               "deben ser abundancias relativas")
})

test_that("check 5: sample_id keys that disagree with the metadata are reported", {
  cfg <- make_cfg(); dict <- make_dict(); m <- make_matrix()
  meta <- make_metadata(); meta$index <- paste0(meta$index, "X")   # break sample_id only
  expect_error(validate_input_data(m, meta, cfg, dict = dict),
               "sample_id.*no casan entre matriz y metadata")
})

test_that("check 5: subject_id keys that disagree with the metadata are reported", {
  cfg <- make_cfg(); dict <- make_dict(); m <- make_matrix()
  meta <- make_metadata(); meta$id <- paste0(meta$id, "X")         # break subject_id only
  expect_error(validate_input_data(m, meta, cfg, dict = dict),
               "subject_id.*no casan entre matriz y metadata")
})
