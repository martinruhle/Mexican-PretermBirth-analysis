# microbiome_columns() + split_domains(): the microbiome/clinical split is driven
# by NAME (role == "microbiome" in the dictionary), never by column position or a
# `^g__` regex. This is the invariant that lets a foreign dataset with a different
# column order run unchanged.

test_that("microbiome_columns returns the dictionary's microbiome variables, in dictionary order", {
  cfg <- make_cfg(); dict <- make_dict()
  expect_identical(microbiome_columns(cfg, dict = dict),
                   dict$variable[dict$role == "microbiome"])
})

test_that("split_domains selects taxa by name and carries the keys into both domains", {
  cfg <- make_cfg(); dict <- make_dict(); m <- make_matrix()
  taxa <- dict$variable[dict$role == "microbiome"]

  sp <- split_domains(m, cfg, dict = dict)

  # microbiome domain = taxa + the sample_id / subject_id keys
  expect_true(all(taxa %in% names(sp$microbiome)))
  expect_true(all(c("index", "id") %in% names(sp$microbiome)))

  # taxa are NOT duplicated into the clinical domain ...
  expect_false(any(taxa %in% names(sp$clinical)))
  # ... but the keys live in BOTH, and outcome / GA / clinical live in clinical
  expect_true(all(c("index", "id") %in% names(sp$clinical)))
  expect_true(all(c("preterm", "sdg_parto", "edad_cronologicamujer") %in% names(sp$clinical)))
})

test_that("split_domains is positional-order-independent (name-based, not by index)", {
  cfg <- make_cfg(); dict <- make_dict(); m <- make_matrix()
  taxa <- dict$variable[dict$role == "microbiome"]

  set.seed(1)
  m_shuffled <- m[, sample(names(m))]                 # scramble the column order

  sp      <- split_domains(m,          cfg, dict = dict)
  sp_shuf <- split_domains(m_shuffled, cfg, dict = dict)

  expect_setequal(setdiff(names(sp_shuf$microbiome), c("index", "id")),
                  setdiff(names(sp$microbiome),      c("index", "id")))
  expect_setequal(setdiff(names(sp_shuf$microbiome), c("index", "id")), taxa)
})

test_that("split_domains preserves hyphenated and family-prefixed taxa names exactly", {
  cfg <- make_cfg(); dict <- make_dict(); m <- make_matrix()
  sp <- split_domains(m, cfg, dict = dict)
  expect_true("Escherichia-Shigella" %in% names(sp$microbiome))   # hyphen not mangled to a dot
  expect_true("f__Rhizobiaceae"      %in% names(sp$microbiome))   # family prefix intact
})
