# Standard testthat harness for R CMD check / CI (wired up in Chat 8, once renv
# is restored and the package installs). Today the suite runs via
#   Rscript -e 'devtools::load_all(); testthat::test_dir("tests/testthat")'
# (see tests/testthat/setup.R), because the package is not yet installed.
library(testthat)
library(ptbpredict)

test_check("ptbpredict")
