# Installation Guide

Step-by-step setup for running the analysis. See the [README](../README.md) for what the pipeline
does and how to run it once installed.

## Requirements

- **R** 4.4.2 or newer ([download](https://cran.r-project.org/))
- **pandoc**, for rendering the analysis report. RStudio and Quarto both bundle one;
  `scripts/run_baseline.R` looks in the usual install locations and sets `RSTUDIO_PANDOC` itself, so
  in most cases there is nothing to do.
- **RStudio** — optional, convenient for interactive work
- **RAM:** 8 GB minimum, 16 GB recommended
- **Disk:** ~500 MB for the repository plus its R library

## Install

```bash
git clone https://github.com/martinruhle/Mexican-PretermBirth-analysis.git
cd Mexican-PretermBirth-analysis
```

```r
install.packages("renv")
renv::restore()
```

`renv::restore()` installs the exact package versions recorded in `renv.lock`, including the
Bioconductor packages (ANCOMBC, phyloseq). The first install is slow, since some of those are
compiled from source. The project `.Rprofile` activates `renv` automatically whenever you open the
project, so no further setup is needed.

## Check that it works

The quickest check is the test suite:

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

For the full suite, including the slow end-to-end no-leakage test (this is what CI runs):

```bash
PTB_RUN_SLOW_TESTS=1 Rscript -e 'testthat::test_dir("tests/testthat")'
```

Then run the pipeline end to end on the synthetic example data (≈12 minutes):

```bash
PTB_PROFILE=example Rscript scripts/run_baseline.R
```

Outputs land in `analysis/_output/`.

## Troubleshooting

**`renv::restore()` stops while building a Bioconductor package.** These build from source and need
system libraries. On Debian/Ubuntu:

```bash
sudo apt-get install -y libglpk-dev libxml2-dev zlib1g-dev libbz2-dev liblzma-dev \
  libgsl-dev libcurl4-openssl-dev libssl-dev libudunits2-dev libfontconfig1-dev \
  libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```

The [CI workflow](../.github/workflows/check-code.yml) installs exactly this set, so it is a
reliable reference for a clean Linux machine.

**Rendering stops with a pandoc message.** Install RStudio or Quarto, or point `RSTUDIO_PANDOC` at
an existing pandoc directory.

**`renv::status()` reports the project is out of sync.** Expected: `config` is declared under
`Suggests` and is deliberately not installed, and two packages (`microbiome`, `Rtsne`) are pinned by
hand. See the notes in [`CONTRIBUTING.md`](../CONTRIBUTING.md) before running `renv::snapshot()`.

**Something else.** Open a [GitHub issue](https://github.com/martinruhle/Mexican-PretermBirth-analysis/issues)
with your `sessionInfo()` and the full console output.

---

**Last updated:** 2026-08-10
