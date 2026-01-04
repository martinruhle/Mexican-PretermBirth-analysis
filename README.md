# Machine Learning Models for Preterm Birth Prediction Using Vaginal Microbiome Profiles in a Mexican Cohort

[![R Code Verification](https://github.com/martinruhle/Mexican-PretermBirth-analysis/actions/workflows/check-code.yml/badge.svg)](https://github.com/martinruhle/Mexican-PretermBirth-analysis/actions/workflows/check-code.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Status:** Manuscript submitted to *International Journal of Molecular Sciences* - Special Issue on Pregnancy Complications

---

## Citation

```
Ruhle, M., et al. (2025). Machine Learning Models for Preterm Birth Prediction 
Using Vaginal Microbiome Profiles in a Mexican Cohort. 
International Journal of Molecular Sciences. [Manuscript submitted]
```

---

## Overview

This repository contains the complete analysis code, documentation, and public data for a machine learning study investigating preterm birth prediction using vaginal microbiome profiles in a Mexican cohort. This represents the **first machine learning-based preterm birth prediction study conducted specifically in a Mexican population**.

### Background

Preterm birth (PTB) affects approximately 11% of pregnancies globally and is the leading cause of infant mortality and morbidity. The vaginal microbiome has emerged as a modifiable risk factor, with distinct microbial signatures associated with preterm delivery. However, most research has focused on populations of European or African ancestry, creating a critical gap in understanding PTB risk factors among Latin American women.

### Study Design

- **Population:** 43 pregnant women from Instituto Nacional de Perinatología, Mexico City
- **Samples:** 110 longitudinal vaginal microbiome samples (16S rRNA sequencing)
- **Outcome:** Preterm birth (<37 weeks gestation); 14 PTB cases (32.6%)
- **Approach:** Machine learning models (Random Forest, Elastic Net) with rigorous nested cross-validation
- **Data:** Integration of microbiome profiles, clinical variables, and nutritional data

### Key Findings

- **Best performance:** Random Forest with microbiome + clinical data achieved AUROC 0.849 (±0.091)
- **Microbial associations:** Differential abundance analysis identified 14 genera with PTB associations
- **Feature importance:** Nutritional variables (folic acid, vitamins) and specific taxa showed predictive value
- **Clinical potential:** Models demonstrate feasibility of microbiome-based PTB risk assessment in Mexican women

---

## Table of Contents

- [Repository Structure](#repository-structure)
- [System Requirements](#system-requirements)
- [Quick Start](#quick-start)
- [Data](#data)
- [Methodology](#methodology)
- [Results](#results)
- [Reproducibility](#reproducibility)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)
- [Acknowledgments](#acknowledgments)

---

## Repository Structure

```
Mexican-PretermBirth-analysis/
│
├── README.md                          # This file
├── LICENSE                            # MIT License
├── .gitignore                         # Git ignore rules
├── renv.lock                          # R package dependencies (reproducibility)
│
├── .github/
│   └── workflows/
│       └── check-code.yml             # GitHub Actions CI/CD
│
├── code/
│   └── integrated_preterm_prediction_workflow.Rmd    # Main analysis (R Markdown)
│
├── data/
│   ├── metadata/
│   │   ├── diccionario_variables_completo.csv        # Variable dictionary
│   │   ├── metadata_eugenia_long.csv                 # Longitudinal clinical data
│   │   └── participant_data_clean.csv                # Participant-level data
│   ├── example/
│   │   └── abundance_example.csv                     # Example microbiome data structure
│   └── README_DATA.md                                # Data documentation
│
├── results/
│   ├── integrated_preterm_prediction_workflow.pdf    # Rendered analysis report
│   └── figures/
│       └── [Generated figures from analysis]
│
└── docs/
    ├── INSTALL.md                     # Detailed installation instructions
    └── DATA_ACCESS.md                 # Data access request information
```

---

## System Requirements

### Software
- **R:** Version 4.4.2 or higher ([download](https://cran.r-project.org/))
- **RStudio:** Recommended for interactive analysis ([download](https://posit.co/download/rstudio-desktop/))
- **Git:** For repository cloning ([download](https://git-scm.com/downloads))

### Hardware
- **RAM:** 8 GB minimum, 16 GB recommended
- **Storage:** ~500 MB for repository + dependencies
- **OS:** Windows 10+, macOS 10.15+, or Linux (Ubuntu 20.04+)

### Key R Packages
The analysis depends on several R packages (full list in `renv.lock`):
- **Core ML:** tidymodels, ranger, glmnet
- **Microbiome:** ANCOMBC, compositions, vegan
- **Tidyverse:** dplyr, ggplot2, tidyr, purrr
- **Performance:** pROC, yardstick, probably

---

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/martinruhle/Mexican-PretermBirth-analysis.git
cd Mexican-PretermBirth-analysis
```

### 2. Install Dependencies

**Option A: Using renv (Recommended)**
```r
# In R console
install.packages("renv")
renv::restore()  # Installs exact versions used in original analysis
```

**Option B: Manual installation**
```r
# Install key packages
install.packages(c("tidyverse", "tidymodels", "ranger", "glmnet"))

# Bioconductor packages
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("ANCOMBC")
```

For detailed installation instructions and troubleshooting, see **[docs/INSTALL.md](docs/INSTALL.md)**.

### 3. View the Analysis

**Option 1:** Open the rendered PDF:
```
results/integrated_preterm_prediction_workflow.pdf
```

**Option 2:** Run the R Markdown notebook (requires full microbiome data):
```r
# Open in RStudio
file.edit("code/integrated_preterm_prediction_workflow.Rmd")

# Or render from command line
rmarkdown::render("code/integrated_preterm_prediction_workflow.Rmd")
```

**Note:** The full analysis requires restricted microbiome abundance data not included in this repository. The example data structure is provided in `data/example/`. See [Data Access](#data) for information on requesting complete data.

---

## Data

### Publicly Available Data

This repository includes **publicly available metadata** (n=43 participants, 110 samples):

1. **diccionario_variables_completo.csv:** Variable dictionary (Spanish/English)
2. **metadata_eugenia_long.csv:** Longitudinal clinical measurements
   - Demographics, anthropometrics, nutritional intake
   - Clinical complications, laboratory results
   - Pregnancy outcomes
3. **participant_data_clean.csv:** Participant-level summary data

### Restricted Data (Not Included)

**Microbiome abundance tables** are not included due to:
- Ethical review board restrictions
- Participant privacy considerations  
- Data sharing agreements with Instituto Nacional de Perinatología

An **example data structure** is provided in `data/example/abundance_example.csv` to demonstrate the format.

### Requesting Full Data Access

Researchers interested in accessing the complete dataset (including microbiome abundances) should follow the process outlined in **[docs/DATA_ACCESS.md](docs/DATA_ACCESS.md)**.

**Requirements:**
- Institutional affiliation
- Research ethics approval
- Data use agreement
- Scientific justification

**Expected timeline:** 2-4 weeks for review

For detailed data documentation, see **[data/README_DATA.md](data/README_DATA.md)**.

---

## Methodology

### Machine Learning Approach

**Models:**
- **Random Forest:** Ensemble of 500 decision trees (ranger package)
- **Elastic Net:** Logistic regression with L1+L2 regularization (glmnet package)

**Feature Sets:**
1. Clinical variables only
2. Microbiome only (CLR-transformed genus abundances)
3. Combined microbiome + clinical

**Validation Strategy:**
- **Nested 5-fold cross-validation** to prevent data leakage
- **Outer loop:** Performance estimation on held-out test folds
- **Inner loop:** Threshold optimization on independent validation set
- **Stratified sampling** to maintain outcome prevalence

### Data Processing

**Microbiome:**
- Centered log-ratio (CLR) transformation for compositional data
- Differential abundance analysis (ANCOM-BC2) within CV folds
- Rare taxa filtering (prevalence thresholds)

**Clinical:**
- Univariate feature selection (p<0.20 threshold)
- Completeness filtering (≥80% data availability)
- Collinearity management (correlation threshold |r|>0.95)

**Missing Data:**
- Multiple imputation for sporadic missingness
- Case weighting for longitudinal samples (inverse number of visits)

### Performance Metrics

- **AUROC:** Primary metric (discrimination)
- **PRAUC:** Secondary metric (imbalanced data sensitivity)
- **Sensitivity/Specificity:** Clinical interpretation at optimized threshold
- **Threshold optimization:** Youden's Index on independent validation data

For complete methodological details, see the manuscript and `code/integrated_preterm_prediction_workflow.Rmd`.

---

## Results

### Model Performance

| Model | Features | AUROC | Sensitivity | Specificity |
|-------|----------|-------|-------------|-------------|
| Random Forest | Combined | **0.849** ± 0.091 | 0.643 ± 0.134 | 0.828 ± 0.124 |
| Random Forest | Microbiome | 0.772 ± 0.116 | 0.571 ± 0.179 | 0.752 ± 0.156 |
| Random Forest | Clinical | 0.730 ± 0.120 | 0.586 ± 0.161 | 0.724 ± 0.181 |
| Elastic Net | Combined | 0.827 ± 0.088 | 0.671 ± 0.134 | 0.738 ± 0.131 |

*Values reported as mean ± standard deviation across 5 outer CV folds*

### Key Findings

1. **Combination is superior:** Microbiome + clinical data outperforms either alone
2. **Significant improvement:** Combined model AUROC 0.849 exceeds clinical-only baseline (0.730)
3. **Feature importance:** Nutritional variables (folic acid, B vitamins) and specific microbial taxa most predictive
4. **Differential abundance:** 14 genera associated with PTB, including depletion of *Lactobacillus* and enrichment of anaerobic taxa

For complete results with figures and statistical details, see:
- **Rendered report:** `results/integrated_preterm_prediction_workflow.pdf`
- **Figures:** `results/figures/`

---

## Reproducibility

### Exact Version Control

This analysis uses **renv** for dependency management:

```r
# Restore exact package versions
renv::restore()
```

The `renv.lock` file documents:
- R version: 4.4.2
- Package versions with exact commit hashes
- Repository sources (CRAN, Bioconductor, GitHub)

### Computational Environment

**Reproducibility considerations:**
- Fixed random seeds (seed=123) for stochastic operations
- Stratified cross-validation folds for consistent data splits
- Deterministic algorithms where possible

**Limitations:**
- Minor numerical differences may occur across platforms/R versions
- Random Forest uses parallel computation (set `num.threads=1` for exact reproduction)
- Full reproduction requires access to restricted microbiome data

### Session Information

Complete session info (package versions, system details) is included at the end of the rendered analysis report.

---

## Contributing

### Reporting Issues

If you encounter problems or have questions:

1. **Check existing issues:** [GitHub Issues](https://github.com/martinruhle/Mexican-PretermBirth-analysis/issues)
2. **Open new issue:** Provide:
   - Clear description of the problem
   - Steps to reproduce
   - Your system information (R version, OS)
   - Relevant error messages

### Code of Conduct

We are committed to providing a welcoming and respectful environment. Please:
- Be respectful and constructive in discussions
- Focus on scientific merit and reproducibility
- Acknowledge the sensitive nature of health data

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

**Summary:**
- ✅ Free to use, modify, and distribute
- ✅ Can be used for commercial purposes
- ⚠️ Attribution required
- ⚠️ No warranty provided

**Data Usage:** Note that while the code is MIT licensed, access to the full dataset requires separate data use agreement with Instituto Nacional de Perinatología.

---

## Contact

**Principal Investigator:**
- **Martin Ruhle**
- Doctoral Program in Biomedical Sciences
- Instituto Nacional de Perinatología, Mexico City
- Email: [your_institutional_email@inper.edu.mx]

**Repository Issues:**
- GitHub Issues: https://github.com/martinruhle/Mexican-PretermBirth-analysis/issues

**Data Access Inquiries:**
- See [docs/DATA_ACCESS.md](docs/DATA_ACCESS.md)

---

## Acknowledgments

### Institutional Support
- **Instituto Nacional de Perinatología (INPer)**, Mexico City
- Biomedical Sciences Doctoral Program

### Participants
We extend our deepest gratitude to the 43 pregnant women who participated in this study. Their contribution advances our understanding of preterm birth risk factors in Mexican populations.

### Scientific Community
This work builds upon methodological foundations established by:
- The DREAM Preterm Birth Prediction Challenge
- The broader microbiome and machine learning research communities
- Open-source software developers (R, tidymodels, Bioconductor ecosystems)

### Funding
[Add funding sources if applicable]

---

## Additional Resources

### Related Publications
- [Link to manuscript when published]
- [Related work from the research group]

### Learn More
- **Preterm birth:** [March of Dimes](https://www.marchofdimes.org/)
- **Microbiome methods:** [ANCOM-BC2](https://doi.org/10.3389/fmicb.2023.1207829)
- **Machine learning in medicine:** [tidymodels](https://www.tidymodels.org/)

---

**Last Updated:** January 2025

**Repository maintained by:** Martin Ruhle ([GitHub](https://github.com/martinruhle))
