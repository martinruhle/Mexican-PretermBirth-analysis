# Data Documentation

## Overview

This directory contains publicly available metadata from the preterm birth prediction study conducted at Instituto Nacional de Perinatología, Mexico City. The complete dataset includes restricted microbiome abundance data that is not publicly available due to ethical and privacy considerations.

---

## Available Data

### Metadata Files (Public)

All metadata files are located in `data/metadata/` and are publicly available:

#### 1. `diccionario_variables_completo.csv`

**Description:** Complete variable dictionary providing definitions and descriptions of all variables used in the study.

**Structure:**
- **Rows:** 67 variables
- **Columns:** 2
  - `Diccionario de variables`: Variable name (Spanish/technical)
  - Second column: Variable description/definition

**Purpose:** Reference guide for understanding variable names and their meanings in the longitudinal metadata file.

**Encoding:** CP1252 (Windows Latin-1)

---

#### 2. `metadata_eugenia_long.csv`

**Description:** Longitudinal clinical, nutritional, and outcome data for all study participants and visits.

**Structure:**
- **Rows:** 110 observations (samples/visits)
- **Participants:** 43 unique women
- **Columns:** 70 variables

**Key Variable Categories:**

**Identifiers:**
- `index`: Sample/visit index
- `id`: Participant unique identifier  
- `visita`: Visit designation

**Gestational Age:**
- `dg_visita`: Gestational day at visit
- `sdg_visita`: Gestational weeks + days at visit
- `dg_parto`: Gestational day at delivery
- `sdg_parto`: Gestational weeks + days at delivery

**Maternal Characteristics:**
- `edad_cronologicamujer`: Maternal age (years)
- `peso_pregestacional_kg`: Pre-pregnancy weight (kg)
- `talla_mujer_cm`: Height (cm)
- `imc_pregestacional`: Pre-pregnancy BMI
- `imc_pregest_categ`: Pre-pregnancy BMI category
- `peso_kg`: Weight at visit (kg)

**Nutritional Intake (24-hour recall):**
- `calories`: Total caloric intake
- `calories_fat`: Calories from fat
- `carbohydrates`: Carbohydrate intake (g)
- `protein`: Protein intake (g)
- `fat`: Total fat intake (g)
- `vitamin_b1`, `vitamin_b2`, `vitamin_b6`, `vitamin_b12`: B vitamin intake
- `folic_ac_correg`: Corrected folic acid intake
- `folate`: Folate intake
- `folate_dfe`: Dietary folate equivalents
- `folate_food`: Folate from food sources
- `choline`: Choline intake
- Amino acids: `cystine`, `glycine`, `methionine`, `serine`

**Dietary Supplements:**
- `dietsuppl3mon`: Dietary supplement use in past 3 months
- `vitaminsup`: Vitamin supplement use
- `supintakefreq`: Supplement intake frequency

**Complications (Binary: 0/1):**
- `comp1tribleed`: First trimester bleeding
- `compvaginf`: Vaginal infection
- `compsexualinf`: Sexual infection  
- `comppreeclam`: Preeclampsia
- `rpm_preterm`: Preterm premature rupture of membranes
- `rpm`: Premature rupture of membranes
- `diabetes_gest`: Gestational diabetes
- `obito`: Fetal death
- `oligohidramnios`: Oligohydramnios
- `rciu`: Intrauterine growth restriction

**Fetal Outcomes:**
- `bajo_peso_nac`: Low birth weight
- `sex_baby`: Baby's sex
- `birthweightgr`: Birth weight (grams)
- `peso_nacimiento`: Birth weight category

**Fetal Biometry:**
- `pfetal`: Fetal weight estimate
- `fcf`: Fetal heart rate
- `ccef`: Cephalic circumference
- `dbip`: Biparietal diameter
- `cabd`: Abdominal circumference
- `longfe`: Femur length

**Sociodemographic:**
- `nivel_academico`: Educational level
- `maritalstat`: Marital status
- `workouthome`: Work outside home

**Laboratory Values:**
- `hemoglobin_g_dl`: Hemoglobin (g/dL)
- `hemoglobin_alti_adj`: Altitude-adjusted hemoglobin
- `anemia_visita`: Anemia status at visit

**Outcomes:**
- `desenlace_parto`: Delivery outcome
- `preterm`: Preterm birth (<37 weeks) - **PRIMARY OUTCOME**
- `early_preterm`: Early preterm birth (<34 weeks)

**Derived Variables:**
- `imc_visita`: BMI at visit
- `lag_peso_kg`: Lagged weight (previous visit)
- `lag_hemoglobin`: Lagged hemoglobin
- `lag_imc_visita`: Lagged BMI
- `rolling_avg_peso`: Rolling average weight
- `rolling_avg_imc`: Rolling average BMI
- `time_since_first`: Time since first visit

**Encoding:** UTF-8

**Missing Data:**
- Variable missingness varies by visit timing and data availability
- Nutritional data collected at select visits
- Fetal biometry timing dependent on gestational age
- Not all complications occur in all participants

---

#### 3. `participant_data_clean.csv`

**Description:** Participant-level summary data (one row per participant).

**Structure:**
- **Rows:** 43 participants
- **Columns:** 13 variables

**Variables:**
- `id`: Participant unique identifier
- `edad`: Age (years)
- `peso_pregest`: Pre-pregnancy weight (kg)
- `talla`: Height (cm)
- `imc_cat`: BMI category
- `educ_level`: Education level
- `marital_status`: Marital status
- `workouthome`: Work outside home
- `preterm_status`: Preterm birth outcome (0/1)
- `early_preterm_status`: Early preterm outcome (0/1)
- `anemia`: Anemia status (0/1)
- `diabetes`: Diabetes status (0/1)
- `n_visits`: Number of study visits

**Purpose:** Participant-level characteristics and outcomes for descriptive statistics and baseline comparisons.

**Encoding:** UTF-8

---

### Example Data

#### `data/example/abundance_example.csv`

**Description:** Small synthetic example demonstrating the structure of microbiome abundance data.

**Structure:**
- **Rows:** 10 simulated samples
- **Columns:** ~50 (sample ID + genus-level CLR-transformed abundances)

**Important Notes:**
- **NOT real data** - simulated for demonstration only
- Shows expected data format and column structure
- CLR-transformed values (centered log-ratio)
- Column names format: `Genus_GenusName` or `CLR_GenusName`

**Purpose:** Enable users to understand data format and test code without access to restricted data.

---

## Restricted Data (Not Included)

### Microbiome Abundance Tables

**What is restricted:**
- Raw 16S rRNA gene sequencing counts
- Genus-level abundance tables
- Species-level abundance profiles (if available)
- Sample-level taxonomic profiles
- Alpha diversity metrics
- Beta diversity distance matrices

**Why restricted:**
1. **Ethical approval:** IRB restrictions on public data sharing
2. **Privacy:** Microbiome profiles can be re-identifying
3. **Data sharing agreements:** Institutional agreements with INPer
4. **Participant consent:** Limited to specified research uses

**What we can share:**
- Aggregated statistics (differential abundance results)
- Model coefficients for microbial features
- Visualizations of group-level patterns
- Example data structure (synthetic)

**Access process:** See [docs/DATA_ACCESS.md](../docs/DATA_ACCESS.md)

---

## Data Dictionary

### Cross-Reference with Variable Dictionary

To understand any variable in `metadata_eugenia_long.csv`:

1. Find the variable name in the file
2. Look up the corresponding entry in `diccionario_variables_completo.csv`
3. Read the description (may be in Spanish)

**Example:**
```r
# Load dictionary
dictionary <- read.csv("data/metadata/diccionario_variables_completo.csv", 
                       fileEncoding = "CP1252")

# Look up a variable
dictionary[dictionary$`Diccionario.de.variables` == "edad_cronologicamujer", ]
```

---

## Data Quality

### Completeness

**Longitudinal structure:**
- Mean visits per participant: 2.6 (range: 1-5)
- Not all variables collected at all visits
- Visit timing varies by participant and clinical indication

**Missing data patterns:**
- **Nutritional data:** Collected at 1-2 visits per participant
- **Fetal biometry:** Timing dependent on gestational age
- **Complications:** Only recorded when they occur
- **Laboratory values:** May not be available at all visits

### Quality Control

**Data cleaning performed:**
- Outlier detection for anthropometric measurements
- Biological plausibility checks for gestational ages
- Consistency checks for delivery outcomes
- Duplicate record resolution

**Variables derived:**
- BMI calculated from weight and height
- Gestational age standardized to weeks + days
- Lagged and rolling average variables for longitudinal analysis
- Altitude-adjusted hemoglobin for Mexico City elevation

---

## Data Usage Guidelines

### Ethical Considerations

1. **Respect participant privacy:** Do not attempt to re-identify individuals
2. **Use for stated purposes:** Scientific research on preterm birth
3. **Acknowledge contributors:** Cite the manuscript and acknowledge INPer
4. **Share responsibly:** Do not redistribute raw data without permission

### Statistical Considerations

1. **Longitudinal structure:** Account for repeated measures within participants
2. **Missing data:** Use appropriate imputation or complete-case analysis
3. **Sample size:** Limited power for rare outcomes (n=43 participants)
4. **Compositional data:** Microbiome abundances require appropriate transformations (CLR)
5. **Class imbalance:** 14/43 preterm births (~33%) - stratification important

### Reproducibility

**Loading data in R:**
```r
# Variable dictionary (CP1252 encoding)
dictionary <- read.csv("data/metadata/diccionario_variables_completo.csv",
                       fileEncoding = "CP1252")

# Longitudinal metadata (UTF-8 encoding)
metadata_long <- read.csv("data/metadata/metadata_eugenia_long.csv",
                          fileEncoding = "UTF-8")

# Participant-level data (UTF-8 encoding)
participant_data <- read.csv("data/metadata/participant_data_clean.csv",
                             fileEncoding = "UTF-8")

# Example abundance structure
example_abundance <- read.csv("data/example/abundance_example.csv")
```

---

## Data Access

### Public Data
The metadata files in this repository are publicly available under the repository license (MIT). You may use them freely with appropriate attribution.

### Restricted Data
To request access to complete microbiome abundance data:

1. Review requirements in [docs/DATA_ACCESS.md](../docs/DATA_ACCESS.md)
2. Prepare required documents:
   - Research protocol
   - Ethics approval
   - Data use agreement
   - Scientific justification
3. Submit request to Principal Investigator

**Expected timeline:** 2-4 weeks for review and approval

---

## Updates and Versioning

**Current version:** v1.0 (January 2025)
- Initial public release with manuscript submission
- Metadata reflects final analysis dataset used in manuscript

**Future updates:**
- Additional follow-up data (if collected)
- Corrections or clarifications (if needed)
- Will be documented in this file with version numbers

---

## Contact

**Data-related questions:**
- Principal Investigator: Martin Ruhle
- Institution: Instituto Nacional de Perinatología
- Email: [your_institutional_email@inper.edu.mx]

**Technical issues with files:**
- GitHub Issues: https://github.com/martinruhle/Mexican-PretermBirth-analysis/issues

---

## Additional Resources

### Related Documentation
- [Installation Guide](../docs/INSTALL.md)
- [Data Access Request Process](../docs/DATA_ACCESS.md)
- [Main Repository README](../README.md)

### Data Standards
This study follows:
- TRIPOD+AI guidelines for prediction model reporting
- STROBE guidelines for observational studies
- MIxS standards for microbiome metadata

---

**Last Updated:** January 2025
