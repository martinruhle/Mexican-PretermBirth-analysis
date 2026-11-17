# Clinical variable availability relative to sample collection

**Last updated:** 2026-10-07

A prediction made at a prenatal visit can only use information that exists at that visit. This
document classifies every clinical variable in
[`config/data_dictionary.csv`](../config/data_dictionary.csv) by **when its value is known relative
to the visit at which the vaginal sample was taken**, records the evidence for each classification,
and lists the variables used by the models before this classification that are **not** known at the
visit, with the number of cross-validation folds in which each one entered a model.

The pipeline now reads the classification: only the classes listed under
`clinical_availability$allowed` in [`config/config.yml`](../config/config.yml) can enter a model.
The re-run of the 12 combinations with visit-time variables only, with the before/after table, is
reported in [`CLINICAL_AVAILABILITY_RESULTS.md`](CLINICAL_AVAILABILITY_RESULTS.md). The re-run uses
two variable sets:

- **Main analysis:** `at_visit` variables only.
- **Sensitivity analysis:** `at_visit` plus `after_visit` variables, to measure how much of the
  current performance depends on information that may not exist at the visit. In the clinical pool
  this adds back the undated obstetric complications (`comppreeclam`, `compvaginf`,
  `diabetes_gest`, `oligohidramnios`, `rciu`) and `sex_baby`.

`outcome_defined` variables are not used in either set.

---

## 1. Classes

The dictionary has three new columns: `availability`, `availability_source` (the evidence, one
sentence per variable) and `derived_from` (inputs of derived variables).

| `availability` | Meaning |
|---|---|
| `at_visit` | Known before or at the visit: measured at the visit, defined as prior to the pregnancy, or recorded at enrolment |
| `after_visit` | Becomes known after the visit: pregnancy events or diagnoses without a date in the cohort data, and characteristics of the birth |
| `outcome_defined` | Defined by the delivery itself: gestational age at delivery, the outcome classification, or the 37-week cut-off that defines the outcome |

Only `at_visit` variables are eligible as predictors.

## 2. Rules

1. **Measured at the visit → `at_visit`.** Its value changes between visits of the same participant
   (weight, 24-hour dietary recall, fetal biometry, hemoglobin, anemia). For weight, the dietary
   recall and fetal weight, the cohort's variable dictionary also says the value is taken "at the
   visit" or "on the consultation date".
2. **Defined before the pregnancy or recorded at enrolment → `at_visit`.** Pre-pregnancy weight and
   BMI, height, and maternal age (age at enrolment is an inclusion criterion of the study).
3. **Defined by the delivery → `outcome_defined`.**
4. **Anything the cohort data cannot place at or before the visit → `after_visit`.** The cohort's
   original clinical database has one row per visit and no diagnosis dates. A variable that takes the
   same value at every visit of a participant is a single value per pregnancy, with no record of when
   it was entered; unless its definition places it before the visit, it is treated as known after the
   visit.
5. **Derived variables inherit the latest class among their inputs** (order: `at_visit` <
   `after_visit` < `outcome_defined`). Inputs are listed in `derived_from`. The variables the pipeline
   creates itself (chunk `feature_engineering` of the analysis notebook) are now dictionary rows with
   `role = clinical_derived`; they are not columns of the input matrix.

Three classifications rest on an assumption rather than on a date in the data. All three are
recorded in `availability_source`, and all can be checked against the cohort's case report form:

- **Supplement use** (`dietsuppl3mon`, `vitaminsup`, `supintakefreq`) is classed `at_visit`: the
  study team reports that it is asked at the prenatal visits, and its constant value across a
  participant's visits is consistent with supplementation maintained throughout pregnancy. Unlike
  the obstetric complications, these variables are not part of the per-pregnancy block described in
  §3.
- **Sociodemographic variables** (`nivel_academico`, `maritalstat`, `workouthome`) are classed
  `at_visit` as enrolment characteristics. The data do not record when they were collected.
  `workouthome` is the weakest of the three, because employment can change during pregnancy.
- **First-trimester bleeding** (`comp1tribleed`) is classed `at_visit` because, by definition, it
  occurs before week 14, and 102 of the 110 samples were taken at week 14 or later. The 8 first-trimester
  samples carry a value that could still change after the visit. In the cohort's dictionary, the
  entry for this variable carries the note *Desenlace* (outcome); it is the first entry of the
  obstetric complications (`comp1tribleed` to `bajo_peso_nac`).

## 3. Evidence from the data

All figures below are computed by
[`scripts/check_variable_availability.R`](../scripts/check_variable_availability.R) from the cohort
matrix (110 samples, 43 participants: 29 term, 14 preterm; 31 participants have two or more samples).

- **The obstetric complications are stored as one value per pregnancy.** The question for these
  variables is not whether they are assessed during prenatal care, but whether the value stored in
  each sample's row is the status *at that visit*. In the cohort's dictionary, the obstetric
  complications (`comp1tribleed` to `bajo_peso_nac`) form one block, and that block also holds
  events of the delivery itself: `rpm` (rupture of membranes after week 37) equals 1 at visits taken
  at weeks 14+4 and 19+5, and `rpm_preterm` equals 1 at a 19+5-week visit of a pregnancy that ended
  at 32+5 weeks. The block therefore stores the pregnancy's final status, copied onto every visit;
  a diagnosis made late in pregnancy appears in the rows of earlier samples. Consistent with this,
  every variable of the block takes the same value at every visit of every participant with two or
  more samples, and `compvaginf` (active infection) equals 1 at all three visits of the two
  participants with the condition (weeks 11 to 27, and 20 to 38). For `oligohidramnios`,
  `comppreeclam` and `rciu`, the participants with a value of 1 have a single visit each, so the
  constancy itself says nothing about them; their class rests on the block's structure.
- **Other variables with a single value per pregnancy.** The supplement variables, the birth
  variables (`birthweightgr`, `peso_nacimiento`, `sex_baby`) and the sociodemographic variables are
  also constant across a participant's visits. The cohort's original database (121 visits, 46
  participants) shows the same pattern for all of them.
- **Measured at each visit.** Weight, the dietary recall variables, fetal biometry and hemoglobin
  change between visits in all or nearly all participants with two or more samples (for example,
  weight in 30 of 30, hemoglobin in 26 of 27); `anemia_visita` changes in 7 of 27.
- **Preeclampsia was diagnosed after the sampled visit.** `comppreeclam` equals 1 in one participant
  (preterm, delivery at 34+6 weeks), whose only visit is visit 1, the enrolment visit, at 17+2
  weeks. The study's exclusion criteria removed participants with complications present at
  enrolment (article §2.1.2), so the diagnosis came after that visit. The same participant is the
  only one with `rciu` = 1.
- **The rupture-of-membranes flags follow the outcome.** `rpm` (rupture after week 37) equals 1 in 5
  of 29 term participants and in none of the 14 preterm participants: a value of 1 implies the
  pregnancy reached 37 weeks. `rpm_preterm` (rupture before week 37) equals 1 in 3 of 14 preterm
  participants and in no term participant; it uses the same cut-off as the outcome, and preterm
  rupture of membranes is one of the clinical presentations of preterm birth.
- **Derived columns.** `hemoglobin_alti_adj` = `hemoglobin_g_dl` − 0.95 in every sample;
  `bajo_peso_nac` = `birthweightgr` < 2500 in every sample with data; `lag_peso_kg`,
  `lag_hemoglobin` and `rolling_avg_peso` use the current and previous visit only (none uses a later
  visit); `imc_visita` is identical to `imc_pregestacional` in all 110 samples.
- **Gestational ages use weeks.days notation.** `sdg_visita` and `sdg_parto` are written as
  weeks.days (24.4 = 24 weeks and 4 days), and they match the day counts `dg_visita` and `dg_parto`
  in 110 of 110 samples.

## 4. The clinical variables available to the models

Approach 1 uses `sdg_visita` and `edad_cronologicamujer`, both `at_visit`. The Approach 2
candidates (17, from `config/config.yml`) are all part of the Approach 3 pool, which is every
clinical variable the notebook prepares (chunks `extract_clinical` and `feature_engineering`): 37
variables before the completeness filter.

| Class | Variables in the Approach 3 pool |
|---|---|
| `at_visit` (26) | `sdg_visita`, `edad_cronologicamujer`, `nivel_academico`, `maritalstat`, `workouthome`, `peso_pregestacional_kg`, `talla_mujer_cm`, `imc_pregestacional`, `imc_pregest_categ`, `comp1tribleed`, `hemoglobin_g_dl`, `hemoglobin_alti_adj`, `anemia_visita`, `dietsuppl3mon`, `vitaminsup`, `folic_ac_correg`, `peso_kg`, `imc_visita`, and the derived `age_risk_category`, `extreme_age`, `underweight`, `obese`, `extreme_bmi`, `low_education`, `unmarried`, `ses_risk_score` |
| `after_visit` (6) | `comppreeclam`, `compvaginf`, `diabetes_gest`, `oligohidramnios`, `rciu`, `sex_baby` |
| `outcome_defined` (5) | `rpm`, `rpm_preterm`, and the derived `complication_count`, `any_complication`, `multiple_complications` (their inputs include `rpm` and `rpm_preterm`) |

## 5. Variables that are not known at the visit, and how often the earlier models used them

Counts come from the baseline run before this restriction (engine at commit `ac9af7a`). The Approach 2 variable list
is fixed and is the same in all 5 outer folds; Approach 3 selects its variables in each outer fold,
and the selection is identical in its four combinations (2 models × 2 microbiome inputs).

| Variable | Class | Why | Approach 2 (final list of 10) | Approach 3 (folds of 5) |
|---|---|---|---|---|
| `rpm` | `outcome_defined` | Rupture after week 37 implies a term delivery | **Used, all 5 folds** | 0 |
| `comppreeclam` | `after_visit` | Diagnosed after the enrolment visit, the only sampled visit | **Used, all 5 folds** | 0 |
| `oligohidramnios` | `after_visit` | Per-pregnancy status without a diagnosis date | **Used, all 5 folds** | **3** |
| `compvaginf` | `after_visit` | Per-pregnancy status without a date | Candidate, not selected | **3** |
| `complication_count` | `outcome_defined` | Inherits from `rpm` and `rpm_preterm` | — | **1** |
| `any_complication` | `outcome_defined` | Inherits from `complication_count` | — | **1** |
| `diabetes_gest` | `after_visit` | Per-pregnancy status without a date (0 in every participant) | Candidate, not selected | 0 |
| `rpm_preterm` | `outcome_defined` | Same 37-week cut-off as the outcome | — | 0 |
| `rciu` | `after_visit` | Per-pregnancy status; same participant as `comppreeclam` | — | 0 |
| `sex_baby` | `after_visit` | Known at birth or from an undated ultrasound | — | 0 |
| `multiple_complications` | `outcome_defined` | Inherits from `complication_count` | — | 0 |

In summary, **8 of the 12 combinations used at least one variable that is not known at the visit,
in every outer fold**: the four Approach 2 combinations (three such variables in each fold) and the
four Approach 3 combinations (one to three per fold: 3, 1, 1, 1 and 2 in folds 1 to 5). The four
Approach 1 combinations use none.

## 6. The delivery at 20.0 weeks

One of the 14 participants counted as preterm has a pregnancy that ended at 20+0 weeks (140 days);
the only sample from this participant was taken at 19+4 weeks. The analysis outcome `preterm`
(gestational age at delivery < 37 weeks) codes this pregnancy as 1, while the cohort's own outcome
variable `desenlace_parto` codes it as `Abortion`. The published article kept this participant in
the preterm group, citing a WHO threshold of 20 weeks, and reported a sensitivity analysis without
this participant (42 participants, 13 preterm).

The reference definitions place the boundary at 22 weeks:

- The Mexican standard for pregnancy, delivery and puerperium care,
  [NOM-007-SSA2-2016](https://www.gob.mx/cms/uploads/attachment/file/512098/NOM-007-SSA2-2016.pdf),
  defines delivery (*parto*) as the expulsion of a fetus of 22 weeks or more (§3.31), abortion as the
  expulsion of an embryo or fetus under 500 g, "a weight reached at approximately 22 completed weeks"
  (§3.1), and the perinatal period as starting at 22 weeks (§3.35).
- ICD-11 starts the perinatal period at 154 days of gestation, 22+0 completed weeks
  ([Blencowe et al. 2025, *Int J Gynecol Obstet* 168:1–9](https://pmc.ncbi.nlm.nih.gov/articles/PMC11649847/)),
  and WHO defines preterm as babies born alive before 37 completed weeks
  ([WHO fact sheet](https://www.who.int/news-room/fact-sheets/detail/preterm-birth)).

Under the cohort's own coding and the 22-week definitions, this pregnancy would be classed as an
abortion; under the 20-week threshold cited in the article, as a preterm birth. **The participant is
retained as preterm, as in the published article.** With 14 preterm participants, removing one is a
substantial loss of information, and the aim of the study is to predict and characterize preterm
birth, to which a pregnancy ending at 20 weeks is relevant. The cohort's coding is recorded here so
that readers can apply the stricter definition. The article reported a sensitivity analysis
without this participant; it will be repeated on the updated analysis (visit-time variables only),
so that the effect of this choice is reported with the current pipeline.

## 7. Reproducing these numbers

```bash
PTB_RUN_TAG=thresholdfix_2026-09-30 Rscript scripts/check_variable_availability.R
```

The script needs the cohort data in `data/raw/` and the `<tag>_cv_folds.rds` that
`scripts/run_baseline.R` writes for a run. The dictionary columns are checked by
[`tests/testthat/test-data-dictionary.R`](../tests/testthat/test-data-dictionary.R): every row has a
class and a source, the outcome and rupture-of-membranes flags are `outcome_defined`, every variable
named in the approach definitions has a row, and every derived variable carries the latest class of
its inputs.
