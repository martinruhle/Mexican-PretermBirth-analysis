# Precision of the AUROC and the smallest detectable difference in this cohort

**Last updated:** 2026-11-02

With 43 participants (14 preterm), every AUROC in this repository carries a wide margin. This
document quantifies it, for this cohort and for larger ones:

1. **Precision:** how far an AUROC estimated with these participants can fall from the model's
   true value.
2. **Smallest detectable difference:** the smallest difference in AUROC between two models
   evaluated on the same participants that would be detected with a two-sided test at α = 0.05 and
   80% power. The comparison of interest is the same model with and without microbiome features, as
   in [`CLINICAL_ONLY_ARM.md`](CLINICAL_ONLY_ARM.md).
3. **Participants by gestational-week window:** the counts that a comparison with the DREAM
   Preterm Birth Prediction Challenge (Golob et al. 2024) or a landmark analysis would use.

The pipeline and every model number are unchanged. The inputs are the versioned results of the
permutation test and of the clinical-only arm, and the cohort data for the counts.

**In short:**

- With the metric the pipeline reports (the mean AUROC across the 5 outer folds, with the model
  refitted in each fold), a model whose true AUROC is 0.70 would be measured within about **±0.21**
  (95%).
- **The microbiome features would have to add about 0.28 or more** to a model at 0.70, taking it to
  about 0.98, for the difference to be detected with 80% power. Starting from 0.80, no attainable
  difference reaches that power. For a single fixed model evaluated once on the 43 participants (the
  most favourable case), the figure is about 0.20–0.22.
- The results of [`CLINICAL_ONLY_ARM.md`](CLINICAL_ONLY_ARM.md) and the permutation test are
  therefore compatible both with no contribution of the microbiome features and with a moderate one.
- With this design, differences of about 0.10 require about **300 participants**, and about 0.08
  require about 500.
- 42 of the 43 participants have a sample within the window of the challenge's sub-challenge 1. Sub-challenge 2
  is **not evaluable (≤ 1 case)**.

---

## 1. What is computed and why

The question is about this cohort and this cross-validation, so the calculation follows how the
pipeline measures discrimination:

- **The unit is the participant.** The engine averages the predicted probabilities of each
  participant's samples, and the folds are built by participant, so the 110 samples are not 110
  independent observations. The counts of preterm (n1) and term (n0) participants come from the data
  (section 8).
- **The metric is the one the engine reports:** the mean of the per-fold AUROCs over the 5 outer
  folds. In each fold the variables are selected, the model is fitted on the training participants of
  that fold, and the AUROC is computed on its 7–9 test participants. Its uncertainty has three parts:
  1. which 43 participants were recruited;
  2. averaging five within-fold AUROCs instead of computing one AUROC over all participants;
  3. refitting the model in each fold.

  The first two follow from the variance of the AUROC of Hanley & McNeil (1982). The third cannot be
  estimated without bias from a single cross-validation run (Bengio & Grandvalet 2004), so it is
  measured here from the permutation test, which refits everything in each permutation.
- **The comparison is paired.** Both models are evaluated on the same participants and folds, so
  part of the chance in their AUROCs is shared and cancels in the difference. How much is given by
  the correlation r between the two AUROC estimates (Hanley & McNeil 1983). Here r is estimated from
  the 12 model–set pairs with and without microbiome (section 6) rather than assumed.
- **Each model gets the variance of its own AUROC.** The variance shrinks as the AUROC approaches 1,
  so a model at 0.95 varies less than one at 0.70, and no difference can take a model above 1.
- **The detectable difference Δ** solves (z₁₋α/₂ + z₁₋β)·√V₁ = Δ, with
  V₁ = V(A) + V(A+Δ) − 2r·√(V(A)·V(A+Δ)) and A + Δ ≤ 1. This is the Wald form: the variance of the
  difference is the one under the alternative, as the paired DeLong test (DeLong et al. 1988, as in
  `pROC::roc.test`) estimates it from the data.
- **Two shapes for the variance of the AUROC:** Hanley–McNeil (exponential model) and binormal
  scores (by simulation). Both are reported where they differ.
- **The reference:** the same quantities for a fixed model evaluated once on the 43 participants,
  with a single AUROC over all of them. This is the most favourable case and the one that applies
  to an already-trained external model.

**On the Wald form.** The formula written down before the calculation was the planning form of
Obuchowski (1998), z₁₋α/₂·√V₀ + z₁₋β·√V₁ = Δ, with the critical value taken from the variance under
no difference, V₀ = 2·V(A)·(1 − r). The simulation of section 3 showed that the DeLong test reaches
80% power at smaller differences than that formula predicts. The Wald form was adopted after seeing
the simulation; it agrees with it within the tolerance set beforehand (0.02) at AUROC 0.60 and
0.70. Both forms are in the CSV files.

## 2. The preliminary planning estimate, reproduced

A preliminary estimate used for planning gave about 0.22–0.25 as the smallest detectable
difference with 43 participants. It was computed as 2.8 × SE(A), with the Hanley–McNeil standard
error (SE) of a single AUROC, that is, its standard deviation across repeated cohorts of the same
size:

- **The factor 2.8** is 1.96 + 0.84. The 1.96 sets the threshold that a difference with no true
  effect exceeds 2.5% of the time on each side (α = 0.05, two-sided). The 0.84 places the true
  difference far enough beyond that threshold for 80% of the observed differences to exceed it.
- **The comparison of two models** is implicit. With the same SE for both models and r = 0.5, the SE
  of the difference, √(SE² + SE² − 2·0.5·SE²), equals the SE of one AUROC.

Recomputed with the same method, all seven values match to three decimals:

| Cohort | Participants | Preterm / term | Preterm share | Model AUROC | SE | 95% CI half-width | 2.8 × SE |
|---|---|---|---|---|---|---|---|
| This cohort | 43 | 14 / 29 | 33% | 0.70 | 0.090 | ±0.176 | 0.251 |
| This cohort | 43 | 14 / 29 | 33% | 0.80 | 0.078 | ±0.154 | 0.220 |
| This cohort + DREAM validation Hispanic subgroup | 81 | 29 / 52 | 36% | 0.70 | 0.063 | ±0.123 | 0.176 |
| Hypothetical | 130 | 43 / 87 | 33% | 0.70 | 0.051 | ±0.100 | 0.142 |
| Hypothetical | 200 | 67 / 133 | 34% | 0.70 | 0.041 | ±0.080 | 0.114 |
| Hypothetical | 300 | 100 / 200 | 33% | 0.70 | 0.033 | ±0.065 | 0.093 |
| Hypothetical | 500 | 167 / 333 | 33% | 0.70 | 0.026 | ±0.051 | 0.072 |

The 81-participant row adds to this cohort the 38 Hispanic participants of the challenge's
validation set (15 preterm).

How each of its elements compares with the calculation of section 1:

- **r = 0.5** is consistent with the median r of the 12 pairs (0.48–0.55, section 6). The range
  across pairs is wide, from 0.09 to 0.81.
- **The preterm share of the hypothetical rows (one third)** matters little. With this cohort's
  share (14/43) the values change by 0.001–0.002. The 81-participant row changes more, 0.176 to
  0.182, because its 29/52 has a higher preterm share than the cohort.
- **The same SE for both models** gives a slightly larger value from a model at 0.70, because the
  model with the added contribution would vary less. From 0.80, 0.80 + 0.220 = 1.02 lies above the
  maximum AUROC of 1, so 0.220 is not an attainable difference.
- **The question it answers** is the fixed-model one (a single AUROC over all participants). For
  that question it is conservative: the simulated test detects 0.202 from 0.70 (section 3). For the
  pipeline's metric the detectable difference is larger (section 4).

## 3. Check by simulation: fixed model, AUROC over all participants

Paired binormal scores for 14 preterm and 29 term participants were simulated, with the correlation
between the two AUROC estimates set to about 0.5. The paired DeLong test was implemented in the
script and checked against `pROC`. Power was estimated with 2,000 replicates per difference. The
table gives the difference at which each method reaches 80% power:

| Model AUROC | Simulated DeLong test | Wald form, Hanley–McNeil | Wald form, binormal | Obuchowski form, Hanley–McNeil | 2.8 × SE (section 2) |
|---|---|---|---|---|---|
| 0.60 | 0.237 | 0.240 | 0.242 | 0.257 | 0.266 |
| 0.70 | 0.202 | 0.219 | 0.214 | 0.241 | 0.251 |
| 0.80 | 0.165 | not attainable | 0.192 | not attainable | 0.220 (to 1.02) |

- **False-positive rate of the simulated DeLong test** with no true difference: 0.061, 0.058 and
  0.048 at AUROC 0.60, 0.70 and 0.80.
- **Monte Carlo uncertainty:** about ±0.01–0.02 for the simulated values.
- **Agreement with the simulation:** the Wald form is within 0.02 at 0.60 and 0.70. The Obuchowski
  form is 0.039 above it at 0.70.
- **At 0.80:** with the Hanley–McNeil variance, no attainable difference reaches 80% power. With
  binormal scores the simulated test reaches it at 0.165.

## 4. The pipeline's metric: five folds and refitting

The outer folds hold 3 preterm and 6 term participants (four folds) and 2 and 5 (one fold). The
within-fold AUROCs compare 82 of the 406 preterm–term pairs of the cohort.

- **Averaging five within-fold AUROCs:** the SE is 1.058 times that of a single AUROC over all
  participants (1.064 with binormal scores), at AUROC 0.70.
- **Refitting in each fold.** This is measured with the permutation test (199 permutations, all 18
  combinations refitted in each).
  - **With the pipeline:** the SD of the permuted AUROCs is 0.125 (median over the 18 combinations).
  - **Without refitting:** with random scores on the same folds, the same permutations and the same
    rule for skipped folds, the SD is 0.109.
  - **The ratio, f = 1.153** (bootstrap 95% CI 1.10–1.20; 1.06 to 1.23 across the 18 combinations),
    is the variance that refitting adds.
- **SE of the pipeline's AUROC** (Hanley–McNeil × 1.058 × f): 0.115, 0.109 and 0.096 at 0.60, 0.70
  and 0.80, that is, 95% half-widths of **±0.226, ±0.214 and ±0.188**.
- **Smallest detectable difference** with the pipeline's metric (Wald form, r = 0.5):

| Model AUROC | Hanley–McNeil | Binormal |
|---|---|---|
| 0.60 | 0.285 | 0.304 |
| 0.70 | not attainable | 0.275 |
| 0.80 | not attainable | not attainable |

"Not attainable" means that no difference with A + Δ ≤ 1 reaches 80% power.

**Assumption:** f was measured with shuffled labels, where the true AUROC is 0.5. It is assumed to
hold for models at 0.60–0.80 and for the larger cohorts of section 7.

## 5. How optimistic are the fixed-model figures?

Three reasons were considered, each checked separately:

- **The spread between folds (SD about 0.2) is not an additional source.** The between-fold
  variance of the 18 combinations is 1.045 times the one expected from sampling the 7–9 test
  participants of each fold alone (95% range under that sampling alone: 0.72–1.32). That spread is
  what folds of this size produce, and the cohort-level SE already includes it: 0.2 / √5 ≈ 0.09.
- **Refitting is an additional source.** f = 1.15 (section 4). The five models are fitted on
  overlapping participants (any two folds share three quarters of their training participants), so
  their deviations from the true AUROC are correlated across folds and the mean over folds varies
  more than the per-fold spread suggests. For the pipeline's metric, the fixed-model figures are therefore optimistic.
- **The choice among combinations is an additional source.**
  - **Under the permutation null,** one combination has a mean AUROC of 0.493 and the best of the 18
    has 0.662 ± 0.083: the best combination is expected to be 0.17 higher with no signal at all.
  - **If the 12 microbiome pairs were tested together** with a Bonferroni correction, the detectable
    difference would be 1.32 times larger.

For a fixed model and a single AUROC over all participants, the preliminary 0.251 is conservative:
the simulated test detects 0.202 (section 3).

## 6. Correlation between the AUROCs of the paired models

r was estimated for the 12 pairs of [`CLINICAL_ONLY_ARM.md`](CLINICAL_ONLY_ARM.md) (each clinical-only
model against the same model with ANCOM-BC2 taxa or with the full microbiome):

| Estimate | Median r | Range across the 12 pairs |
|---|---|---|
| DeLong covariance, AUROC over all participants | 0.546 | 0.108–0.803 |
| DeLong covariance, mean of the folds | 0.480 | 0.087–0.808 |
| Permutation null (same permutation for both models) | 0.443 | 0.122–0.803 |

The tables use r = 0.5. The CSV files also give r = 0.3, 0.7 and the three medians. A lower r
means the two models share less of the chance, so a larger difference is needed.

## 7. Precision and detectable difference by cohort size

Model AUROC 0.70, r = 0.5, Wald form, Hanley–McNeil variance:

| Cohort | Participants | Preterm (share) | 95% CI half-width: fixed model / pipeline | Detectable difference: fixed model / pipeline |
|---|---|---|---|---|
| This cohort | 43 | 14 (33%) | ±0.176 / ±0.214 | 0.219 / not attainable (binormal 0.275) |
| DREAM sub-challenge 1 window (section 8) | 42 | 13 (31%) | ±0.181 / ±0.218 | 0.225 / not attainable |
| Landmark week 32 (section 8) | 40 | 12 (30%) | ±0.188 / ±0.227 | 0.233 / not attainable |
| This cohort + DREAM validation Hispanic subgroup | 81 | 29 (36%) | ±0.123 / ±0.146 | 0.159 / 0.185 |
| Hypothetical | 81 | 27 (33%) | ±0.126 / ±0.149 | 0.163 / 0.189 |
| Hypothetical | 130 | 43 (33%) | ±0.100 / ±0.117 | 0.131 / 0.152 |
| Hypothetical | 200 | 67 (34%) | ±0.080 / ±0.093 | 0.107 / 0.123 |
| Hypothetical | 300 | 100 (33%) | ±0.065 / ±0.076 | 0.089 / 0.102 |
| Hypothetical | 500 | 167 (33%) | ±0.051 / ±0.058 | 0.070 / 0.080 |

- **The pipeline column for hypothetical cohorts** assumes the same design: 5 folds of N/5
  participants and the same f.
- **The fixed-model column applies to an already-trained external model.** For example, a challenge
  model with AUROC 0.60–0.70 applied to the 42 participants of the sub-challenge 1 window would be
  measured within ±0.18–0.19.
- **The CSV files** also give AUROC 0.60 and 0.80, the hypothetical rows with this cohort's preterm
  share (they differ by at most 0.002), and the Obuchowski form.

## 8. Participants by gestational-week window

Gestational ages are recorded as weeks+days: 32+6 is 32 weeks and 6 days. Participants are counted
when they have at least one sample in the window.

| Window | Outcome | Samples counted | Participants (preterm / term) | Eligible without a sample in the window |
|---|---|---|---|---|
| All participants | preterm (< 37+0) | all | 43 (14 / 29) | — |
| DREAM sub-challenge 1 | preterm (< 37+0) | before 32+0 | 41 (13 / 28) | 2 (1 preterm) |
| **DREAM sub-challenge 1 (chosen)** | preterm (< 37+0) | **through 32+6** | **42 (13 / 29)** | 1 (preterm) |
| DREAM sub-challenge 2 | delivery < 32+0 | through 28+6 | **not evaluable (≤ 1 case)** | |
| Landmark week 24, 42 who reached 24+0 | preterm (< 37+0) | before 24+0 | 31 (8 / 23) | 11 (5 preterm) |
| Landmark week 28, 42 who reached 28+0 | preterm (< 37+0) | before 28+0 | 38 (11 / 27) | 4 (2 preterm) |
| **Landmark week 32, 42 who reached 32+0 (chosen)** | preterm (< 37+0) | **before 32+0** | **40 (12 / 28)** | 2 (1 preterm) |

- **Timing of the preterm deliveries:** ≤ 1 before 32+0 and 4 before 34+0. The preterm births of
  this cohort are mostly late preterm.
- **The windows were fixed by definition, not by the counts.**
  - **Sub-challenge 1 and 2 windows:** the challenge's task description (Synapse syn26133770) admits specimens "collected
    no later than 32 weeks" (sub-challenge 1) and "no later than 28 weeks" (sub-challenge 2), and
    gives the collection week in whole weeks, so a sample taken at 32+6 is in week 32. How the
    challenge rounded the collection week is not documented, so this is an assumption. The other
    reading, before 32+0, removes one term participant (42 → 41).
  - **Landmark windows:** they count samples taken before week X+0 among the participants still
    pregnant at X+0, as in [`SAMPLE_COUNT_REFERENCE.md`](SAMPLE_COUNT_REFERENCE.md). Every sample
    then precedes the landmark week.
- **Both readings are in the CSV** for sub-challenge 1 and the landmarks.

## 9. How to read these numbers

- **They are planning values.** They depend on approximations to the variance of the AUROC
  (Hanley–McNeil or binormal) and on a simulation with a Monte Carlo uncertainty of about ±0.01–0.02.
  The conclusions do not depend on the third decimal: in this cohort only very large differences are
  detectable, and differences of about 0.10 need a few hundred participants.
- **The refitting factor f** was measured under shuffled labels, on one partition of the 43
  participants (seed 123), and is carried over to other AUROCs and cohort sizes.
- **They do not say whether the microbiome features contribute.** They say how large a contribution
  would have to be for this cohort to show it.
- **The selection among combinations** (section 5) is not included in the tables. A model chosen as
  the best of several needs the selection-adjusted test of
  [`results/permutation_test/`](../results/permutation_test/).

## 10. Files and reproducing these numbers

[`results/auroc_precision/`](../results/auroc_precision/) holds the numbers of this document; its
README describes each file. They are written by
[`scripts/auroc_precision.R`](../scripts/auroc_precision.R):

```bash
PTB_RESULTS_DIR=results/auroc_precision Rscript scripts/auroc_precision.R
```

- **What the script reads.** The formulas, the simulation and the permutation-based factors use only
  the versioned results of [`results/permutation_test/`](../results/permutation_test/) and
  [`results/clinical_only_arm/`](../results/clinical_only_arm/).
- **What needs the cohort data.** The window counts and the correlations of section 6 need the cohort
  data in `data/raw/` and the saved folds of the clinical-only arm run (`clinonly_2026-10-16`, see
  [`CLINICAL_ONLY_ARM.md`](CLINICAL_ONLY_ARM.md) §6). Without them the script skips that part and
  the tables use the fixed values of r.
- **Checks.** The script stops unless the rebuilt outer folds are identical to those of the 18
  combinations of that run and its per-fold AUROC reproduces the engine's in all 18.
- **Runtime.** It took 1.0 and 1.8 minutes in two runs (Windows 11, R 4.4.2). The seeds are fixed,
  and both runs gave the same numbers.

## References

- DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more correlated receiver
  operating characteristic curves: a nonparametric approach. *Biometrics* 1988;44:837–845.
- Bengio Y, Grandvalet Y. No unbiased estimator of the variance of K-fold cross-validation.
  *J Mach Learn Res* 2004;5:1089–1105.
- Golob JL, et al. Microbiome preterm birth DREAM challenge: crowdsourcing machine learning
  approaches to advance preterm birth research. *Cell Rep Med* 2024.
- Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic
  (ROC) curve. *Radiology* 1982;143:29–36.
- Hanley JA, McNeil BJ. A method of comparing the areas under receiver operating characteristic
  curves derived from the same cases. *Radiology* 1983;148:839–843.
- Obuchowski NA. Sample size calculations in studies of test accuracy. *Stat Methods Med Res*
  1998;7:371–392.
