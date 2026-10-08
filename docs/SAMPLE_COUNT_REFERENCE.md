# Number of samples per participant as a reference for the AUROC

**Last updated:** 2026-10-08

A participant who delivers preterm stops contributing samples earlier, because the pregnancy ends
earlier. In this cohort the 14 preterm participants have 1.71 samples on average and the 29 term
participants 2.97; 8 of the 14 preterm participants have a single sample, against 4 of the 29 term
participants. The number of samples per participant therefore discriminates the outcome without
using the microbiome or any clinical variable.

This document measures that discrimination with the same outer folds and the same AUROC as the
models, and asks how much of it remains when only the weeks in which every participant was still
pregnant are counted.

---

## 1. What was computed

- **Score.** The number of samples per participant. The orientation is fixed in advance: fewer
  samples = preterm, the direction implied by an earlier end of pregnancy. Nothing is fitted (no
  parameter, no threshold), so there is no inner split: the score is evaluated directly on each
  outer test fold.
- **Folds.** The 5 outer folds of the models (subject level, stratified by outcome, seed 123). The
  script rebuilds them with the pipeline's code and checks that they are identical, subject by
  subject, to the folds of the 12 combinations in the main analysis of
  [`CLINICAL_AVAILABILITY_RESULTS.md`](CLINICAL_AVAILABILITY_RESULTS.md).
- **Metric.** AUROC with the fixed orientation of the engine (`roc_ptb()` in
  [`R/nested_cv.R`](../R/nested_cv.R)), reported as the mean ± SD across the 5 outer test folds,
  as for the models. The script checks that the same per-fold computation, applied to the saved
  predictions of the 12 combinations, reproduces the per-fold AUROC that the engine reported.
  Participants with the same number of samples are ties and count as half a correctly ordered
  pair. A test fold holds at most 18 preterm–term pairs, so the AUROC pooled over all participants
  of a variant is also given, with a DeLong 95% confidence interval.

Variants, with X = 24, 28 and 32 weeks:

| Variant | Participants | Samples counted |
|---|---|---|
| All samples | all 43 | every sample |
| All samples, participants who reached week X | delivery at week X or later | every sample |
| Samples before week X | delivery at week X or later | samples taken before week X |

Every preterm participant but one delivered at 32+5 weeks or later, so the participants who reached
weeks 24, 28 and 32 are the same 42 (13 preterm, 29 term): all except the participant whose
pregnancy ended at 20+0 weeks. Within the window before week X all 42 were still pregnant, so the
count before week X does not reflect the earlier end of the preterm pregnancies. A participant
whose first sample was taken at week X or later counts as 0 samples.

Gestational ages are recorded in weeks.days notation (24.4 = 24 weeks and 4 days). Comparisons with
a whole week, as here, give the same result in that notation as in decimal weeks.

## 2. Results

The second row is the same for weeks 24, 28 and 32, because the same 42 participants reached the
three weeks.

| Samples counted | Participants (preterm / term) | Samples | Mean samples per participant (preterm / term) | Participants with no sample in the window | AUROC, mean ± SD across folds | Pooled AUROC (95% CI) |
|---|---|---|---|---|---|---|
| All samples, all participants | 43 (14 / 29) | 110 | 1.71 / 2.97 | 0 | 0.771 ± 0.160 | 0.781 (0.63–0.93) |
| All samples, the 42 who reached week 24, 28 or 32 | 42 (13 / 29) | 109 | 1.77 / 2.97 | 0 | 0.738 ± 0.199 | 0.769 (0.61–0.92) |
| Before week 24 | 42 (13 / 29) | 48 | 0.85 / 1.28 | 11 | 0.636 ± 0.220 | 0.619 (0.45–0.79) |
| Before week 28 | 42 (13 / 29) | 64 | 1.23 / 1.66 | 4 | 0.622 ± 0.224 | 0.626 (0.45–0.80) |
| Before week 32 | 42 (13 / 29) | 82 | 1.46 / 2.17 | 2 | 0.664 ± 0.193 | 0.678 (0.51–0.84) |

AUROC by outer test fold. The participant whose pregnancy ended at 20+0 weeks is in fold 1, which
holds 9 participants (3 preterm) with all samples and 8 (2 preterm) in the other variants. The last
column is the best model of the main analysis (elastic net, data-driven clinical set, ANCOM-BC2
taxa), evaluated with all samples.

| Fold | Test participants (preterm) | All samples | All samples, the 42 | Before week 24 | Before week 28 | Before week 32 | Best model |
|---|---|---|---|---|---|---|---|
| 1 | 9 (3) / 8 (2) | 0.667 | 0.500 | 0.583 | 0.500 | 0.500 | 0.944 |
| 2 | 9 (3) | 0.889 | 0.889 | 0.806 | 0.889 | 0.944 | 0.667 |
| 3 | 9 (3) | 0.556 | 0.556 | 0.472 | 0.389 | 0.500 | 1.000 |
| 4 | 9 (3) | 0.944 | 0.944 | 0.917 | 0.833 | 0.778 | 0.778 |
| 5 | 7 (2) | 0.800 | 0.800 | 0.400 | 0.500 | 0.600 | 0.400 |

Participants by number of samples (term / preterm):

| Samples | All samples (43) | Before week 24 (42) | Before week 28 (42) | Before week 32 (42) |
|---|---|---|---|---|
| 0 | — | 6 / 5 | 2 / 2 | 1 / 1 |
| 1 | 4 / 8 | 13 / 5 | 12 / 7 | 8 / 7 |
| 2 | 5 / 3 | 6 / 3 | 10 / 3 | 10 / 3 |
| 3 | 12 / 2 | 4 / 0 | 4 / 1 | 6 / 2 |
| 4 | 5 / 1 | — | 1 / 0 | 3 / 0 |
| 5 | 2 / 0 | — | — | 1 / 0 |
| 6 | 1 / 0 | — | — | — |

## 3. What remains when every participant is still pregnant

Among the same 42 participants, counting only the samples taken while all of them were still
pregnant lowers the AUROC from 0.738 ± 0.199 (pooled 0.769) to 0.636, 0.622 and 0.664 before weeks
24, 28 and 32 (pooled 0.619, 0.626 and 0.678). **Between 44% and 69% of the margin above 0.5
remains**, depending on the week and on whether the per-fold mean or the pooled AUROC is used. The
pooled 95% intervals include 0.5 before weeks 24 and 28 and start at 0.51 before week 32, and the
per-fold SDs are about 0.2. With 42 participants, the remaining signal cannot be told apart from
chance, but it does not disappear.

Part of the discrimination of the full count is therefore the earlier end of the preterm
pregnancies. The point estimates also indicate that preterm participants had fewer samples during
weeks in which every participant was pregnant (1.46 against 2.17 samples before week 32). This
analysis does not establish why. The median gestational age at the first sample is 19+5 weeks for
the 13 preterm participants and 18+5 weeks for the 29 term participants.

Leaving out the participant whose pregnancy ended at 20+0 weeks (first row against second) changes
the pooled AUROC little (0.781 to 0.769). The per-fold mean falls more (0.771 to 0.738) because that
participant is in fold 1, which goes from 0.667 to 0.500 without them.

## 4. Comparison with the models

The models are evaluated with all the samples of every participant, so their reference is the
first row: **0.771 ± 0.160**. None of the 12 combinations of the main analysis has a higher mean
AUROC. The best, elastic net with the data-driven clinical set and ANCOM-BC2 taxa, is 0.758 ±
0.240, and the mean of the 12 combinations is 0.527. Fold by fold, the best model is higher in
folds 1 and 3 and the count in folds 2, 4 and 5.

How to read this comparison:

- The total number of samples is known only once sampling has ended, which for a preterm
  participant is when the pregnancy ends. It is a reference for how much discrimination the
  sampling pattern alone yields, not a predictor that could be used at a visit.
- The models do not receive the number of samples. But the samples of preterm participants cover a
  shorter part of pregnancy, and gestational age at the visit is one of the models' inputs
  (Approaches 1 and 2, and one outer fold of Approach 3). Whether the models' AUROC depends on the
  same sampling pattern is not measured here. Measuring it requires evaluating the models on the
  samples taken before a fixed week, among the participants who reached it (a landmark design).
  That is a separate analysis.

## 5. Files and reproducing these numbers

[`results/sample_count_reference/`](../results/sample_count_reference/) holds the numbers of this
document:

- `auroc_summary.csv`: one row per variant, with participants, samples, the per-fold mean and SD
  of the AUROC and the pooled AUROC with its 95% CI;
- `auroc_by_fold.csv`: one row per variant and outer test fold;
- `samples_per_participant.csv`: participants by number of samples and outcome, per variant.

The script needs the cohort data in `data/raw/`:

```bash
PTB_RESULTS_DIR=results/sample_count_reference Rscript scripts/sample_count_reference.R
```

It takes about 10 seconds (Windows 11, R 4.4.2) and prints the checks of section 1 and the tables
of this document. The fold check uses the saved output of the main-analysis run
(`PTB_RUN_TAG=atvisit_2026-10-07`, see [`CLINICAL_AVAILABILITY_RESULTS.md`](CLINICAL_AVAILABILITY_RESULTS.md)
§7); without that run the script says that the check was skipped. `PTB_COUNT_WEEKS` changes the
weeks (default `24,28,32`).
