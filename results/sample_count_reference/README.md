# Number of samples per participant as a reference

Numbers behind [`docs/SAMPLE_COUNT_REFERENCE.md`](../../docs/SAMPLE_COUNT_REFERENCE.md): the AUROC
of the number of samples per participant (fewer samples = preterm, orientation fixed in advance),
evaluated on the same 5 outer folds as the models. Column `variant`:

| `variant` | Participants | Samples counted |
|---|---|---|
| `all` | all 43 | every sample |
| `all_reached_X` | delivery at week X or later | every sample |
| `before_X` | delivery at week X or later | samples taken before week X |

`week` is X (24, 28 or 32; empty for `all`).

- `auroc_summary.csv`: one row per variant. Number of participants (preterm, term), samples
  counted, mean samples per participant by outcome, participants with no sample in the window, the
  mean and SD of the AUROC across the outer folds, and the AUROC pooled over all participants of
  the variant with its DeLong 95% CI.
- `auroc_by_fold.csv`: one row per variant and outer fold, with the participants, preterm and
  term participants and samples in the test fold, and its AUROC.
- `samples_per_participant.csv`: number of term and preterm participants with each number of
  samples, per variant.

All three files are written by
[`scripts/sample_count_reference.R`](../../scripts/sample_count_reference.R); the command is in
section 5 of the document.
