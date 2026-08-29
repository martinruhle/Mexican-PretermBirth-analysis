# Data Access

How to obtain the data behind this analysis, and what is already available without a request.

## Citation

Work using these data should cite:

```
Ruhle, M., et al. (2025). Leakage-aware machine learning reveals structured clinical and vaginal
microbiome patterns associated with preterm birth in a Mexican cohort.
Front. Glob. Women's Health 7:1799518. doi: 10.3389/fgwh.2026.1799518
```

Note that the analysis code has been updated since publication and some reported numbers changed;
see [`UPDATE_SINCE_PUBLICATION.md`](UPDATE_SINCE_PUBLICATION.md).

## Openly available

| What | Where |
|---|---|
| Raw 16S rRNA sequencing reads | NCBI Sequence Read Archive, BioProject **PRJNA1440471** |
| De-identified longitudinal clinical metadata | `data/metadata/metadata_eugenia_long.csv` in this repository |
| De-identified participant-level summary | `data/metadata/participant_data_clean.csv` |
| Variable dictionaries | `data/metadata/diccionario_variables_completo.csv` and `config/data_dictionary.csv` |
| Synthetic example dataset (runs the full pipeline) | `data/example/` |

The clinical metadata are de-identified and carry no direct identifiers.

## Not currently in this repository

The **genus-level abundance tables** (relative abundances and integer counts) are not currently
included. The pipeline expects them under `data/raw/`, which is untracked. `data/example/` exists so
that the pipeline can be run end to end without them.

To request them, contact the maintainer.

**Contact:** Martin Ruhle · martinruhle@gmail.com
Instituto Nacional de Medicina Genómica, Mexico City

**Please include:**

1. A short research proposal
2. Ethics approval from your institution
3. A signed data use agreement
4. CV of the principal investigator

**Expected timeline:** 4–8 weeks for review.

---

**Last updated:** 2026-08-10
