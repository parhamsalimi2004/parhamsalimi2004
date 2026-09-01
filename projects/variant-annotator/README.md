# Variant Annotation Pipeline

A from-scratch variant annotator that maps raw genomic positions to genes and predicts
variant type, built and validated against a professionally annotated dataset used as ground
truth.

## Overview

Raw variant callers output only genomic coordinates and alleles (e.g. `chr17:7578403, C→T`)
— they have no knowledge of gene boundaries or functional consequence. Translating that into
something interpretable is the job of an annotation pipeline (production tools: VEP,
ANNOVAR). This project implements a lightweight version of that process and evaluates it
honestly against the same dataset's original professional annotations, which happen to be
available as ground truth.

**Design constraint (by intent):** the input is deliberately stripped down to only what a
real variant caller would output — chromosome, position, reference/alternate allele, sample
ID. Gene name and functional consequence are hidden from the pipeline and used only for
post-hoc validation.

## Dataset

| **Cohort** | TCGA-LAML — Acute Myeloid Leukemia, The Cancer Genome Atlas (same cohort as the other MAF-based projects in this portfolio) |
| **Source file** | `tcga_laml.maf`, from [`maftools`](https://github.com/PoisonAlien/maftools) |
| **Split** | By patient (not by mutation row): 135 patients / 1,495 mutations for building the reference table; 58 patients / 712 mutations held out for evaluation |
| **Driver gene panel** (for the driver-flag stage) | Same 25-gene curated AML panel used in the [Driver Mutation Classifier](../driver-classifier/) project, from Papaemmanuil et al. 2016 |

## Pipeline Stages

1. **Gene mapping** — build a `(chromosome, gene) → position range` lookup table from the
   training-split patients only (min/max observed mutation position per gene, padded 2kb),
   then classify each held-out variant's gene by which range its position falls inside.
2. **Variant typing** — classify SNP / insertion / deletion purely from comparing
   reference vs. alternate allele string lengths.
3. **Frame-consequence prediction** — for indels, whether the length change is a multiple
   of 3 (in-frame) or not (frameshift).
4. **Driver flag** — cross-reference the predicted gene against the driver gene panel.

## Results

| Metric | Value |
|---|---|
| Reference coverage (test variants whose gene had a training-set entry) | 21.3% (152/712) |
| Gene assignment accuracy, **given** a reference match was found | **100%** |
| Overall gene assignment accuracy (uncovered counted as wrong) | 21.3% |
| Variant type classification accuracy | 100% (after a bug fix — see below) |
| Driver-flag agreement with ground truth | 99.3% |

**Interpretation:** the low headline gene-accuracy figure is a **coverage problem, not a
logic problem**. This cohort has ~1,600 unique mutated genes across only 2,207 total
mutations, so most genes are hit exactly once in the entire dataset — a patient-level
train/test split will always leave many test-set genes that training patients never
mutated. When the reference table did have an entry for a gene, position-based lookup was
correct 100% of the time. This is the direct, disclosed limitation of approximating a gene
reference from observed mutation positions rather than using a genome-wide transcript
database (as VEP/ANNOVAR do) — the root cause is identified and stated explicitly rather
than the headline number reported in isolation.

## A Bug Found and Fixed

Initial variant-type classification misclassified 15 deletions and 14 insertions as SNPs.
Root cause: MAF format represents a **zero-length allele** (pure insertion or deletion) with
the literal character `"-"`, not an empty string — `len("-")` evaluates to `1` in Python, so
a deletion like `ref="G", alt="-"` was read as two equal-length one-character alleles and
misclassified as a same-length substitution. Fixed by explicitly normalizing `"-"` to a
zero-length allele before comparing lengths, which brought variant-type accuracy to 100%.

## Repository Structure

```
variant-annotator/
├── analyze.py           # gene-mapping, variant-typing, and validation pipeline
├── requirements.txt
└── data/
    └── tcga_laml.maf    # TCGA-LAML mutation data (MAF format)
```

Output (`variant_annotator.png`: annotation outcome breakdown + variant-type confusion
matrix) is written to the shared `assets/img/` directory used by the portfolio site.

## Reproduction

```bash
pip install -r requirements.txt
python analyze.py
```

## Limitations & Future Work

- Gene reference is built from observed mutation positions in the training split, not a
  real genome-wide transcript annotation (RefSeq/Ensembl) — this is the direct cause of the
  21.3% coverage ceiling and the primary planned improvement.
- Frame-consequence prediction uses raw genomic position spans, with no awareness of actual
  exon/intron boundaries.
- Planned extension: substitute a real genome-wide coordinate reference (RefSeq or Ensembl
  GTF) to eliminate the coverage gap, and add exon-aware frame prediction.

## Data Citation

Mutation data originally generated by The Cancer Genome Atlas Research Network (*N Engl J
Med* 2013); redistributed as example data by the `maftools` package (Mayakonda et al.,
Genome Research 2018). Driver gene panel: Papaemmanuil E, et al. *Genomic Classification
and Prognosis in Acute Myeloid Leukemia.* N Engl J Med. 2016.

## Author

Parham Salimi