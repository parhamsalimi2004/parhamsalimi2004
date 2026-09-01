# Bioinformatics Portfolio — Parham Salimi

A set of self-contained bioinformatics projects in cancer genomics. Each project uses real public data, states its dataset and methods explicitly, and validates its results against an independent, checkable reference (a published finding or the dataset's own ground-truth annotations) rather than reporting
numbers in isolation.

**GitHub:** [github.com/parhamsalimi2004](https://github.com/parhamsalimi2004/parhamsalimi2004)
**LinkedIn:** [linkedin.com/in/parham-salimi-13a781217](https://www.linkedin.com/in/parham-salimi-13a781217)

## Projects

| Project | Focus | Key result |
|---|---|---|
| [Mutation Landscape Explorer](projects/mutation-landscape.html) | Oncoplot of TCGA-LAML somatic mutations | Recovers FLT3/DNMT3A/NPM1 as top AML driver genes, matching known biology |
| [Survival Analysis Dashboard](projects/survival-analysis.html) | Kaplan-Meier survival curves + log-rank tests | TP53 mutation status significantly predicts survival (p < 0.0001) |
| [Driver Mutation Classifier](projects/driver-classifier.html) | ML: driver-gene vs. passenger mutation prediction, gene-agnostic features | Random Forest ROC-AUC 0.894; hotspot recurrence is the strongest predictor |
| [Variant Annotation Pipeline](projects/variant-annotator.html) | From-scratch gene/variant-type annotator | 100% variant-type accuracy; 100% gene accuracy where reference coverage exists |

Each project folder contains its own detailed `README.md` with full methods, dataset
provenance, and results — the table above is a summary; the per-project READMEs are the
primary technical documentation.

All four projects share the same underlying dataset (TCGA-LAML, Acute Myeloid Leukemia),
allowing findings from one to inform and cross-validate another — e.g. the top genes
identified in the Mutation Landscape project are the same genes tested in the Survival
Analysis and used to define the label in the Driver Mutation Classifier.

## Tech Stack

- **Languages:** Python
- **Core libraries:** pandas, NumPy, matplotlib, scikit-learn, lifelines
- **Data formats:** MAF (Mutation Annotation Format), TSV/CSV
- **Data sources:** TCGA (via the `maftools` R package's redistributed example data)

## Repository Structure

```
.
├── index.html                       # homepage
├── style.css                        # shared styles
├── assets/img/                      # generated figures
└── projects/
    ├── mutation-landscape/
    │   ├── README.md                 # project-specific documentation
    │   ├── analyze.py                # analysis script
    │   ├── requirements.txt
    │   └── data/tcga_laml.maf        # public TCGA-LAML MAF dataset
    ├── mutation-landscape.html       # project write-up page (site)
    ├── survival-analysis/
    │   ├── README.md
    │   ├── analyze.py
    │   ├── requirements.txt
    │   └── data/                     # MAF + clinical annotation file
    ├── survival-analysis.html
    ├── driver-classifier/
    │   ├── README.md
    │   ├── analyze.py
    │   ├── requirements.txt
    │   └── data/tcga_laml.maf
    ├── driver-classifier.html
    ├── variant-annotator/
    │   ├── README.md
    │   ├── analyze.py
    │   ├── requirements.txt
    │   └── data/tcga_laml.maf
    └── variant-annotator.html
```

## Running the Analyses Locally

Each project is independently reproducible:

```bash
cd projects/mutation-landscape && pip install -r requirements.txt && python analyze.py
cd ../survival-analysis          && pip install -r requirements.txt && python analyze.py
cd ../driver-classifier          && pip install -r requirements.txt && python analyze.py
cd ../variant-annotator          && pip install -r requirements.txt && python analyze.py
```

Each script regenerates its own figure(s) in `assets/img/`.

## Data Sources & Citation

All four projects use the TCGA-LAML (Acute Myeloid Leukemia) cohort, originally generated
by The Cancer Genome Atlas Research Network (*N Engl J Med* 2013), redistributed as public
example data by the [`maftools`](https://github.com/PoisonAlien/maftools) R/Bioconductor
package (Mayakonda et al., *Genome Research*, 2018). The driver gene panel used in two
projects is from Papaemmanuil et al., *N Engl J Med*, 2016. Full citations and per-dataset
details are in each project's own README.

## Roadmap

Two additional projects — a differential gene expression analysis and an end-to-end
FASTQ-to-VCF variant calling pipeline (Snakemake) are coming soon...

## Author

**Parham Salimi**
[GitHub](https://github.com/parhamsalimi2004/parhamsalimi2004) ·
[LinkedIn](https://www.linkedin.com/in/parham-salimi-13a781217)
