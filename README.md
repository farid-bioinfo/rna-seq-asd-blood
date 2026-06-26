# Blood Transcriptomics in Autism Spectrum Disorder
### Differential Expression Analysis — GSE212645

Differential expression analysis of human whole-blood RNA-seq data comparing 20
individuals with Autism Spectrum Disorder (ASD) to their 20 unaffected siblings.
The paired sibling design controls for shared family environment and genetic
background, isolating condition-specific transcriptional signal.

**Dataset:** [GSE212645](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212645)
&nbsp;·&nbsp; **Author:** [Farid Hakimi](https://linkedin.com/in/farid-hakimi) · London, UK

---

## Biological question

Blood transcriptomics is an accessible but noisy tissue for ASD research.
This analysis asks: are there reproducible differential expression signals in
whole blood between ASD individuals and neurotypical siblings, after correcting
for within-family correlation?

The paired design (one ASD + one unaffected sibling per family) is critical.
Ignoring family structure in a paired study inflates false positives by treating
related samples as independent. The final model ~ family + condition explicitly
accounts for this.

---

## Dataset

| Property | Detail |
|---|---|
| GEO accession | GSE212645 |
| Tissue | Human whole blood |
| Samples | 40 (20 ASD · 20 unaffected siblings) |
| Study design | Matched sibling pairs (1 ASD + 1 sibling per family) |
| Genes (raw) | 29,536 |
| Genes (post-filter) | 13,439 (>=10 counts in >=20 samples) |
| Covariates | Family ID, sex (collinear with family — excluded from final model) |

---

## Pipeline
data/raw/GSE212645_CountsMatrix.txt.gz

|

v

01_load_data.R

Reads count matrix · fetches GEO metadata via GEOquery

Aligns sample names · builds DESeqDataSet (design ~ condition)

|

v

02_deseq2_analysis.R

Adds family covariate · applies low-count filter

Final model: ~ family + condition

Contrast: ASD vs SIB · saves results CSV

|

v

03_visualisation.R

Volcano plots (nominal p < 0.05 threshold)

PCA coloured by condition and family

|

v

04_go_enrichment.R

Gene symbol to Entrez ID conversion via bitr

GO Biological Process enrichment (clusterProfiler)

BH-adjusted p < 0.05 · dotplot output
---

## Key analytical decisions

**Why ~ family + condition and not ~ sex + condition?**
Sex and family ID are collinear in this dataset — each family has one male and
one female sibling. Including both would cause rank deficiency in the design
matrix. Family ID captures the sex effect and additionally controls for shared
genetic background, so it is the stronger covariate.

**Why no FDR-significant genes?**
799 genes reach nominal significance (raw p < 0.05); none survive
Benjamini-Hochberg correction (adjusted p < 0.1). This is expected:
blood transcriptomics for neurological phenotypes has a low signal-to-noise
ratio, and n = 20 pairs is modest. The result is biologically interpretable,
not a pipeline failure.

**Why use p < 0.05 without FDR correction for GO enrichment?**
With no FDR-significant genes, applying a strict dual threshold (padj < 0.05,
|LFC| > 1) reduces the gene set to a small but high-confidence list. GO
enrichment on this set is exploratory and should be interpreted cautiously —
it identifies directional hypotheses for future larger studies, not confirmed
pathways.

---

## Results

| Metric | Value |
|---|---|
| Genes tested | 13,439 |
| Nominally significant (p < 0.05) | 799 |
| FDR-significant (padj < 0.1) | 0 |
| Top hits (nominal) | ANK1, FKBP5, ESPN, HBQ1 |
| GO enrichment input | 53 genes (p < 0.05, |LFC| > 0.5) |
| Enriched GO BP terms | 12 |

The absence of FDR hits is consistent with published blood transcriptomics
studies in ASD at comparable sample sizes. Results are exploratory and would
require replication in an independent cohort.

---

## Repository structure
rna-seq-asd-blood/

├── data/

│   └── raw/

│       └── GSE212645_CountsMatrix.txt.gz

├── scripts/

│   └── R/

│       ├── 01_load_data.R

│       ├── 02_deseq2_analysis.R

│       ├── 03_visualisation.R

│       └── 04_go_enrichment.R

├── results/

│   ├── tables/

│   │   └── DESeq2_results_family_corrected.csv

│   └── figures/

│       ├── volcano_plot.png

│       ├── volcano_plot_labelled.png

│       └── 04_go_dotplot.png

├── docs/

├── notebooks/

├── CLAUDE.md

└── README.md
---

## How to reproduce

**Requirements:** R >= 4.4, internet access (GEO metadata fetch)

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "GEOquery", "clusterProfiler",
                       "org.Hs.eg.db", "enrichplot"))
install.packages(c("ggplot2", "ggrepel"))
```

```r
setwd("path/to/rna-seq-asd-blood")
source("scripts/R/01_load_data.R")
source("scripts/R/02_deseq2_analysis.R")
source("scripts/R/03_visualisation.R")
source("scripts/R/04_go_enrichment.R")
```

All outputs are written to results/tables/ and results/figures/.

---
## Reproducibility — Docker

The full pipeline is containerised and available on Docker Hub.

```bash
docker pull faridbioinfo2026/rna-seq-asd-pipeline:latest
docker run --rm -v $(pwd)/results:/pipeline/results faridbioinfo2026/rna-seq-asd-pipeline:latest
```

[![Docker Hub](https://img.shields.io/badge/Docker%20Hub-faridbioinfo2026-blue?logo=docker)](https://hub.docker.com/r/faridbioinfo2026/rna-seq-asd-pipeline)
Pipeline outputs (figures and results tables) are archived on AWS S3:
`s3://rna-seq-asd-blood-farid/results/` (eu-west-2)
---

## Stack

| Tool | Role |
|---|---|
| R 4.6 | Analysis language |
| DESeq2 | Differential expression (negative binomial GLM) |
| GEOquery | Programmatic metadata retrieval from NCBI GEO |
| ggplot2 · ggrepel | Volcano and PCA visualisation |
| clusterProfiler | GO and pathway enrichment |
| org.Hs.eg.db | Human gene annotation (symbol to Entrez mapping) |
| Snakemake | Workflow manager (pipeline orchestration) |
| Docker | Containerisation and reproducible execution |
| AWS S3 | Cloud storage for pipeline outputs (eu-west-2) |

---

## Background

This project sits within a broader research interest in the gut-brain axis and
neurodevelopmental conditions. Blood transcriptomics is one of several data
modalities used to characterise ASD biology; this analysis is a computational
complement to a parallel systematic review examining probiotic and
microbiota-targeted interventions in ASD, currently under review at the
Journal of Autism and Developmental Disorders.

---

## Author

**Farid Hakimi** — MSc Bioinformatics & Biostatistics (9.43/10) · MSc Biomedical
Science (Medical Microbiology) · Former Healthcare Scientist, Public Health England

[LinkedIn](https://linkedin.com/in/farid-hakimi) · [GitHub](https://github.com/farid-bioinfo) · London, UK

*Open to bioinformatics analyst and computational biology roles at research
institutions and biotech.*
