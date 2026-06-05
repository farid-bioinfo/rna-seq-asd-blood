# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Differential expression analysis of human blood RNA-seq data (GEO: GSE212645) comparing individuals with Autism Spectrum Disorder (ASD) to unaffected siblings (40 samples). The final model corrects for family structure using `~ family + condition`.

## Running the Analysis

All scripts must be run with the **working directory set to the project root** (`rna-seq-asd-blood/`). Scripts are numbered and meant to be sourced in order:

```r
setwd("C:/Users/localuser/rna-seq-asd-blood")

source("scripts/R/01_load_data.R")   # loads counts + builds dds object
source("scripts/R/02_deseq2_analysis.R")  # filters, runs DESeq2, saves results
# visualisation script: scripts-R-03_visualisation.R (in project root)
```

To run an individual script in RStudio or Rscript:
```bash
Rscript -e 'setwd("C:/Users/localuser/rna-seq-asd-blood"); source("scripts/R/02_deseq2_analysis.R")'
```

## Key R Packages

- **DESeq2** (Bioconductor) — differential expression
- **GEOquery** (Bioconductor) — fetches sample metadata from NCBI GEO
- **ggplot2**, **ggrepel** — visualisation
- **metafor** — meta-analysis (separate exploratory work in `01_load_data.R`)

Install Bioconductor packages:
```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "GEOquery"))
install.packages(c("ggplot2", "ggrepel"))
```

## Pipeline Architecture

```
data/raw/GSE212645_CountsMatrix.txt.gz
        │
        ▼
01_load_data.R
  - Reads count matrix (genes × samples)
  - Fetches metadata from GEO (sex, puberty, family covariates)
  - Aligns sample names between count matrix and GEO metadata
  - Builds DESeqDataSet (dds) with design ~ condition
        │
        ▼
02_deseq2_analysis.R
  - Filters low-count genes (≥10 counts in ≥20 samples)
  - Iterates on model design: plain → sex-corrected → family-corrected
  - Final model: dds3 with design ~ family + condition
  - Extracts results: contrast = c("condition", "SIB", "ASD")
  - Saves: results/tables/DESeq2_results_family_corrected.csv
        │
        ▼
scripts-R-03_visualisation.R  (project root)
  - Reads DESeq2_results_family_corrected.csv
  - Volcano plots coloured by nominal significance (p<0.05, no FDR threshold met)
  - Saves: results/figures/volcano_plot.png, volcano_plot_labelled.png
```

## Important Analysis Notes

- **No FDR-significant genes** were found; results use nominal p<0.05 as the significance threshold in visualisations
- The contrast direction is `SIB vs ASD` (positive log2FC = higher in siblings)
- `scripts/R/03_visualisation.R` is empty — the working visualisation script is `scripts-R-03_visualisation.R` in the project root
- `01_load_data.R` also contains unrelated metafor meta-analysis code at the top (exploratory work, not part of the RNA-seq pipeline)
- GEO metadata fetch requires internet access; `gse` object is not cached locally
