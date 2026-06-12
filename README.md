# RNA-seq Differential Expression Analysis
## ASD vs Unaffected Siblings — GSE212645

A bioinformatics portfolio project analysing human blood RNA-seq data to identify
differentially expressed genes between individuals with Autism Spectrum Disorder (ASD)
and their unaffected siblings.

---

## Dataset

| Property | Detail |
|---|---|
| GEO Accession | GSE212645 |
| Tissue | Human whole blood |
| Samples | 40 (20 ASD / 20 unaffected siblings) |
| Design | Matched family pairs (1 ASD + 1 sibling per family) |
| Genes | 29,536 (13,439 after low-count filtering) |

---

## Methods

### 1. Data loading
Count matrix downloaded from GEO. Sample metadata fetched programmatically
via GEOquery to ensure accurate group assignment.

### 2. Differential expression
DESeq2 with a paired family design to account for within-family correlation:
Family and sex are collinear in this dataset — only family was included.

### 3. Results
- No FDR-significant genes (expected for blood tissue at n=20 pairs)
- 799 nominally significant genes (raw p < 0.05)
- Top hits: ANK1, FKBP5, ESPN, HBQ1

### 4. GO enrichment
clusterProfiler applied to 53 high-confidence genes (p < 0.05, |LFC| > 0.5).
12 enriched Biological Process terms identified.

---

## Repository Structure
---

## Stack

R 4.6 · DESeq2 · GEOquery · ggplot2 · ggrepel · clusterProfiler · org.Hs.eg.db

---

## Interpretation

The absence of FDR-significant hits is consistent with the known low
signal-to-noise ratio of blood transcriptomics for ASD phenotype and the
modest sample size. Results should be interpreted as exploratory.
