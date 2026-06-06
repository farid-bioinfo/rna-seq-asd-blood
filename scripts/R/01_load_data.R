# =========================================================================
# 01_load_data.R
# Purpose: Load GSE212645 counts matrix and fetch GEO metadata
# Author: Farid Hakimi
# Date: June 2026
# =========================================================================

# --- Load Libraries ---
library(DESeq2)
library(GEOquery)

# --- Load counts matrix ---
counts <- read.table(
    gzfile("data/raw/GSE212645_CountsMatrix.txt.gz"),
    header = TRUE,
    row.names = 1,
    sep = "\t"
)

# --- Inspect counts matrix ---
dim(counts)
head(rownames(counts))
colnames(counts)
sum(is.na(counts))

# --- Fetch metadata from GEO ---
gse <- getGEO("GSE212645", GSEMatrix = TRUE)
metadata <- pData(gse[[1]])

# --- Build coldata for DESeq2 ---
coldata <- data.frame(
    row.names = metadata$title,
    condition = factor(metadata$"genotype:ch1")
)

# --- Reorder to match counts ---
coldata <- coldata[colnames(counts), , drop = FALSE]

# --- Verify ---
table(coldata$condition)
all(rownames(coldata) == colnames(counts))
