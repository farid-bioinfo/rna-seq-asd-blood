# 01_load_data.R
# Input:  data/raw/GSE212645_CountsMatrix.txt.gz
# Output: data/processed/dds_initial.rds
#         data/processed/metadata.csv

library(DESeq2)
library(GEOquery)

# ── 1. Load count matrix ──────────────────────────────────────────────────────
counts <- read.table(
  gzfile("data/raw/GSE212645_CountsMatrix.txt.gz"),
  header   = TRUE,
  row.names = 1,
  sep      = "\t"
)
cat("Counts loaded:", nrow(counts), "genes x", ncol(counts), "samples\n")

# ── 2. Fetch metadata from GEO ────────────────────────────────────────────────
gse      <- getGEO("GSE212645", GSEMatrix = TRUE)
metadata <- pData(gse[[1]])

# ── 3. Build coldata ──────────────────────────────────────────────────────────
coldata <- data.frame(
  row.names = metadata$title,
  condition = factor(metadata$"genotype:ch1")
)
coldata <- coldata[colnames(counts), , drop = FALSE]

cat("Samples:", ncol(counts), "| Genes:", nrow(counts), "\n")
print(table(coldata$condition))

# ── 4. Build DESeqDataSet ─────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ condition
)
cat("DESeq2 object created\n")

# ── 5. Save outputs to disk ───────────────────────────────────────────────────
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

saveRDS(dds, "data/processed/dds_initial.rds")
cat("Saved: data/processed/dds_initial.rds\n")

write.csv(metadata, "data/processed/metadata.csv", row.names = FALSE)
cat("Saved: data/processed/metadata.csv\n")

