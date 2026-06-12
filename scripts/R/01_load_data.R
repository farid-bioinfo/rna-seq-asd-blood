# 01_load_data.R
library(DESeq2)
library(GEOquery)

counts <- read.table(
  gzfile("data/raw/GSE212645_CountsMatrix.txt.gz"),
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

gse      <- getGEO("GSE212645", GSEMatrix = TRUE)
metadata <- pData(gse[[1]])

coldata <- data.frame(
  row.names = metadata$title,
  condition = factor(metadata$"genotype:ch1")
)

coldata <- coldata[colnames(counts), , drop = FALSE]

cat("Samples:", ncol(counts), "| Genes:", nrow(counts), "\n")
cat("Conditions:\n")
print(table(coldata$condition))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ condition
)

cat("DESeq2 object created successfully\n")
