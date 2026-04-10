library(metafor)
library(dplyr)
dat <- data.frame(
  study     = c("Liu 2019","Liu 2019","Santocchi 2020","Santocchi 2020",
                "Mazzone 2024","Mazzone 2024","Narula Khanna 2025"),
  outcome   = c("ABC","SRS","ADOS","Adaptive","SRS","Adaptive","SRS"),
  scale     = c("ABC-T","SRS","ADOS-CSS","VABS-II",
                "SRS Total","ABAS-2 Social","SRS-2 Total T"),
  direction = c("lower_better","lower_better","lower_better","higher_better",
                "lower_better","higher_better","lower_better"),
  n1i  = c(36,36,31,31,21,21,90),
  m1i  = c(14.67,132.77,6.19,67.39,80.70,69.70,70.49),
  sd1i = c(8.97,22.99,1.56,22.29,16.75,14.45,8.51),
  n2i  = c(35,35,32,32,22,22,90),
  m2i  = c(16.21,135.79,7.00,59.72,85.09,63.86,72.39),
  sd2i = c(10.11,25.79,1.80,16.38,17.14,13.50,7.03)
)

colnames(counts)
dat <- escalc(measure="SMD",
              m1i=m1i, sd1i=sd1i, n1i=n1i,
              m2i=m2i, sd2i=sd2i, n2i=n2i,
              data=dat, append=TRUE)

dat$yi <- ifelse(dat$direction == "higher_better", dat$yi * -1, dat$yi)

print(dat[, c("study","scale","yi","vi")])
res_overall <- rma(yi, vi, data=dat, method="REML",
                   slab=paste(study, scale))
summary(res_overall)

counts <- read.table(
  gzfile("C:/Users/localuser/rna-seq-asd-blood/data/raw/GSE212645_CountsMatrix.txt.gz"),
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

dim(counts)
# Check sample names (columns)
colnames(counts)

# Check first 6 gene names (rows)
head(rownames(counts))

# Check for any missing values
sum(is.na(counts))

head(rownames(counts))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
"n"
"n"
library(DESeq2)

# Build metadata table
sample_names <- colnames(counts)

condition <- ifelse(grepl("A$", sample_names), "ASD", "Sibling")

coldata <- data.frame(
  row.names = sample_names,
  condition = factor(condition)
)

# Check it
coldata

table(coldata$condition)

# Install GEOquery to fetch metadata
BiocManager::install("GEOquery")
library(GEOquery)

# Fetch metadata for GSE212645
gse <- getGEO("GSE212645", GSEMatrix = TRUE)
metadata <- pData(gse[[1]])

# Check what columns are available
colnames(metadata)
metadata <- pData(gse[[1]])
colnames(metadata)

table(metadata$"genotype:ch1")

# Build correct coldata from GEO metadata
coldata <- data.frame(
  row.names = metadata$title,
  condition = factor(metadata$"genotype:ch1")
)

# Verify
table(coldata$condition)

all(rownames(coldata) == colnames(counts))

rownames(coldata)
colnames(counts)

coldata <- coldata[colnames(counts), , drop = FALSE]

# Verify
all(rownames(coldata) == colnames(counts))


dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

dds
