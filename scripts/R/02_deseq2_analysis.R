# 02_deseq2_analysis.R
# Goal: Filter low counts, run DESeq2, extract DE results

library(DESeq2)

# ── 1. Load the dds object from script 01 ────────────────────────────────────
source("scripts/R/01_load_data.R")

# ── 2. Filter low-count genes ─────────────────────────────────────────────────
# Keep genes with at least 10 counts across at least 20 samples (half the dataset)
keep <- rowSums(counts(dds) >= 10) >= 20
dds <- dds[keep, ]

cat("Genes remaining after filtering:", nrow(dds), "\n")

# ── 3. Run DESeq2 ─────────────────────────────────────────────────────────────
dds <- DESeq(dds)

# ── 4. Extract results (ASD vs sibling) ───────────────────────────────────────
res <- results(dds, contrast = c("condition", "ASD", "sibling"))
res <- res[order(res$padj), ]  # sort by adjusted p-value

summary(res)

getwd()

setwd("C:/path/to/rna-seq-asd-blood")  # adjust to your actual path

list.dirs("C:/Users/localuser", recursive = FALSE)

setwd("C:/Users/localuser/Documents/rna-seq-asd-blood")  # adjust if needed
getwd()  # confirm it
source("scripts/R/01_load_data.R")

setwd("C:/Users/localuser/rna-seq-asd-blood")
getwd()  # should confirm the path
source("scripts/R/01_load_data.R")
n

dds

keep <- rowSums(counts(dds) >= 10) >= 20
dds <- dds[keep, ]
cat("Genes remaining after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "ASD", "sibling"))
res <- res[order(res$padj), ]
summary(res)

# Check what DESeq2 named the comparison
resultsNames(dds)

# Extract results with correct order
res <- results(dds, contrast = c("condition", "sibling", "ASD"))
res <- res[order(res$padj), ]
summary(res)

res <- results(dds, contrast = c("condition", "SIB", "ASD"))
res <- res[order(res$padj), ]
summary(res)

# Check what condition labels actually exist in your dds object
colData(dds)$condition |> table()

# Look at the raw (unadjusted) p-values
sum(res$pvalue < 0.05, na.rm = TRUE)
sum(res$padj < 0.05, na.rm = TRUE)

# Look at the top genes regardless of significance
head(res, 10)

# Plot p-value distribution
hist(res$pvalue, breaks = 50, main = "P-value distribution", xlab = "p-value")

# Pull full metadata from GEO to see what covariates exist
library(GEOquery)
gse <- getGEO("GSE212645", GSEMatrix = TRUE)
pdata <- pData(phenoData(gse[[1]]))

# See all available columns
colnames(pdata)

# Look at key characteristics
head(pdata[, c("title", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")])

# Extract key covariates
meta_ext <- pdata[, c("Sex:ch1", "puberty:ch1", "family:ch1")]
colnames(meta_ext) <- c("sex", "puberty", "family")

# Clean up values
meta_ext$sex     <- gsub("Sex: ", "", meta_ext$sex)
meta_ext$puberty <- gsub("puberty: ", "", meta_ext$puberty)

# Align with dds sample order
meta_ext <- meta_ext[colnames(dds), ]

# Add to dds
colData(dds)$sex     <- meta_ext$sex
colData(dds)$puberty <- meta_ext$puberty

# Rebuild model with sex as covariate
design(dds) <- ~ sex + condition

# Rerun DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("condition", "SIB", "ASD"))
res <- res[order(res$padj), ]
summary(res)

# Replot p-value distribution
hist(res$pvalue, breaks = 50, main = "P-value distribution (sex-corrected)", xlab = "p-value")
# Check sex actually got added
colnames(colData(dds))

colData(dds)$sex    <- factor(meta_ext$sex)
colData(dds)$family <- factor(meta_ext$family)
colnames(colData(dds))

meta_ext <- pData(phenoData(gse[[1]]))[, c("Sex:ch1", "puberty:ch1", "family:ch1")]
colnames(meta_ext) <- c("sex", "puberty", "family")
meta_ext$sex    <- gsub("Sex: ", "", meta_ext$sex)
meta_ext$family <- gsub("family: ", "", meta_ext$family)
meta_ext <- meta_ext[colnames(dds), ]
head(meta_ext)

rownames(meta_ext)[1:5]
colnames(dds)[1:5]

meta_ext <- pData(phenoData(gse[[1]]))[, c("title", "Sex:ch1", "family:ch1")]
colnames(meta_ext) <- c("sample", "sex", "family")
meta_ext$sex    <- gsub("Sex: ", "", meta_ext$sex)
meta_ext$family <- gsub("family: ", "", meta_ext$family)
rownames(meta_ext) <- meta_ext$sample
head(meta_ext)

all(colnames(dds) %in% rownames(meta_ext))
colData(dds)$sex    <- factor(meta_ext[colnames(dds), "sex"])
colData(dds)$family <- factor(meta_ext[colnames(dds), "family"])
colnames(colData(dds))

dds3 <- DESeqDataSetFromMatrix(
  countData = counts(dds),
  colData   = colData(dds),
  design    = ~ family + sex + condition
)
dds3 <- dds3[keep, ]
dds3 <- DESeq(dds3)

  +   countData = counts(dds),
  +   colData   = colData(dds),
  +   design    = ~ family + sex + condition
  + )


table(meta_ext$family, meta_ext$sex)

dds3 <- DESeqDataSetFromMatrix(
  countData = counts(dds),
  colData   = colData(dds),
  design    = ~ family + condition
)
dds3 <- dds3[keep, ]
dds3 <- DESeq(dds3)

nrow(dds3)

res3 <- results(dds3, contrast = c("condition", "SIB", "ASD"))
res3 <- res3[order(res3$padj), ]
summary(res3)
hist(res3$pvalue, breaks = 50, main = "P-value distribution (family-corrected)", xlab = "p-value")

sum(res3$pvalue < 0.05, na.rm = TRUE)
sum(res3$padj < 0.05, na.rm = TRUE)
head(res3, 10)

res3_df <- as.data.frame(res3)
res3_df$gene <- rownames(res3_df)
write.csv(res3_df, "results/tables/DESeq2_results_family_corrected.csv", row.names = FALSE)


