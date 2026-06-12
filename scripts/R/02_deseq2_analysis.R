# 02_deseq2_analysis.R
library(DESeq2)
library(GEOquery)

gse <- getGEO("GSE212645", GSEMatrix = TRUE)
meta <- pData(phenoData(gse[[1]]))[, c("title", "Sex:ch1", "family:ch1")]
colnames(meta) <- c("sample", "sex", "family")
meta$sex    <- gsub("Sex: ", "", meta$sex)
meta$family <- gsub("family: ", "", meta$family)
rownames(meta) <- meta$sample

colData(dds)$family <- factor(meta[colnames(dds), "family"])

keep <- rowSums(counts(dds) >= 10) >= 20
dds  <- dds[keep, ]
cat("Genes remaining after filtering:", nrow(dds), "\n")

design(dds) <- ~ family + condition
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "ASD", "SIB"))
res <- res[order(res$padj), ]
summary(res)

res_df      <- as.data.frame(res)
res_df$gene <- rownames(res_df)
write.csv(res_df, "results/tables/DESeq2_results_family_corrected.csv", row.names = FALSE)
cat("Done\n")
