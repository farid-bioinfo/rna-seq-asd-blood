# 02_deseq2_analysis.R
# Input:  data/processed/dds_initial.rds
#         data/processed/metadata.csv
# Output: data/processed/dds_final.rds
#         results/tables/DESeq2_results_family_corrected.csv

library(DESeq2)

# ── 1. Load from disk ─────────────────────────────────────────────────────────
dds      <- readRDS("data/processed/dds_initial.rds")
metadata <- read.csv("data/processed/metadata.csv")
cat("Loaded dds:", nrow(dds), "genes x", ncol(dds), "samples\n")

# ── 2. Add family covariate ───────────────────────────────────────────────────
metadata$family <- gsub("family: ", "", metadata$family.ch1)
rownames(metadata) <- metadata$title
colData(dds)$family <- factor(metadata[colnames(dds), "family"])

# ── 3. Filter low-count genes ─────────────────────────────────────────────────
keep <- rowSums(counts(dds) >= 10) >= 20
dds  <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

# ── 4. Run DESeq2 with family-corrected model ─────────────────────────────────
design(dds) <- ~ family + condition
dds <- DESeq(dds)

# ── 5. Extract results ────────────────────────────────────────────────────────
res    <- results(dds, contrast = c("condition", "ASD", "SIB"))
res    <- res[order(res$padj), ]
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
summary(res)

# ── 6. Save outputs to disk ───────────────────────────────────────────────────
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

saveRDS(dds, "data/processed/dds_final.rds")
cat("Saved: data/processed/dds_final.rds\n")

write.csv(res_df, "results/tables/DESeq2_results_family_corrected.csv",
          row.names = FALSE)
cat("Saved: results/tables/DESeq2_results_family_corrected.csv\n")