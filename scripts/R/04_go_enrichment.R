# 04_go_enrichment.R
# Input:  results/tables/DESeq2_results_family_corrected.csv
# Output: results/figures/04_go_dotplot.png

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# ── 1. Load results from disk ─────────────────────────────────────────────────
res_df <- read.csv("results/tables/DESeq2_results_family_corrected.csv")
cat("Loaded:", nrow(res_df), "genes\n")

# ── 2. Apply dual threshold filter ───────────────────────────────────────────
sig_genes <- res_df[
  !is.na(res_df$pvalue) &
    res_df$pvalue < 0.05 &
    abs(res_df$log2FoldChange) > 0.5, ]

cat("Genes passing threshold:", nrow(sig_genes), "\n")

# ── 3. Convert gene symbols to Entrez IDs ────────────────────────────────────
entrez_ids <- bitr(sig_genes$gene,
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db)

cat("Genes mapped to Entrez IDs:", nrow(entrez_ids), "\n")

# ── 4. Run GO enrichment (Biological Process) ─────────────────────────────────
go_results <- enrichGO(
  gene          = entrez_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

cat("Enriched GO terms:", nrow(go_results@result[
  go_results@result$p.adjust < 0.05, ]), "\n")

# ── 5. Save dotplot ───────────────────────────────────────────────────────────
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

if (nrow(go_results) > 0) {
  p <- dotplot(go_results, showCategory = 20,
               title = "GO Enrichment — Biological Process")
  ggsave("results/figures/04_go_dotplot.png",
         plot = p, width = 10, height = 8, dpi = 150)
  cat("Saved: results/figures/04_go_dotplot.png\n")
} else {
  cat("No enriched GO terms found at this threshold\n")
}