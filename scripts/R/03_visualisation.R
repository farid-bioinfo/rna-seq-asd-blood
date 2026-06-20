# 03_visualisation.R
# Input:  results/tables/DESeq2_results_family_corrected.csv
# Output: results/figures/volcano_plot.png
#         results/figures/volcano_plot_labelled.png

library(ggplot2)
library(ggrepel)

# ── 1. Load results from disk ─────────────────────────────────────────────────
res_df <- read.csv("results/tables/DESeq2_results_family_corrected.csv")
cat("Loaded:", nrow(res_df), "genes\n")

# ── 2. Add significance labels ────────────────────────────────────────────────
res_df$significance <- "NS"
res_df$significance[res_df$pvalue < 0.05 & res_df$log2FoldChange > 0] <- "Up"
res_df$significance[res_df$pvalue < 0.05 & res_df$log2FoldChange < 0] <- "Down"
res_df$significance <- factor(res_df$significance, levels = c("Up", "Down", "NS"))

top_genes <- head(res_df[order(res_df$pvalue), ], 15)

# ── 3. Volcano plot ───────────────────────────────────────────────────────────
p1 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue),
                         colour = significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_colour_manual(values = c("Up" = "#E41A1C",
                                 "Down" = "#377EB8",
                                 "NS"  = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  labs(title    = "ASD vs Unaffected Siblings — GSE212645",
       subtitle = "Model: ~ family + condition  |  Nominal p < 0.05",
       x        = "Log2 Fold Change (ASD vs SIB)",
       y        = "-log10(p-value)",
       colour   = NULL) +
  theme_classic(base_size = 12)

# ── 4. Labelled volcano plot ──────────────────────────────────────────────────
p2 <- p1 +
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  size = 3, max.overlaps = 20)

# ── 5. Save outputs ───────────────────────────────────────────────────────────
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

ggsave("results/figures/volcano_plot.png",
       plot = p1, width = 8, height = 6, dpi = 150)
cat("Saved: results/figures/volcano_plot.png\n")

ggsave("results/figures/volcano_plot_labelled.png",
       plot = p2, width = 8, height = 6, dpi = 150)
cat("Saved: results/figures/volcano_plot_labelled.png\n")