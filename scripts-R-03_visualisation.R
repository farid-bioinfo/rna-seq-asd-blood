library(DESeq2)
library(ggplot2)

# Load results
res3_df <- read.csv("results/tables/DESeq2_results_family_corrected.csv")

# Add significance labels
res3_df$significance <- "Not significant"
res3_df$significance[res3_df$pvalue < 0.05] <- "Nominal (p<0.05)"

# Volcano plot
ggplot(res3_df, aes(x = log2FoldChange, y = -log10(pvalue), colour = significance)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_colour_manual(values = c("grey70", "firebrick")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  labs(
    title = "ASD vs Sibling — Differential Expression",
    x = "Log2 Fold Change (SIB vs ASD)",
    y = "-log10(p-value)",
    colour = ""
  ) +
  theme_minimal()

ggsave("results/figures/volcano_plot.png", width = 8, height = 6, dpi = 300)


scale_colour_manual(values = c("Not significant" = "grey70", "Nominal (p<0.05)" = "firebrick")) +

  table(res3_df$significance)  
table(res3_df$significance)


ggplot(res3_df, aes(x = log2FoldChange, y = -log10(pvalue), colour = significance)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_colour_manual(values = c("Not significant" = "grey70", "Nominal (p<0.05)" = "firebrick")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  labs(
    title = "ASD vs Sibling — Differential Expression",
    x = "Log2 Fold Change (SIB vs ASD)",
    y = "-log10(p-value)",
    colour = ""
  ) +
  theme_minimal()

ggsave("results/figures/volcano_plot.png", width = 8, height = 6, dpi = 300)

library(ggrepel)

# Top 15 genes by p-value
top15 <- head(res3_df[order(res3_df$pvalue), ], 15)

ggplot(res3_df, aes(x = log2FoldChange, y = -log10(pvalue), colour = significance)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_colour_manual(values = c("Not significant" = "grey70", "Nominal (p<0.05)" = "firebrick")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_label_repel(data = top15, aes(label = gene), size = 3, max.overlaps = 20) +
  labs(
    title = "ASD vs Sibling — Differential Expression",
    x = "Log2 Fold Change (SIB vs ASD)",
    y = "-log10(p-value)",
    colour = ""
  ) +
  theme_minimal()

ggsave("results/figures/volcano_plot_labelled.png", width = 8, height = 6, dpi = 300)

install.packages("ggrepel")
library(ggrepel)
