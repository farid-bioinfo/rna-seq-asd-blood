rule all:
    input:
        "data/processed/dds_final.rds",
        "results/tables/DESeq2_results_family_corrected.csv",
        "results/figures/volcano_plot.png",
        "results/figures/volcano_plot_labelled.png",
        "results/figures/04_go_dotplot.png"


rule load_data:
    output:
        "data/processed/dds_initial.rds",
        "data/processed/metadata.csv"
    shell:
        "Rscript scripts/R/01_load_data.R"


rule deseq2_analysis:
    input:
        "data/processed/dds_initial.rds",
        "data/processed/metadata.csv"
    output:
        "data/processed/dds_final.rds",
        "results/tables/DESeq2_results_family_corrected.csv"
    shell:
        "Rscript scripts/R/02_deseq2_analysis.R"


rule visualisation:
    input:
        "data/processed/dds_final.rds",
        "results/tables/DESeq2_results_family_corrected.csv"
    output:
        "results/figures/volcano_plot.png",
        "results/figures/volcano_plot_labelled.png"
    shell:
        "Rscript scripts/R/03_visualisation.R"


rule go_enrichment:
    input:
        "results/tables/DESeq2_results_family_corrected.csv"
    output:
        "results/figures/04_go_dotplot.png"
    shell:
        "Rscript scripts/R/04_go_enrichment.R"

