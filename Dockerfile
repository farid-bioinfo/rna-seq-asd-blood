FROM rocker/r-ver:4.6.0

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install snakemake==9.23.1 --break-system-packages

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(c('DESeq2', 'clusterProfiler', 'org.Hs.eg.db'), ask=FALSE)"
RUN R -e "install.packages(c('ggplot2', 'ggrepel', 'RColorBrewer', 'pheatmap', 'GEOquery', 'tidyverse'), repos='https://cloud.r-project.org')"

WORKDIR /pipeline

COPY . .
