FROM condaforge/mambaforge:latest

# Install all bioinformatics tools
RUN mamba install -y -c conda-forge -c bioconda \
    snakemake-minimal \
    bwa-mem2 \
    samtools \
    bcftools \
    freebayes \
    fastp \
    fastqc \
    bash \
    && mamba clean -afy

# Replace system bash with conda bash
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

WORKDIR /pipeline

# Copy pipeline files
COPY Snakefile .
COPY config.yaml .
