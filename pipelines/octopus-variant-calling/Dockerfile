FROM condaforge/mambaforge:latest

# Install main tools first
RUN mamba install -y -c conda-forge -c bioconda \
    snakemake-minimal \
    bwa-mem2 \
    samtools \
    bcftools \
    fastp \
    fastqc \
    bash \
    && mamba clean -afy

# Install octopus in a separate step with flexible deps
RUN mamba create -y -n octopus -c conda-forge -c bioconda octopus \
    && mamba clean -afy

# Create wrapper script so octopus is accessible from main env
RUN echo '#!/bin/bash' > /usr/local/bin/octopus && \
    echo 'exec /opt/conda/envs/octopus/bin/octopus "$@"' >> /usr/local/bin/octopus && \
    chmod +x /usr/local/bin/octopus

# Replace system bash with conda bash
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

WORKDIR /pipeline

COPY Snakefile .
COPY config.yaml .
