FROM condaforge/miniforge3:latest

RUN conda install -y -c conda-forge -c bioconda \
        snakemake-minimal \
        python=3.11 \
        bash \
    && conda clean -afy

# Ensure conda bash is used to avoid GLIBC mismatch
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash \
 && ln -sf /opt/conda/bin/bash /bin/sh

WORKDIR /pipeline
COPY Snakefile   /pipeline/Snakefile
COPY config.yaml /pipeline/config.yaml
COPY scripts/    /pipeline/scripts/