FROM condaforge/miniforge3:latest

RUN conda install -y -c bioconda -c conda-forge \
    snakemake-minimal \
    python=3.11 \
    bash \
    && conda clean -afy

RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

COPY Snakefile /pipeline/Snakefile
COPY config.yaml /pipeline/config.yaml
COPY scripts/ /pipeline/scripts/

WORKDIR /pipeline