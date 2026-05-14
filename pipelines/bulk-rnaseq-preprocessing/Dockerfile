FROM condaforge/miniforge3:latest

# Install snakemake and bash via conda
RUN conda install -y -c bioconda -c conda-forge \
    snakemake-minimal \
    bash \
    curl \
    openjdk=17 \
    && conda clean -afy

# Replace system bash with conda bash (prevents GLIBC mismatch)
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

# Install Nextflow (pinned to 25.04.3 for nf-core/rnaseq 3.22.2 compatibility)
ENV NXF_VER=25.04.3
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# Setup pipeline
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .

CMD ["snakemake", "--help"]