FROM condaforge/miniforge3:latest

# Install system dependencies and Snakemake
RUN conda install -y -c bioconda -c conda-forge \
    snakemake-minimal \
    bash \
    curl \
    openjdk=17.0.10 \
    nextflow=25.04.3 \
    && conda clean -afy

# Replace system bash with conda bash (prevents GLIBC mismatch)
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

# Setup pipeline
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .

CMD ["snakemake", "--help"]