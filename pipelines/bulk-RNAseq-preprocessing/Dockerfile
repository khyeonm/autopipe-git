FROM condaforge/miniforge3:latest

# Install snakemake, nextflow, and basic tools via conda
RUN conda install -y -c bioconda -c conda-forge \
    snakemake-minimal \
    bash \
    openjdk=17 \
    nextflow=24.10.4 \
    && conda clean -afy

# Replace system bash with conda bash (prevents GLIBC mismatch)
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

# Setup pipeline
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .

CMD ["snakemake", "--help"]