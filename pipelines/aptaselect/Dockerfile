FROM condaforge/miniforge3:latest

# Install snakemake and bash
RUN conda install -y -c bioconda -c conda-forge \
    snakemake-minimal \
    bash \
    python=3.11 \
    && conda clean -afy

# Replace system bash with conda bash (prevents GLIBC mismatch)
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

# Setup pipeline
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/

CMD ["snakemake", "--help"]