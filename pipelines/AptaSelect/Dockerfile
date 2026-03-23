FROM condaforge/miniforge3:latest

# Install Snakemake and Python (Python 3.11+ ships with str | None union types)
RUN conda install -y -c conda-forge -c bioconda \
    snakemake-minimal \
    bash \
    python=3.11 \
    && conda clean -afy

# Replace system bash with conda bash (prevents GLIBC mismatch)
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

# Pipeline files
WORKDIR /pipeline
COPY Snakefile  .
COPY config.yaml .
COPY scripts/   scripts/

CMD ["snakemake", "--cores", "4", "--nolock"]