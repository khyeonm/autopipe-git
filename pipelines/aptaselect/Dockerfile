FROM condaforge/miniforge3:latest

# Install snakemake and bash via conda
RUN conda install -y -c bioconda -c conda-forge \
    snakemake \
    bash \
    && conda clean -afy

# Replace system bash with conda bash (prevents GLIBC mismatch)
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

# Install Python packages via pip
RUN pip install uv && uv pip install --system pyyaml

# Setup pipeline
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/

CMD ["snakemake", "--cores", "1"]