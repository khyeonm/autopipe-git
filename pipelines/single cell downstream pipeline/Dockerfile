FROM condaforge/miniforge3:latest

# Install core tools via conda
RUN conda install -y -c conda-forge -c bioconda \
    snakemake \
    bash \
    python=3.11 \
    && conda clean -afy

# Replace system bash with conda bash (prevents GLIBC mismatch)
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

# Install uv for fast Python package installation
RUN pip install uv

# Install scanpy and dependencies via uv
RUN uv pip install --system \
    "scanpy[leiden]" \
    leidenalg \
    scikit-misc \
    pyyaml \
    matplotlib \
    anndata \
    h5py

# Setup pipeline
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/

CMD ["snakemake", "--cores", "4", "-p"]