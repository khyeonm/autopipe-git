FROM condaforge/miniforge3:latest

# Install core dependencies and Python 3.11 (more compatible with bioinformatics packages)
RUN conda install -y -c bioconda -c conda-forge \
    snakemake-minimal \
    bash \
    curl \
    python=3.11 \
    && conda clean -afy

# Replace system bash with conda bash (prevents GLIBC mismatch)
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

# Install uv for fast Python package installation
RUN pip install uv

# Install scanpy and all dependencies (including leidenalg)
RUN uv pip install --system \
    scanpy==1.9.5 \
    anndata==0.9.2 \
    numpy==1.24.3 \
    pandas==2.0.3 \
    scikit-learn==1.3.0 \
    umap-learn==0.5.7 \
    igraph==0.10.8 \
    pynndescent==0.5.10 \
    matplotlib==3.7.3 \
    seaborn==0.12.2 \
    h5py==3.9.0 \
    scipy==1.11.2 \
    tqdm==4.66.1 \
    natsort==8.4.0 \
    leidenalg

# Setup pipeline
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ ./scripts/

CMD ["snakemake", "--help"]