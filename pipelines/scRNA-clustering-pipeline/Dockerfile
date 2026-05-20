# =============================================================================
# Dockerfile for scRNA-clustering-pipeline
# Scanpy-based single-cell RNA-seq analysis
# =============================================================================
FROM condaforge/miniforge3:latest

# Install all Python packages via uv
RUN pip install uv && \
    uv pip install --system \
    snakemake \
    scanpy \
    anndata \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    umap-learn \
    igraph \
    leidenalg \
    h5py \
    scikit-misc && \
    conda clean -afy

# Copy pipeline files
WORKDIR /pipeline
COPY Snakefile config.yaml /pipeline/
COPY scripts/ /pipeline/scripts/

CMD ["bash"]