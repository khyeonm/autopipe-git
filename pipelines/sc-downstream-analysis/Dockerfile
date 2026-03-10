FROM condaforge/mambaforge:latest

RUN conda install -y -c bioconda -c conda-forge \
    snakemake-minimal \
    python=3.10 \
    && conda clean -afy

RUN pip install \
    "scanpy>=1.9" \
    leidenalg \
    python-igraph \
    "umap-learn>=0.5" \
    matplotlib \
    pandas \
    numpy \
    anndata

WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/

CMD ["snakemake", "--cores", "4"]