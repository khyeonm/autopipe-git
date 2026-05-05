FROM condaforge/miniforge3:latest

RUN conda install -y -c conda-forge -c bioconda \
    python=3.11 \
    snakemake-minimal \
    bash \
    scanpy \
    squidpy \
    anndata \
    pandas \
    numpy \
    scipy \
    matplotlib \
    seaborn \
    openpyxl \
    && conda clean -afy

RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/

CMD ["snakemake", "--help"]