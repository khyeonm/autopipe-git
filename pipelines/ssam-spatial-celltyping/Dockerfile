FROM condaforge/mambaforge:latest

RUN mamba install -y -c conda-forge -c bioconda \
    snakemake-minimal \
    python=3.10 \
    numpy \
    pandas \
    scipy \
    scikit-learn \
    matplotlib \
    scanpy \
    gxx_linux-64 \
    bash \
    && mamba clean -afy

RUN ln -sf $(which x86_64-conda-linux-gnu-g++) /usr/local/bin/g++ \
    && ln -sf $(which x86_64-conda-linux-gnu-gcc) /usr/local/bin/gcc

RUN pip install uv

RUN uv pip install --system \
    "anndata==0.10.9" \
    "git+https://github.com/pnucolab/ssam.git" \
    spatialdata

# Replace system bash and sh with conda bash (avoids host GLIBC bind-mount issue)
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/