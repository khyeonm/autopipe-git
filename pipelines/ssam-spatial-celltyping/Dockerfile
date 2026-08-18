FROM condaforge/mambaforge:latest

# Step 1: Install packages (mamba clean은 mamba를 망가뜨리므로 conda clean 사용)
RUN mamba install -y -c conda-forge -c bioconda \
    snakemake-minimal \
    python=3.10 \
    numpy \
    pandas \
    scipy \
    scikit-learn \
    matplotlib \
    scanpy \
    leidenalg \
    python-igraph \
    gxx_linux-64 \
    bash \
    && /opt/conda/bin/conda clean -afy

# Step 2: Setup C++ compiler symlinks for ssam build
RUN ln -sf $(which x86_64-conda-linux-gnu-g++) /usr/local/bin/g++ 2>/dev/null || true \
    && ln -sf $(which x86_64-conda-linux-gnu-gcc) /usr/local/bin/gcc 2>/dev/null || true

# Step 3: Install ssam from pnucolab GitHub
RUN pip install uv \
    && uv pip install --system \
        "anndata==0.10.9" \
        "git+https://github.com/pnucolab/ssam.git"

# Step 4: Use conda bash
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash \
    && ln -sf /opt/conda/bin/bash /bin/sh

WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/

# The container runs as a non-root user, so snakemake/matplotlib/numba cannot
# create caches under the default HOME ("/"). Point them at /tmp and make the
# working dir writable so snakemake can create .snakemake/ there.
ENV HOME=/tmp \
    XDG_CACHE_HOME=/tmp/.cache \
    MPLCONFIGDIR=/tmp/.config/matplotlib \
    NUMBA_CACHE_DIR=/tmp/.cache/numba
RUN mkdir -p /tmp/.cache /tmp/.config/matplotlib /tmp/.cache/numba \
    && chmod -R a+rwX /pipeline /tmp/.cache /tmp/.config
