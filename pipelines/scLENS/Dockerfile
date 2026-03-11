FROM condaforge/mambaforge:latest

# Install snakemake and system deps (NO conda julia - use official julia binary instead)
RUN mamba install -y -c conda-forge -c bioconda \
    snakemake-minimal \
    bash \
    curl \
    && mamba clean -afy

# Replace system bash with conda bash (prevents GLIBC mismatch)
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

# Install official Julia 1.11 binary (1.11+ supports 'public' keyword used by UMAP.jl)
RUN curl -fsSL -o /tmp/julia.tar.gz \
    https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-1.11.4-linux-x86_64.tar.gz && \
    tar -xzf /tmp/julia.tar.gz -C /usr/local --strip-components=1 && \
    rm /tmp/julia.tar.gz

# Install uv for Python packages (anndata for h5ad output)
RUN pip install uv
RUN uv pip install --system anndata

WORKDIR /pipeline

# Copy pipeline files
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/

# Install Julia packages using official binary (no libcurl conflict)
RUN julia --project=/pipeline -e " \
    import Pkg; \
    Pkg.add(url=\"https://github.com/Mathbiomed/scLENS\"); \
    Pkg.add([\"CSV\", \"DataFrames\", \"CairoMakie\", \"CUDA\", \"JLD2\", \"HDF5\"]); \
    Pkg.instantiate(); \
    Pkg.precompile()"