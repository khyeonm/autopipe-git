FROM condaforge/miniforge3:latest

# Dedicated env pinned to Python 3.12 so that edlib installs from a prebuilt wheel
# (the base image ships Python 3.14, for which edlib has no wheel and no compiler is present).
RUN CREATE="conda create -y"; \
    command -v mamba >/dev/null 2>&1 && CREATE="mamba create -y"; \
    $CREATE -n pipeline -c bioconda -c conda-forge \
    python=3.12 \
    snakemake-minimal \
    bash \
    && conda clean -afy

ENV PATH=/opt/conda/envs/pipeline/bin:$PATH

# Replace system bash with conda bash (prevents GLIBC mismatch), guarded
RUN if [ -x /opt/conda/envs/pipeline/bin/bash ]; then \
        ln -sf /opt/conda/envs/pipeline/bin/bash /usr/bin/bash && \
        ln -sf /opt/conda/envs/pipeline/bin/bash /bin/sh; \
    fi

# edlib (Myers bit-vector edit distance, used for Stage 1 LCR/RCR search), prebuilt wheel only
RUN pip install uv && \
    uv pip install --python /opt/conda/envs/pipeline/bin/python --only-binary edlib edlib==1.3.9.post1

# Fail the build early if edlib is not importable
RUN python -c "import edlib, sys; print('edlib ok on', sys.version)"

WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/

CMD ["snakemake", "--help"]