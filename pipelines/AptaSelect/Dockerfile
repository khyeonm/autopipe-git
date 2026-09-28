FROM condaforge/miniforge3:latest

# Baseline AptaSelect needs only Python's standard library + Snakemake
# (no compiled / external string-matching libraries by design).
# Installed into a dedicated env to avoid the base-env re-link bug
# ("Cannot find a valid extracted directory cache") on rootless/overlayfs hosts.
RUN CREATE="conda create -y"; \
    command -v mamba >/dev/null 2>&1 && CREATE="mamba create -y"; \
    $CREATE -n pipeline -c bioconda -c conda-forge \
    python=3.11 \
    snakemake-minimal=8.20.5 \
    bash \
    && conda clean -afy

ENV PATH=/opt/conda/envs/pipeline/bin:$PATH

RUN if [ -x /opt/conda/envs/pipeline/bin/bash ]; then \
        ln -sf /opt/conda/envs/pipeline/bin/bash /usr/bin/bash && \
        ln -sf /opt/conda/envs/pipeline/bin/bash /bin/sh; \
    fi

WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/

CMD ["snakemake", "--help"]