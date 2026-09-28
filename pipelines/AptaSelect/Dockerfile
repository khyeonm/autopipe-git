FROM condaforge/miniforge3:latest

# AptaSelect needs Python + Snakemake + edlib (compiled C/C++ edit-distance
# library used for Stage 1 approximate matching).
# Installed into a dedicated env to avoid the base-env re-link bug
# ("Cannot find a valid extracted directory cache") on rootless/overlayfs hosts.
RUN CREATE="conda create -y"; \
    command -v mamba >/dev/null 2>&1 && CREATE="mamba create -y"; \
    $CREATE -n pipeline -c conda-forge -c bioconda \
    python=3.11 \
    snakemake-minimal=8.20.5 \
    python-edlib \
    bash \
    && conda clean -afy

ENV PATH=/opt/conda/envs/pipeline/bin:$PATH

RUN python -c "import edlib; print('edlib', edlib.align('ACGT', 'TTACGTT', mode='HW', task='distance', k=1))"

RUN if [ -x /opt/conda/envs/pipeline/bin/bash ]; then \
        ln -sf /opt/conda/envs/pipeline/bin/bash /usr/bin/bash && \
        ln -sf /opt/conda/envs/pipeline/bin/bash /bin/sh; \
    fi

WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/

CMD ["snakemake", "--help"]