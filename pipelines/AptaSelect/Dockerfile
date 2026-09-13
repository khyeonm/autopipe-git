FROM condaforge/miniforge3:latest

# AptaSelect is pure-Python sequence processing (no external bioinformatics
# tools required). We install into a DEDICATED conda env ("pipeline") rather
# than the base env: on rootless / overlayfs Docker build hosts, re-linking the
# base env can trigger a "Cannot find a valid extracted directory cache" error,
# which the dedicated-env pattern avoids.
RUN INSTALL="conda create -y -n pipeline"; \
    command -v mamba >/dev/null 2>&1 && INSTALL="mamba create -y -n pipeline"; \
    command -v micromamba >/dev/null 2>&1 && INSTALL="micromamba create -y -n pipeline"; \
    $INSTALL -c bioconda -c conda-forge \
    snakemake-minimal \
    && conda clean -afy

# Put the dedicated env on PATH so `python` and `snakemake` resolve to it.
ENV PATH=/opt/conda/envs/pipeline/bin:$PATH

# Replace system bash with the env's bash (prevents GLIBC mismatch), guarded so
# a missing bash does not leave a dangling /bin/sh.
RUN if [ -x /opt/conda/envs/pipeline/bin/bash ]; then \
        ln -sf /opt/conda/envs/pipeline/bin/bash /usr/bin/bash && \
        ln -sf /opt/conda/envs/pipeline/bin/bash /bin/sh; \
    fi

# Setup pipeline
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts ./scripts

CMD ["snakemake", "--help"]