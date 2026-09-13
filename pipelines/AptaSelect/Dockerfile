FROM condaforge/miniforge3:latest

# AptaSelect runs on the Python standard library only (no external aligners or
# bioinformatics tools), so we just need snakemake + a Python interpreter.
#
# Install into a DEDICATED env (conda create) rather than the base env. This
# avoids the "Cannot find a valid extracted directory cache" linking error that
# base-env installs can trigger on overlayfs/rootless Docker build hosts.
RUN INSTALL="conda create -y -n pipeline"; \
    command -v mamba >/dev/null 2>&1 && INSTALL="mamba create -y -n pipeline"; \
    command -v micromamba >/dev/null 2>&1 && INSTALL="micromamba create -y -n pipeline"; \
    $INSTALL -c bioconda -c conda-forge \
    snakemake-minimal \
    python=3.11 \
    bash \
    && conda clean -afy

# Put the dedicated env on PATH so snakemake/python resolve to it.
ENV PATH=/opt/conda/envs/pipeline/bin:$PATH

# Replace system bash with the env's conda bash (prevents GLIBC mismatch),
# guarded so a missing bash does not leave a dangling /bin/sh.
RUN if [ -x /opt/conda/envs/pipeline/bin/bash ]; then \
        ln -sf /opt/conda/envs/pipeline/bin/bash /usr/bin/bash && \
        ln -sf /opt/conda/envs/pipeline/bin/bash /bin/sh; \
    fi

# Setup pipeline
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ ./scripts/

CMD ["snakemake", "--help"]