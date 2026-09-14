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

# ---------------------------------------------------------------------------
# APPENDED: MEME Suite for the motif-analysis step.
#
# MEME is installed into its OWN dedicated conda env ("meme"), separate from the
# "pipeline" (snakemake) env, so its dependencies never clash with the workflow
# engine. The motif step calls the MEME binary by its FULL PATH
# (/opt/meme/bin/meme -> symlinked to the env binary) and runs it as the plain
# SERIAL binary (no -p / MPI, which fails in the container).
#
# ghostscript is included so MEME can render sequence logos as PNG (not only
# EPS). It is installed into the same "meme" env alongside MEME.
# ---------------------------------------------------------------------------
RUN INSTALL="conda create -y -n meme"; \
    command -v mamba >/dev/null 2>&1 && INSTALL="mamba create -y -n meme"; \
    command -v micromamba >/dev/null 2>&1 && INSTALL="micromamba create -y -n meme"; \
    $INSTALL -c bioconda -c conda-forge \
    meme \
    ghostscript \
    && conda clean -afy

# Expose the MEME binaries at a stable full path used by the motif step.
# The step invokes /opt/meme/bin/meme explicitly (never relying on PATH order,
# and never adding the meme env to the pipeline env's PATH so the two stay
# isolated).
RUN mkdir -p /opt/meme/bin && \
    for b in meme meme-chip dreme fimo tomtom ceqlogo; do \
        if [ -x /opt/conda/envs/meme/bin/$b ]; then \
            ln -sf /opt/conda/envs/meme/bin/$b /opt/meme/bin/$b; \
        fi; \
    done && \
    ln -sf /opt/conda/envs/meme/bin/gs /opt/meme/bin/gs 2>/dev/null || true

# Verify the serial MEME binary is present and runnable at build time.
RUN /opt/meme/bin/meme -version

# Setup pipeline
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts ./scripts

CMD ["snakemake", "--help"]