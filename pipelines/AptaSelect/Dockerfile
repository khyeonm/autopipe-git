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

# ---- motif analysis add-on: MEME in its OWN env ------------------------------
# Kept separate from the Snakemake env so their dependencies cannot clash.
# The pipeline calls it by full path (/opt/conda/envs/meme/bin/meme) and never
# adds this env to the global PATH. ghostscript lets MEME write PNG logos
# (not only EPS). MEME is always run serially: the parallel option (-p) is not used.
RUN CREATE="conda create -y"; \
    command -v mamba >/dev/null 2>&1 && CREATE="mamba create -y"; \
    $CREATE -n meme -c conda-forge -c bioconda \
    meme=5.5.9 \
    ghostscript \
    && conda clean -afy

# Build-time check: MEME runs serially on DNA and writes meme.xml, meme.html and a PNG logo.
RUN printf '>s1\nACGTAGGATCCTTAGCAT\n>s2\nTTGCAGGATCCTTAGCAAC\n>s3\nCCGAGGATCCTTAGCGTA\n>s4\nAATAGGATCCTTAGCTTT\n' > /tmp/memetest.fa \
    && PATH=/opt/conda/envs/meme/bin:$PATH OMP_NUM_THREADS=1 \
       /opt/conda/envs/meme/bin/meme /tmp/memetest.fa -dna -oc /tmp/memetest \
       -mod zoops -nmotifs 1 -minw 6 -maxw 10 -nostatus \
    && ls -l /tmp/memetest \
    && test -s /tmp/memetest/meme.xml \
    && test -s /tmp/memetest/meme.html \
    && test -s /tmp/memetest/logo1.png \
    && /opt/conda/envs/meme/bin/meme -version \
    && rm -rf /tmp/memetest /tmp/memetest.fa

WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/

CMD ["snakemake", "--help"]