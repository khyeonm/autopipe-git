FROM condaforge/miniforge3:latest

# Install snakemake and bash
RUN conda install -y -c bioconda -c conda-forge \
    snakemake-minimal \
    bash \
    && conda clean -afy

# Replace system bash with conda bash (prevents GLIBC mismatch)
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

# --- Motif discovery (MEME), isolated in its OWN conda env so it does not
# --- clash with the workflow engine in the base env. ghostscript (+imagemagick)
# --- let MEME rasterise sequence logos to PNG, not only EPS. MEME is always
# --- called by full path: /opt/conda/envs/meme/bin/meme
RUN conda create -y -n meme -c bioconda -c conda-forge \
    meme=5.5.9 \
    ghostscript \
    imagemagick \
    && conda clean -afy

# Build-time self-test: prove MEME installed, runs SINGLE-PROCESS (no -p, so no
# MPI daemons), and emits a PNG logo. If any of that fails, the image build fails.
RUN set -e; \
    printf '>s1\nGGGGGGTGACGTCAGGGGGG\n>s2\nCCCCCCTGACGTCACCCCCC\n>s3\nAAAAAATGACGTCAAAAAAA\n>s4\nTTTTTTTGACGTCATTTTTT\n>s5\nACACACTGACGTCAACACAC\n>s6\nGTGTGTTGACGTCAGTGTGT\n' > /tmp/meme_selftest.fa; \
    export PATH=/opt/conda/envs/meme/bin:$PATH; \
    /opt/conda/envs/meme/bin/meme /tmp/meme_selftest.fa -dna -oc /tmp/meme_selftest_out \
        -mod zoops -nmotifs 1 -minw 6 -maxw 8 -nostatus; \
    test -n "$(find /tmp/meme_selftest_out -name '*.png' -print -quit)" \
        || { echo 'MEME self-test FAILED: no PNG logo produced (ghostscript issue)'; exit 1; }; \
    echo 'MEME self-test OK: single-process run produced PNG logos'; \
    rm -rf /tmp/meme_selftest.fa /tmp/meme_selftest_out

# Setup pipeline
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/

CMD ["snakemake", "--help"]