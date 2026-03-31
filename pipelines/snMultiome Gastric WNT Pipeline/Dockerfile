# =============================================================================
#  Dockerfile — 10X snMultiome Gastric WNT Pipeline
#  Based on the user-provided Dockerfile (v10)
# =============================================================================
FROM condaforge/miniforge3:latest

# System dependencies
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libxml2-dev \
        libpng-dev libbz2-dev liblzma-dev libreadline-dev \
        libx11-dev libxt-dev zlib1g-dev \
        libgl1 libegl1 \
        libhdf5-dev \
        libharfbuzz-dev libfribidi-dev libfreetype6-dev \
        libjpeg-dev libcairo2-dev libtiff5-dev libicu-dev \
        build-essential gfortran procps pkg-config \
    && apt-get clean

# Stage 1: Create conda env with core R + bioconductor foundations
RUN conda create -y -n pipeline -c conda-forge -c bioconda \
        r-base=4.3 python=3.10 snakemake bash \
        r-curl r-xml2 r-httr r-rcurl r-openssl r-png r-data.table r-xml r-httpuv \
        r-igraph r-leiden r-hdf5r r-matrix r-systemfonts r-textshaping \
        zlib \
        bioconductor-rhtslib bioconductor-rsamtools bioconductor-genomicranges \
        bioconductor-biostrings bioconductor-genomeinfodb bioconductor-rtracklayer \
        bioconductor-genomicalignments bioconductor-bsgenome bioconductor-genomicfeatures \
        bioconductor-annotationdbi bioconductor-summarizedexperiment \
        bioconductor-singlecellexperiment bioconductor-xvector \
        bioconductor-cner bioconductor-tfbstools \
    && conda clean --all --yes

# Stage 2: Install remaining R packages into the existing env
RUN conda install -y -n pipeline -c conda-forge -c bioconda \
        r-seurat r-seuratobject r-signac \
        r-ggforce r-tidygraph r-graphlayouts r-ggraph r-clustree \
        r-future r-future.apply \
        r-dplyr r-tidyr r-purrr r-tibble r-ggplot2 \
        r-scales r-ggrepel r-gridextra r-cowplot \
        r-rcolorbrewer r-viridis r-colorspace r-paletteer \
        r-yaml r-hash r-readxl r-ggnewscale r-dendextend \
        r-rcpp r-rcpparmadillo r-rcppeigen \
        r-shiny r-plotly r-dt r-miniui \
        r-patchwork \
    && conda clean --all --yes

ENV PATH="/opt/conda/envs/pipeline/bin:$PATH"
ENV CFLAGS="-I/opt/conda/envs/pipeline/include"
ENV LDFLAGS="-L/opt/conda/envs/pipeline/lib -Wl,-rpath,/opt/conda/envs/pipeline/lib"
ENV CPPFLAGS="-I/opt/conda/envs/pipeline/include"
ENV PKG_CONFIG_PATH="/opt/conda/envs/pipeline/lib/pkgconfig"
SHELL ["conda", "run", "--no-capture-output", "-n", "pipeline", "/bin/bash", "-c"]

# MACS2
RUN pip install macs2

# Bioconductor high-level packages
RUN Rscript -e " \
    cat('Installing Bioconductor packages...\n'); \
    options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/2024-01-01')); \
    install.packages('BiocManager'); \
    BiocManager::install(c( \
        'biovizBase', \
        'EnsDb.Mmusculus.v79', \
        'BSgenome.Mmusculus.UCSC.mm10', \
        'dittoSeq', \
        'motifmatchr', \
        'JASPAR2020', \
        'chromVAR', \
        'msigdbr', \
        'gprofiler2', \
        'MAST' \
    ), ask = FALSE, update = FALSE); \
    if (!require('Seurat', quietly=TRUE)) stop('Seurat FAILED'); \
    if (!require('Signac', quietly=TRUE)) stop('Signac FAILED'); \
    if (!require('MAST', quietly=TRUE)) stop('MAST FAILED'); \
    cat('All packages OK\n')"

# Copy pipeline files
RUN mkdir -p /pipeline/scripts
COPY Snakefile /pipeline/Snakefile
COPY config.yaml /pipeline/config.yaml
COPY scripts/01_load_and_qc.R        /pipeline/scripts/01_load_and_qc.R
COPY scripts/02_macs2_peaks.R        /pipeline/scripts/02_macs2_peaks.R
COPY scripts/03_normalize.R          /pipeline/scripts/03_normalize.R
COPY scripts/04_integrate_cluster.R  /pipeline/scripts/04_integrate_cluster.R
COPY scripts/05_visualization.R      /pipeline/scripts/05_visualization.R
COPY scripts/06_differential.R       /pipeline/scripts/06_differential.R
COPY scripts/07_motif_enrichment.R   /pipeline/scripts/07_motif_enrichment.R

WORKDIR /pipeline