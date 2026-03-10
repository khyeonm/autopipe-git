FROM condaforge/mambaforge:latest

RUN mamba install -y -c bioconda -c conda-forge \
    snakemake-minimal \
    gatk4=4.5.0.0 \
    samtools=1.20 \
    && mamba clean -afy

WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .

CMD ["snakemake", "--help"]