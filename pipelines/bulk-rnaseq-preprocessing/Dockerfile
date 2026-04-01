FROM condaforge/miniforge3:latest
RUN conda install -y -c bioconda -c conda-forge snakemake-minimal bash openjdk=17.0.10 nextflow=25.04.3 && conda clean -afy
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && ln -sf /opt/conda/bin/bash /bin/sh
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
CMD ["snakemake", "--help"]