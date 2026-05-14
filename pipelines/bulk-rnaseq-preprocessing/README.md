# Bulk RNA-seq Preprocessing Pipeline

Preprocesses bulk RNA-seq data using the [nf-core/rnaseq](https://nf-co.re/rnaseq) pipeline (v3.22.2) with the **STAR + Salmon** workflow. Takes raw FASTQ files as input and produces aligned BAMs, gene/transcript-level quantification, and a comprehensive MultiQC quality report.

## Required Inputs

All input files should be placed in the input directory mounted at `/input`.

| File | Description |
|---|---|
| `samplesheet.csv` | nf-core samplesheet with columns: `sample`, `fastq_1`, `fastq_2`, `strandedness`. Paths in the samplesheet must use the **host** paths (the pipeline maps them automatically). |
| `reference/genome.fa` | Reference genome FASTA (required when `genome: null`) |
| `reference/genes.gtf` | Gene annotation GTF (required when `genome: null`) |
| `*.fastq.gz` | Raw paired-end or single-end FASTQ files |

### Samplesheet Example

```csv
sample,fastq_1,fastq_2,strandedness
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,auto
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,auto
```

> **Note:** The `fastq_1` and `fastq_2` paths in the samplesheet must point to the **host** file system paths (not `/input/...`) because Nextflow manages its own container mounts.

## Expected Outputs

Written to the output directory mounted at `/output`:

- `nfcore_rnaseq/multiqc/star_salmon/multiqc_report.html` — MultiQC report
- `nfcore_rnaseq/star_salmon/` — STAR-aligned BAM files
- `nfcore_rnaseq/star_salmon/salmon.merged.gene_counts.tsv` — Merged gene-level counts
- `nfcore_rnaseq/star_salmon/salmon.merged.gene_tpm.tsv` — Merged gene-level TPM
- `nfcore_rnaseq/trimgalore/` — Trimmed FASTQ files
- `nfcore_rnaseq/fastqc/` — FastQC reports

## Configuration

Edit `config.yaml` to adjust parameters:

- `nfcore_rnaseq_version` — nf-core/rnaseq release (default: `3.22.2`)
- `aligner` — Alignment strategy (default: `star_salmon`)
- `genome` — iGenomes key or `null` to use custom fasta/gtf (default: `null`)
- `fasta` / `gtf` — Reference files relative to `/input` (required when genome is null)
- `max_cpus` / `max_memory` / `max_time` — Resource limits
- `extra_args` — Additional arguments passed to `nextflow run`

## How to Run

```bash
# Build the Docker image
docker build -t bulk-rnaseq-preprocessing .

# Run the pipeline (Docker socket required for nf-core)
docker run --rm \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /path/to/input:/input:ro \
    -v /path/to/output:/output \
    -e HOST_INPUT_DIR=/path/to/input \
    -e HOST_OUTPUT_DIR=/path/to/output \
    -e HOST_PIPELINE_DIR=/path/to/pipeline \
    bulk-rnaseq-preprocessing \
    snakemake --cores all
```

> **Important:** The Docker socket mount (`-v /var/run/docker.sock:...`) is required because nf-core/rnaseq pulls and runs its own tool containers via Nextflow's Docker profile.