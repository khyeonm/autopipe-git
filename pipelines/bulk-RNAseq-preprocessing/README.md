# Bulk RNA-seq Preprocessing Pipeline (nf-core/rnaseq)

Preprocesses bulk RNA-seq paired-end FASTQ files using nf-core/rnaseq v3.22.2 with the STAR-Salmon workflow. iGenomes is disabled by default — users provide their own reference FASTA and GTF files.

## What it does

1. **Quality Control** — FastQC on raw reads
2. **Adapter Trimming** — Trim Galore removes adapters and low-quality bases
3. **Alignment** — STAR aligns reads to the provided reference genome
4. **Quantification** — Salmon quantifies transcript/gene-level expression
5. **Post-alignment QC** — RSeQC, dupRadar, Preseq, and more
6. **Summary Report** — MultiQC aggregates all QC metrics into a single HTML report

## Required Inputs

Place the following in your input directory:

| File | Description |
|------|-------------|
| `samplesheet.csv` | nf-core samplesheet with columns: `sample,fastq_1,fastq_2,strandedness` |
| `reads/*.fastq.gz` | Raw paired-end FASTQ files referenced in the samplesheet |
| `reference/genome.fa` | Reference genome FASTA (required when `genome: null`) |
| `reference/genes.gtf` | Gene annotation GTF (required when `genome: null`) |

### Samplesheet example

```csv
sample,fastq_1,fastq_2,strandedness
CONTROL_1,/input/reads/CTRL1_R1.fastq.gz,/input/reads/CTRL1_R2.fastq.gz,auto
CONTROL_2,/input/reads/CTRL2_R1.fastq.gz,/input/reads/CTRL2_R2.fastq.gz,auto
TREATED_1,/input/reads/TRT1_R1.fastq.gz,/input/reads/TRT1_R2.fastq.gz,auto
```

> **Note:** FASTQ paths in the samplesheet must use the `/input/` prefix.
> Set strandedness to `auto` to let Salmon infer it automatically.

## Expected Outputs

| Output | Path |
|--------|------|
| MultiQC report | `nfcore/multiqc/star_salmon/multiqc_report.html` |
| Salmon gene counts | `nfcore/star_salmon/salmon.merged.gene_counts.tsv` |
| Salmon transcript counts | `nfcore/star_salmon/salmon.merged.transcript_counts.tsv` |
| STAR BAM files | `nfcore/star_salmon/<sample>/` |
| FastQC reports | `nfcore/fastqc/` |
| Trim Galore reports | `nfcore/trimgalore/` |

## Configuration

Edit `config.yaml` to adjust:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nfcore_version` | `3.22.2` | nf-core/rnaseq release version |
| `genome` | `null` | Set to iGenomes name to use iGenomes (disabled by default) |
| `aligner` | `star_salmon` | Aligner workflow |
| `fasta` | `null` | Path to reference FASTA (relative to input dir) |
| `gtf` | `null` | Path to gene annotation GTF (relative to input dir) |
| `star_index` | `null` | Pre-built STAR index directory (optional) |
| `max_cpus` | `8` | Maximum CPU cores |
| `max_memory` | `32.GB` | Maximum RAM |
| `extra_params` | `""` | Additional nf-core/rnaseq CLI flags |

## Running

```bash
# Build the Docker image
docker build -t autopipe-bulk-rnaseq-preprocessing .

# Run the pipeline (Docker socket required for nf-core)
docker run --rm \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /path/to/input:/input:ro \
    -v /path/to/output:/output \
    autopipe-bulk-rnaseq-preprocessing \
    snakemake --cores 8
```

> **Important:** This pipeline requires Docker socket access (`-v /var/run/docker.sock:/var/run/docker.sock`) because nf-core/rnaseq pulls and runs its own tool containers internally.