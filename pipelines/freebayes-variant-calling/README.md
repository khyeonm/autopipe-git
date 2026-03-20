# Freebayes Variant Calling Pipeline

A Snakemake pipeline for calling SNPs and indels from paired-end whole-genome sequencing data using freebayes.

## Pipeline Steps

1. **FastQC** — Quality control of raw reads
2. **fastp** — Adapter trimming and quality filtering
3. **BWA-MEM2** — Read alignment to reference genome
4. **freebayes** — Bayesian variant calling (SNPs, indels, MNPs)
5. **bcftools filter** — Quality and depth filtering
6. **bcftools stats** — Variant summary statistics

## Input

- Paired-end FASTQ files (`sample_R1.fastq.gz`, `sample_R2.fastq.gz`)
- Reference genome FASTA (`reference.fa`)

## Output

- `sample.filtered.vcf` — Filtered variant calls
- `sample.stats.txt` — Variant statistics
- `fastqc_R1.html`, `fastqc_R2.html` — QC reports
- `fastp.html` — Trimming report

## Quick Start

```bash
# Build Docker image
docker build -t freebayes-pipeline .

# Run pipeline
docker run --rm \
  -v /path/to/input:/input \
  -v /path/to/output:/output \
  freebayes-pipeline \
  snakemake --cores 4 -s /pipeline/Snakefile
```

## Configuration

Edit `config.yaml` to set:
- `sample_r1` / `sample_r2`: Input FASTQ paths
- `reference`: Reference genome path
- `sample_name`: Sample identifier
- `threads`: Number of CPU threads
- `min_qual`: Minimum variant quality (default: 20)
- `min_depth`: Minimum read depth (default: 5)
