# Octopus Variant Calling Pipeline

A Snakemake pipeline for calling SNPs and indels from paired-end whole-genome sequencing data using Octopus, a haplotype-aware Bayesian variant caller.

## Pipeline Steps

1. **FastQC** — Quality control of raw reads
2. **fastp** — Adapter trimming and quality filtering
3. **BWA-MEM2** — Read alignment to reference genome
4. **Octopus** — Haplotype-aware variant calling (SNPs, indels, complex variants)
5. **bcftools filter** — Quality filtering
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
docker build -t octopus-pipeline .
docker run --rm \
  -v /path/to/input:/input \
  -v /path/to/output:/output \
  octopus-pipeline \
  snakemake --cores 4 -s /pipeline/Snakefile
```
