# AptaSelect

AptaSelect identifies high-frequency aptamer candidate sequences from paired-end FASTQ files produced by SELEX experiments.

## Pipeline Stages

1. **Join** — Reverse-complement Read 2, find optimal overlap, merge reads
2. **Selection Filtering** — Extract region between configurable left/right adapter patterns
3. **1st Sort Filtering** — Validate presence of inner pattern pair
4. **2nd Sort Filtering** — Validate presence of variable-region flanking patterns with exact between-length
5. **Aggregation & Ranking** — Deduplicate and rank sequences by frequency at each stage

## Required Inputs

- `test_1.fq` — Read 1 paired-end FASTQ (plain or gzipped)
- `test_2.fq` — Read 2 paired-end FASTQ (plain or gzipped)

## Outputs

- `stage1_joined.tsv` — Frequency-ranked joined sequences
- `stage2_selection.tsv` — Frequency-ranked sequences after selection filtering
- `stage3_sort1.tsv` — Frequency-ranked sequences after 1st sort filtering
- `stage4_sort2.tsv` — Frequency-ranked sequences after 2nd sort filtering
- `pipeline_stats.json` — Summary statistics (counts at each stage)

## Configuration

All parameters are adjustable in `config.yaml`:

- Library mode (short/long), overlap thresholds
- Pattern sequences and mismatch tolerance for each filtering stage
- Between-length constraint for 2nd Sort
- Processing chunk size

## Running

```bash
docker build -t autopipe-aptaselect .
docker run --rm \
  -v /path/to/input:/input:ro \
  -v /path/to/output:/output \
  autopipe-aptaselect \
  snakemake --cores 1 -s /pipeline/Snakefile
```