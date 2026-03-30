# AptaSelect

AptaSelect identifies high-frequency aptamer candidate sequences from paired-end FASTQ files produced by SELEX experiments.

## Pipeline Stages

1. **Join** — Merges paired-end reads by finding optimal overlap using Hamming distance and Phred quality scores
2. **Selection Filtering** — Extracts sequences flanked by selection adapter patterns
3. **1st Sort Filtering** — Validates extracted sequences contain 1st sort patterns with required spacing (40 bp)
4. **2nd Sort Filtering** — Validates sequences contain 2nd sort patterns with required spacing (20 bp)
5. **Aggregation & Ranking** — Counts unique sequences and reports the top N candidates

## Required Inputs

- `test_1.fq` — Read 1 paired-end FASTQ file
- `test_2.fq` — Read 2 paired-end FASTQ file

## Outputs

- `aptaselect_top_sequences.tsv` — Ranked list of top aptamer candidate sequences with counts
- `aptaselect_stats.json` — Pipeline statistics (reads processed, pass rates per stage)

## Configuration

Edit `config.yaml` to adjust parameters including adapter patterns, overlap settings, mismatch tolerance, and the number of top sequences to report.

## Running

```bash
# Build
docker build -t autopipe-aptaselect .

# Run
docker run --rm \
  -v /path/to/input:/input:ro \
  -v /path/to/output:/output \
  autopipe-aptaselect \
  snakemake --cores 1 -s /pipeline/Snakefile
```