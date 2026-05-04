# AptaSelect

AptaSelect identifies high-frequency aptamer candidate sequences from paired-end FASTQ files produced by SELEX experiments.

## Pipeline Stages

1. **Join** – Reverse-complements Read 2, finds optimal overlap with Read 1, and merges into a single sequence
2. **Selection Filtering** – Extracts the region between configurable selection primer patterns
3. **1st Sort Filtering** – Validates presence of 1st sort primer patterns within the extracted sequence
4. **2nd Sort Filtering** – Validates 2nd sort primer patterns with a required between-length constraint
5. **Aggregation & Ranking** – Deduplicates and ranks sequences by frequency at each stage

## Required Inputs

- `test_1.fq` — Read 1 paired-end FASTQ file
- `test_2.fq` — Read 2 paired-end FASTQ file

Gzipped FASTQ files (`.fq.gz`) are also supported.

## Outputs

- `stage1_joined_ranked.tsv` — Frequency-ranked joined sequences
- `stage2_selection_ranked.tsv` — Frequency-ranked sequences after selection filtering
- `stage3_sort1_ranked.tsv` — Frequency-ranked sequences after 1st sort filtering
- `stage4_sort2_ranked.tsv` — Frequency-ranked sequences after 2nd sort filtering
- `summary.txt` — Statistics for each stage

## How to Run

```bash
# Build the Docker image
docker build -t aptaselect /path/to/aptaselect/

# Run the pipeline
docker run --rm \
  -v /path/to/input:/input:ro \
  -v /path/to/output:/output \
  aptaselect \
  snakemake --cores 4 -s /pipeline/Snakefile
```

## Configuration

Edit `config.yaml` to adjust parameters:

| Parameter | Default | Description |
|---|---|---|
| `long_mode` | `false` | Set `true` for long library inserts |
| `min_overlap` | `6` | Minimum overlap length (bp) |
| `max_mismatch_pct` | `0.08` | Maximum mismatch rate for joining |
| `max_mismatches` | `1` | Pattern matching mismatch tolerance |
| `sel_left` / `sel_right` | see config | Selection primer patterns |
| `sort1_left` / `sort1_right` | see config | 1st Sort primer patterns |
| `sort2_left` / `sort2_right` | see config | 2nd Sort primer patterns |
| `sort2_between_length` | `20` | Required between-length for 2nd sort |
| `chunk_size` | `10000` | Processing chunk size |