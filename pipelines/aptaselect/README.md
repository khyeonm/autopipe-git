# AptaSelect

AptaSelect identifies high-frequency aptamer candidate sequences from paired-end FASTQ files produced by SELEX experiments.

## Pipeline Stages

1. **Join** — Reverse-complements Read 2, finds the best overlap with Read 1, and merges them into a single sequence using quality-aware consensus.
2. **Selection Filtering** — Extracts the region flanked by the selection left/right patterns.
3. **1st Sort Filtering** — Validates that the extracted sequence contains the 1st sort patterns with the correct between-length.
4. **2nd Sort Filtering** — Validates the 2nd sort patterns with the correct between-length.
5. **Aggregation & Ranking** — Counts unique sequences and reports the top N candidates.

## Required Inputs

- `test_1.fq` — Read 1 paired-end FASTQ file (plain or gzipped)
- `test_2.fq` — Read 2 paired-end FASTQ file (plain or gzipped)

## Outputs

- `aptaselect_top_candidates.tsv` — Ranked list of top candidate sequences with counts
- `aptaselect_stats.json` — Statistics for each pipeline stage

## How to Run

```bash
# Build the Docker image
docker build -t autopipe-aptaselect /home/khyeonmin/projects/auto_test/pipelines/aptaselect

# Run the pipeline
docker run --rm \
  -v /path/to/input:/input:ro \
  -v /path/to/output:/output \
  autopipe-aptaselect \
  snakemake --cores 1
```

## Configuration

All parameters in `config.yaml` are adjustable:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_overlap` | 6 | Minimum overlap length (bp) for joining |
| `max_mismatch_rate` | 0.08 | Maximum mismatch rate in overlap |
| `short_mode` | true | Short library mode (true) or long library mode (false) |
| `pattern_mismatches` | 1 | Max mismatches for pattern matching |
| `selection_left_pattern` | CCACTTCTCCTTCCATCCTAAAC | Selection left flanking pattern |
| `selection_right_pattern` | GAGTAGTTTGGAGGGTTGTCTG | Selection right flanking pattern |
| `sort1_between_length` | 40 | Required between-length for 1st sort |
| `sort2_between_length` | 20 | Required between-length for 2nd sort |
| `top_n` | 10 | Number of top candidates to report |