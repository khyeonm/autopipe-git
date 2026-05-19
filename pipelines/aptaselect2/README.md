# AptaSelect2

Multi-sample version of AptaSelect. Automatically detects paired-end FASTQ files from the input directory, runs the AptaSelect pipeline independently per sample, and generates a group comparison table on Stage 4 (2nd Sort) results.

Based on [aptaselect (#80)](https://hub.autopipe.org/pipelines/80).

## Pipeline Stages (per sample)

1. **Join** – Reverse-complements Read 2, finds optimal overlap with Read 1, and merges into a single sequence
2. **Selection Filtering** – Extracts the region between configurable selection primer patterns
3. **1st Sort Filtering** – Validates presence of 1st sort primer patterns within the extracted sequence
4. **2nd Sort Filtering** – Validates 2nd sort primer patterns with a required between-length constraint
5. **Aggregation & Ranking** – Deduplicates and ranks sequences by frequency at each stage
6. **Group Comparison** – Aggregates Stage 4 results across samples by group and produces a comparison table with total count, mean count, and fold change (when exactly 2 groups)

## Required Inputs

Place paired-end FASTQ files in the input directory following this naming convention:

- `{sample}_1.fq[.gz]` — Read 1
- `{sample}_2.fq[.gz]` — Read 2

Samples are auto-detected by filename pattern. Multiple samples are supported.

## Outputs

Per sample (`/output/{sample}/`):
- `stage1_joined_ranked.tsv` — Frequency-ranked joined sequences
- `stage2_selection_ranked.tsv` — Frequency-ranked sequences after selection filtering
- `stage3_sort1_ranked.tsv` — Frequency-ranked sequences after 1st sort filtering
- `stage4_sort2_ranked.tsv` — Frequency-ranked sequences after 2nd sort filtering
- `summary.txt` — Statistics for each stage

Cross-sample:
- `group_comparison.tsv` — Group-level frequency comparison on Stage 4 results

## Configuration

Edit `config.yaml` to assign samples to groups and adjust parameters:

### Group assignment

```yaml
groups:
  sampleA: "round3"
  sampleB: "round3"
  sampleC: "round5"
  sampleD: "round5"
```

Samples not listed in `groups` are assigned to "ungrouped".

### Pipeline parameters

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