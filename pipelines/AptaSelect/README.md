# AptaSelect

Identifies high-frequency aptamer candidate sequences from paired-end FASTQ files produced by SELEX experiments.

## Pipeline Stages

1. **Join** — Overlapping paired-end reads are merged using Hamming-distance scoring. Supports both short-library (`IS_SHORT=true`) and long-library (`IS_SHORT=false`) modes.
2. **Selection Filtering** — Extracts full insert regions bounded by `SEL_READ1` / `SEL_READ2` primer patterns (fuzzy matching, ≤ `MAX_MISMATCHES`).
3. **1st Sort Filtering** — Narrows to intermediate structural regions flanked by `S1_READ1` / `S1_READ2`, enforcing a fixed between-pattern length of `S1_LENGTH` bp.
4. **2nd Sort Filtering** — Extracts the core variable region bounded by `S2_READ1` / `S2_READ2`, enforcing `S2_LENGTH` bp between patterns.
5. **Aggregation & Ranking** — Counts across all parallel chunks are summed and the top `TOP_N` sequences are written to the final output.

## Inputs

| File | Description |
|------|-------------|
| `read1.fastq.gz` | Read 1 FASTQ (plain or gzipped) |
| `read2.fastq.gz` | Read 2 FASTQ (plain or gzipped) |

## Outputs

| File | Description |
|------|-------------|
| `aptaselect_results.tsv` | Tab-separated: rank, sequence, count |

## Configuration (`config.yaml`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_overlap` | 6 | Minimum overlap bp for join |
| `pct_diff` | 8 | Max mismatch % in overlap |
| `is_short` | `true` | Short library mode |
| `sel_read1` | `CCACTTCTCCTTCCATCCTAAAC` | Selection left pattern |
| `sel_read2` | `GAGTAGTTTGGAGGGTTGTCTG` | Selection right pattern |
| `s1_read1` | `TCCTAAAC` | 1st Sort left pattern |
| `s1_read2` | `GAGTAGTT` | 1st Sort right pattern |
| `s1_length` | 40 | Between-pattern length for 1st Sort |
| `s2_read1` | `TCTCTCTCTC` | 2nd Sort left pattern |
| `s2_read2` | `GAGAGAGAGA` | 2nd Sort right pattern |
| `s2_length` | 20 | Between-pattern length for 2nd Sort |
| `max_mismatches` | 1 | Max mismatches per pattern |
| `top_n` | 10 | Top sequences to report |
| `chunk_size` | 10000 | Records per parallel processing chunk |

## Running

```bash
# Build
docker build -t aptaselect:1.0.0 .

# Run
docker run --rm \
  -v /path/to/fastq:/input:ro \
  -v /path/to/output:/output \
  aptaselect:1.0.0 \
  snakemake --cores 4
```

Input FASTQ files must be named `read1.fastq.gz` and `read2.fastq.gz` (or update `r1`/`r2` in `config.yaml`).