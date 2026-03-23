# AptaSelect

Identifies high-frequency aptamer candidate sequences from paired-end FASTQ files produced by SELEX experiments.

## Pipeline Stages

1. **Join** — Overlap-based read joining (short/long library modes)
2. **Selection Filtering** — Extract insert regions flanked by primer sequences
3. **1st Sort Filtering** — Narrow to intermediate structural region (40 bp between-length)
4. **2nd Sort Filtering** — Extract the core 20 bp variable region
5. **Aggregation & Ranking** — Parallel counting, deduplication, and top-N output

## Inputs

| Parameter | Description |
|-----------|-------------|
| `read1` | Read 1 FASTQ file (plain or `.gz`) |
| `read2` | Read 2 FASTQ file (plain or `.gz`) |

## Outputs

| File | Description |
|------|-------------|
| `sort2_top.tsv` | **Final aptamer candidates** ranked by frequency |
| `sort1_top.tsv` | Top sequences after 1st Sort |
| `selection_top.tsv` | Top sequences after Selection |
| `joined_top.tsv` | Top sequences after Join |
| `summary.txt` | Per-stage read counts |

## Configuration (`config.yaml`)

| Key | Default | Description |
|-----|---------|-------------|
| `min_overlap` | 6 | Minimum overlap bp for Join |
| `pct_diff` | 8.0 | Max mismatch % in overlap |
| `is_short` | true | Short library mode |
| `sel_read1/2` | — | Selection primer sequences |
| `s1_read1/2` | — | 1st Sort primer sequences |
| `s1_length` | 40 | Required between-length for 1st Sort |
| `s2_read1/2` | — | 2nd Sort primer sequences |
| `s2_length` | 20 | Required between-length for 2nd Sort |
| `max_mismatches` | 1 | Max mismatches in pattern matching |
| `top_n` | 10 | Number of top candidates output |
| `chunk_size` | 10000 | Records per parallel chunk |
| `threads` | 4 | Worker threads |

## Run

```bash
docker build -t aptaselect .

docker run --rm \
  -v /path/to/data:/input:ro \
  -v /path/to/output:/output \
  aptaselect \
  snakemake --cores 4 --snakefile /pipeline/Snakefile
```