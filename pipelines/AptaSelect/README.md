# AptaSelect

AptaSelect identifies high-frequency aptamer candidate sequences from paired-end
FASTQ files produced by SELEX experiments. It runs three sequential filtering
stages followed by aggregation and ranking. Every pattern sequence, mismatch
tolerance and length is user-configurable, so the same pipeline can be reused for
different library designs by editing only `config.yaml`.

## Insert layout

Each insert is assumed to have the structure:

```
LCR - left stem - random region (N) - right stem - RCR
```

Read 1 and Read 2 are each treated as their own sequence, processed independently
through all three stages, and their survivors are pooled into a single per-stage
result set.

## Stages

- **Stage 0 – Raw.** All raw reads aggregated by identical sequence.
- **Stage 1 – Matching (orient by LCR).** Keep reads containing the LCR pattern
  followed downstream by the RCR pattern, each matched within a configurable
  mismatch tolerance. Both the read and its reverse complement are checked; the
  matching orientation is carried forward.
- **Stage 2 – Flanking extraction (exact match).** Extract the region between a
  short inner-LCR marker and a short inner-RCR marker, both matched exactly.
  Length is not checked at this stage.
- **Stage 3 – N (random region) extraction.** Within the extracted region, find
  the left and right stems by exact match and require the between-region to be
  exactly the target length. Repetitive stems are handled by considering all
  overlapping left/right position pairs and accepting the first that satisfies
  the length constraint.

## Inputs

- `r1` — paired-end Read 1 FASTQ (`.fq`/`.fastq`, optionally gzipped)
- `r2` — paired-end Read 2 FASTQ (`.fq`/`.fastq`, optionally gzipped)

## Outputs (written to `/output`)

- `stage0_raw.ranked.tsv`
- `stage1_matching.ranked.tsv`
- `stage2_flanking.ranked.tsv`
- `stage3_random_region.ranked.tsv` — **main output**: each surviving N sequence
  with its count, highest count first
- `survival_summary.tsv` — read counts surviving each stage, plus
  `total_reads_processed` and `raw_unique_sequences`

Each ranked TSV has a `sequence<TAB>count` header and is sorted by count
descending, with ties broken by first-seen order (stable).

## Configuration (`config.yaml`)

| Key | Description |
|-----|-------------|
| `r1`, `r2` | Paired-end input FASTQ files |
| `lcr`, `rcr` | Left / right constant region sequences |
| `constant_mismatch` | Max mismatches allowed when matching LCR / RCR |
| `left_marker`, `right_marker` | Inner LCR / RCR markers (exact match) |
| `left_stem`, `right_stem` | Left / right stem sequences (exact match) |
| `random_region_length` | Exact length of the random region N |

## Running

```bash
docker build -t aptaselect .
docker run --rm \
    -v /path/to/inputs:/input:ro \
    -v /path/to/outputs:/output \
    aptaselect \
    snakemake --cores 1
```