# AptaSelect

Identifies high-frequency aptamer candidate sequences from paired-end FASTQ files
produced by SELEX experiments. It runs three sequential filtering stages followed
by aggregation and ranking.

Each insert has the layout: **LCR – left stem – N (random region) – right stem – RCR**.
Stem lengths follow from their sequences and the random-region length is
configurable, so the overall insert length is derived automatically and never
hard-coded. Read 1 and Read 2 are each treated as their own sequence, processed
independently through all three stages, and their survivors are pooled into a
single per-stage result set.

## Stages

1. **Matching (orient by LCR).** Keep reads containing the LCR followed downstream
   by the RCR, each within a configurable mismatch tolerance. Both the read and its
   reverse complement are checked; the orientation where the LCR…RCR pair is found
   is carried forward.
2. **Flanking extraction (exact match).** Find a short marker at the inner end of
   the LCR and a short marker at the inner start of the RCR by exact match, and
   extract the sequence between them. Length is not checked here.
3. **N (random region) extraction.** Find the left and right stems by exact match
   and require the region between them to be exactly the target length. All
   overlapping stem position pairs are considered; the first pair satisfying the
   length constraint is accepted.

## Aggregation and ranking

For each stage, identical survivors are counted and sorted by descending count,
with ties broken by first-seen order (stable sort, no secondary alphabetical key).
Each stage writes its own ranked TSV, and a survival summary reports how many
reads pass each stage. The main output is the final-stage ranked TSV.

## Inputs

- `r1` — paired-end R1 FASTQ (`.fq`/`.fastq`, optionally gzipped)
- `r2` — paired-end R2 FASTQ (`.fq`/`.fastq`, optionally gzipped)

## Outputs (written to `/output`)

- `stage1_matching.ranked.tsv` — ranked oriented reads passing Stage 1
- `stage2_flanking.ranked.tsv` — ranked extracted regions passing Stage 2
- `stage3_random_region.ranked.tsv` — ranked N sequences passing Stage 3 (main output)
- `survival_summary.tsv` — read counts surviving each stage + derived insert length

## Configuration (`config.yaml`)

| Key | Meaning |
|-----|---------|
| `r1`, `r2` | Input FASTQ paths |
| `lcr`, `rcr` | Left / right constant region sequences (Stage 1) |
| `lcr_rcr_mismatch` | Max substitutions allowed for LCR/RCR matching (Stage 1) |
| `left_marker`, `right_marker` | Inner LCR / RCR markers, exact match (Stage 2) |
| `left_stem`, `right_stem` | Stem sequences, exact match (Stage 3) |
| `random_region_length` | Exact target length of the random region N (Stage 3) |

## How to run

```bash
docker build -t autopipe-aptaselect .
docker run --rm \
    -v /path/to/inputs:/input:ro \
    -v /path/to/outputs:/output \
    autopipe-aptaselect \
    snakemake --cores 1
```