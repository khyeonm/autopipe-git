# AptaSelect

AptaSelect finds high-frequency aptamer candidate sequences in paired-end SELEX FASTQ files. Each read is passed through three sequential filters (constant-region matching, flanking extraction, random-region extraction), and each stage's survivors are ranked by count.

This is the **baseline** version: reads are processed one after another in a single process, and Stage 1's approximate matching is a plain-Python edit-distance search (no multiprocessing, no compiled or external string-matching libraries).

## Insert layout

```
LCR ─ left stem ─ N (random region) ─ right stem ─ RCR
```

The insert length is derived as `len(LCR) + len(left stem) + N + len(right stem) + len(RCR)` and reported in `survival_summary.tsv`.

## Inputs

| key | description |
|---|---|
| `r1`, `r2` | Read 1 / Read 2 FASTQ (`.fastq` or `.fastq.gz`). Each read is handled as its own sequence, not merged into pairs. |

## Processing

1. **Stage 1: matching.** Keep a read if the full LCR and the full RCR each occur anywhere in it within `max_errors` edits (substitution, insertion or deletion each count as one). Their relative order is not checked. The read is tested as given first; if that fails, its reverse complement is tested. The orientation that carries both regions is kept, so kept reads always read LCR…RCR. Reads where neither orientation carries both are discarded.
2. **Stage 2: flanking extraction.** Find `left_marker` and `right_marker` by exact match and keep the sequence between them. Starting from the leftmost left marker, follow the uninterrupted run of A/C/G/T after it and extend to the farthest right marker inside that run. If an ambiguous base (e.g. N) blocks the way, try the next left-marker occurrence. There is no length limit and no stem check at this stage.
3. **Stage 3: random-region extraction.** Find `left_stem` and `right_stem` by exact match, including overlapping positions. Accept the first left/right pair whose gap is exactly `random_region_length`; that gap is N.

Each stage's survivors, plus the raw input (stage 0), are counted per identical sequence and sorted by count, highest first. The sort is stable, so ties stay in first-seen order.

## Outputs (top level of the output folder)

| file | content |
|---|---|
| `stage0_raw.ranked.tsv` | every input read |
| `stage1_matching.ranked.tsv` | oriented reads passing Stage 1 |
| `stage2_flanking.ranked.tsv` | sequences between the markers |
| `stage3_random_region.ranked.tsv` | **main result**: N sequences |
| `survival_summary.tsv` | `total_reads_processed`, `raw_unique_sequences`, `stage1_matching_survivors`, `stage2_flanking_survivors`, `stage3_random_region_survivors`, `derived_insert_length` (no header) |

Each ranked table has a header line `sequence<TAB>count`, followed by one row per unique sequence. Survivor values count every surviving sequence, not just unique ones. The log is written to `logs/aptaselect.log`.

## Configuration (`config.yaml`)

| key | meaning |
|---|---|
| `lcr`, `rcr` | constant-region sequences (Stage 1) |
| `max_errors` | edit-distance tolerance for each constant region |
| `left_marker`, `right_marker` | exact flanking markers (Stage 2) |
| `left_stem`, `right_stem` | exact stem sequences (Stage 3) |
| `random_region_length` | exact N length |
| `output_dir` | defaults to `/output` |
| `progress_every` | progress log interval |

All sequence values must contain only A/C/G/T.

## Run

```bash
docker build -t aptaselect .
docker run --rm -v /path/to/fastqs:/input:ro -v /path/to/results:/output \
  -v $(pwd)/config.yaml:/pipeline/config.yaml aptaselect \
  snakemake --cores 1
```