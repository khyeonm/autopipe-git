# AptaSelect

AptaSelect finds high-frequency aptamer candidate sequences in paired-end SELEX FASTQ files. Each read is passed through three sequential filters (constant-region matching, flanking extraction, random-region extraction), and each stage's survivors are ranked by count.

Version 1.0.1 is the **optimized** version. It produces byte-for-byte the same output files as 1.0.0, but:

- reads are processed by several worker processes in parallel, and the results are merged back **in the original read order**, so first-seen order (and therefore the tie order of the stable ranking) is unchanged;
- Stage 1 approximate matching runs in compiled code ([edlib](https://github.com/Martinsos/edlib), infix edit distance) instead of a Python loop, with the same edit-distance tolerance;
- a quick exact-match prefilter skips the approximate match when it cannot succeed;
- progress is logged to `logs/aptaselect.log`.

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
   - *Prefilter:* each constant region is split once into `max_errors + 1` non-overlapping pieces. `k` edits can damage at most `k` pieces, so a match within `k` edits always leaves at least one piece intact. If no piece occurs exactly in the read, the approximate match is skipped (it could not succeed). This never changes which reads are kept.
2. **Stage 2: flanking extraction.** Find `left_marker` and `right_marker` by exact match and keep the sequence between them. Starting from the leftmost left marker, follow the uninterrupted run of A/C/G/T after it and extend to the farthest right marker inside that run. If an ambiguous base (e.g. N) blocks the way, try the next left-marker occurrence. There is no length limit and no stem check at this stage.
3. **Stage 3: random-region extraction.** Find `left_stem` and `right_stem` by exact match, including overlapping positions. Accept the first left/right pair whose gap is exactly `random_region_length`; that gap is N.

Each stage's survivors, plus the raw input (stage 0), are counted per identical sequence and sorted by count, highest first. The sort is stable, so ties stay in first-seen order.

### Parallelism

Reads are grouped into batches of `batch_size` reads and sent to a pool of worker processes, which only run the three stages. All counting happens in the main process, which always takes the oldest outstanding batch first, so every counter sees sequences in exactly the original read order regardless of which worker finishes first.

With `workers: 0` (default) the number of workers is detected automatically from the CPUs actually allocated to the process — the CPU affinity mask, capped by any container (cgroup) CPU quota — not the machine total, and it is also capped by Snakemake's `--cores`. Set `workers` to a positive number to override. Any `batch_size` gives identical results; a value so large that there are fewer batches than workers leaves cores idle.

## Outputs (top level of the output folder)

| file | content |
|---|---|
| `stage0_raw.ranked.tsv` | every input read |
| `stage1_matching.ranked.tsv` | oriented reads passing Stage 1 |
| `stage2_flanking.ranked.tsv` | sequences between the markers |
| `stage3_random_region.ranked.tsv` | **main result**: N sequences |
| `survival_summary.tsv` | `total_reads_processed`, `raw_unique_sequences`, `stage1_matching_survivors`, `stage2_flanking_survivors`, `stage3_random_region_survivors`, `derived_insert_length` (no header) |

Each ranked table has a header line `sequence<TAB>count`, followed by one row per unique sequence. Survivor values count every surviving sequence, not just unique ones.

### Progress log

`logs/aptaselect.log` (separate from the result files) records the settings, worker count and container memory limit, the start and end of each input file, periodic progress lines (reads processed, % of input, speed, ETA, per-stage survivors, unique reads so far, main-process memory), the sorting/writing of each output table, and a final timing breakdown.

## Configuration (`config.yaml`)

| key | meaning |
|---|---|
| `lcr`, `rcr` | constant-region sequences (Stage 1) |
| `max_errors` | edit-distance tolerance for each constant region |
| `left_marker`, `right_marker` | exact flanking markers (Stage 2) |
| `left_stem`, `right_stem` | exact stem sequences (Stage 3) |
| `random_region_length` | exact N length |
| `output_dir` | defaults to `/output` |
| `workers` | worker processes; `0` = auto-detect |
| `batch_size` | reads per work unit sent to a worker (default 2000) |
| `progress_every` | log a progress line every this many reads (0 disables) |
| `progress_seconds` | also log a progress line after this many seconds without one (0 disables) |

All sequence values must contain only A/C/G/T.

## Run

```bash
docker build -t aptaselect .
docker run --rm -v /path/to/fastqs:/input:ro -v /path/to/results:/output \
  -v $(pwd)/config.yaml:/pipeline/config.yaml aptaselect \
  snakemake --cores 16
```

Memory: the main process keeps every unique sequence of every stage in memory; on ~100 million 150-nt reads this peaked at about 37 GB.