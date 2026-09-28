# AptaSelect

AptaSelect finds high-frequency aptamer candidate sequences in paired-end SELEX FASTQ files. Each read is passed through three sequential filters (constant-region matching, flanking extraction, random-region extraction), and each stage's survivors are ranked by count. An add-on step writes the top-count random-region cores to FASTA and can run MEME on them to look for a shared motif.

Version 1.0.1 is the **optimized** version. It produces byte-for-byte the same output files as 1.0.0, but:

- reads are processed by several worker processes in parallel, and the results are merged back **in the original read order**, so first-seen order (and therefore the tie order of the stable ranking) is unchanged;
- Stage 1 approximate matching runs in compiled code ([edlib](https://github.com/Martinsos/edlib), infix edit distance) instead of a Python loop, with the same edit-distance tolerance;
- a quick exact-match prefilter skips the approximate match when it cannot succeed;
- progress is logged to `logs/aptaselect.log`.

Version 1.0.2 appends the motif-analysis step (below). The AptaSelect step and its output files are unchanged.

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

### Motif analysis (add-on)

4. **Top-core selection (always runs).** The top-count rows are taken straight from `stage3_random_region.ranked.tsv`, which already holds only the variable core N (primers, constant regions and stems removed), and written unchanged to `meme_input_top_cores.fasta` with headers `>rank<R>_count<C>`. How many: `meme_top_percent` % of the unique sequences in the table (rounded up, at least 1) when set, otherwise the first `meme_top_n` rows. The cut follows the table order; if it falls inside a group of equal counts, a warning is written to `logs/select_top_cores.log`. The FASTA is a regular output, so it can be uploaded to the MEME website even when MEME is not run here.
5. **MEME (off by default, `run_meme: false`).** When turned on, MEME 5.5.9 runs on that FASTA with the DNA alphabet (`-dna`), in plain serial mode (no `-p`, which would launch MPI), writing to `meme_out/`. Only the variable core is given to MEME, never the constant or primer regions, so they cannot bias the motif search.

**Two-step use:** run once with `run_meme: false` to get the count tables and the FASTA; check them; then set `run_meme: true` and run again with the same run name / output folder. Snakemake sees the existing results are up to date and runs only MEME — the count-sorting results are kept, not recomputed or deleted. If you change `meme_top_n` / `meme_top_percent` or any MEME option between runs, only the FASTA and/or MEME steps re-run (their settings are recorded in `logs/.select_top_cores.settings` and `logs/.meme.settings`); the count-sorting step is never re-run by these settings. Re-running MEME rewrites `meme_out/`.

**Interpreting motifs:** a MEME motif is only a candidate. Before trusting it, check the full sequence separately for structure (G-quadruplex / stem-loop) and folding energy, to rule out motifs that are just SELEX artefacts (e.g. sequences complementary to the constant regions). That structural check is not part of this pipeline.

## Outputs (top level of the output folder)

| file | content |
|---|---|
| `stage0_raw.ranked.tsv` | every input read |
| `stage1_matching.ranked.tsv` | oriented reads passing Stage 1 |
| `stage2_flanking.ranked.tsv` | sequences between the markers |
| `stage3_random_region.ranked.tsv` | **main result**: N sequences |
| `survival_summary.tsv` | `total_reads_processed`, `raw_unique_sequences`, `stage1_matching_survivors`, `stage2_flanking_survivors`, `stage3_random_region_survivors`, `derived_insert_length` (no header) |
| `meme_input_top_cores.fasta` | top-count random-region cores for MEME (always written) |
| `meme_out/` | MEME's standard results, only when `run_meme: true`: `meme.xml`, `meme.html`, `meme.txt`, and per motif `logo<N>.png` / `.eps` plus reverse-complement logos `logo_rc<N>.png` / `.eps` (MEME draws these for DNA whether or not `meme_revcomp` is on) |

Each ranked table has a header line `sequence<TAB>count`, followed by one row per unique sequence. Survivor values count every surviving sequence, not just unique ones.

### Logs

`logs/aptaselect.log` (separate from the result files) records the settings, worker count and container memory limit, the start and end of each input file, periodic progress lines (reads processed, % of input, speed, ETA, per-stage survivors, unique reads so far, main-process memory), the sorting/writing of each output table, and a final timing breakdown. `logs/select_top_cores.log` records the selection (how many, count range, core lengths, tie warning) and `logs/meme.log` the MEME run.

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
| `run_meme` | run MEME in the pipeline (default `false`) |
| `meme_top_n` | number of top-count cores written/sent to MEME (default 100) |
| `meme_top_percent` | if set, take this % (0–100] of the ranked list instead of `meme_top_n` |
| `meme_nmotifs` | number of motifs (default 3) |
| `meme_mod` | `zoops` (default), `oops` or `anr` |
| `meme_minw`, `meme_maxw` | motif width range (default 6–15; `meme_maxw` ≤ `random_region_length`) |
| `meme_objfun` | `classic` (default), `de`, `se`, `cd`, `ce`, `nc` |
| `meme_revcomp` | also allow motif sites on the reverse strand (default `false`) |
| `meme_pal` | palindromes only (default `false`) |
| `meme_evt`, `meme_minsites`, `meme_maxsites`, `meme_markov_order`, `meme_seed` | optional; blank = MEME default |
| `meme_extra_args` | other MEME options; `-p`, `-o`/`-oc`, `-text` and alphabet options are rejected |

All sequence values must contain only A/C/G/T.

## Run

```bash
docker build -t aptaselect .
docker run --rm -v /path/to/fastqs:/input:ro -v /path/to/results:/output \
  -v $(pwd)/config.yaml:/pipeline/config.yaml aptaselect \
  snakemake --cores 16
```

MEME is installed in its own conda environment (`/opt/conda/envs/meme`, with ghostscript for PNG logos) separate from the Snakemake environment, and is called by its full path.

Memory: the main process keeps every unique sequence of every stage in memory; on ~100 million 150-nt reads this peaked at about 37 GB.