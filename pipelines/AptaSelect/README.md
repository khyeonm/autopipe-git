# AptaSelect

AptaSelect identifies high-frequency aptamer candidate sequences from paired-end
FASTQ files produced by SELEX experiments. It runs three sequential filtering
stages followed by aggregation and ranking. An appended **MEME motif-analysis
step** can then find the shared motif in the top-ranked variable cores. Every
pattern sequence, mismatch tolerance and length is user-configurable, so the same
pipeline can be reused for different library designs by editing only
`config.yaml`.

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

## Motif analysis (MEME) — appended step

After AptaSelect, an optional step runs MEME on the **variable random core**:

- It reads the final ranked random-region table
  (`stage3_random_region.ranked.tsv`), which is already trimmed to just the
  random core — primers, constant regions and stems removed. Only this variable
  core is fed to MEME; the constant/primer regions are never included, because
  they would bias motif discovery.
- It selects the top-count cores (the ranked table is already count-sorted):
  either a fixed number (`meme_top_n`) or, when `meme_top_percent > 0`, that
  percentage of the ranked list (falling back to `meme_top_n` otherwise).
- It **always** writes the selected cores to `top_cores.fasta` as a real default
  output, so you can upload it to the MEME web server yourself even if you don't
  run MEME here.
- When `run_meme: true`, it runs MEME as the plain **serial** binary (called by
  full path, DNA alphabet, **no** `-p`/MPI, which fails in the container) and
  writes results into `meme_out/`. ghostscript is included so logos are written
  as PNG, not only EPS.

**The MEME step is OFF by default.** Run the count sorting first, check
`stage3_random_region.ranked.tsv`, then set `run_meme: true` to run MEME later.
The MEME step only reads the stage3 output; it does not re-run or delete the
earlier AptaSelect results.

> Note: a MEME motif is only a **candidate**. Afterwards the full sequence should
> be checked separately for structure (G-quadruplex / stem-loop) and folding
> energy to rule out motifs that are merely SELEX bias (e.g. cores complementary
> to the constant regions). That structural check is done separately, **not** in
> this pipeline.

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
- `top_cores.fasta` — selected top-count cores (only when the MEME step is
  enabled)
- `meme_out/` — MEME results including `meme.html` and sequence logos (only when
  the MEME step is enabled)

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
| `run_meme` | Enable the MEME motif step (default `false`) |
| `meme_top_n` | Fixed number of top cores to send to MEME |
| `meme_top_percent` | If `> 0`, take this % of the ranked list instead of `meme_top_n` |
| `meme_nmotifs` | Number of motifs for MEME to find |
| `meme_mod` | MEME site distribution: `oops` / `zoops` / `anr` |
| `meme_minw`, `meme_maxw` | Min / max motif width (`0` = MEME default) |
| `meme_bin` | Full path to the serial MEME binary |
| `meme_extra` | Extra raw args passed through to MEME |

## Running

```bash
docker build -t aptaselect .

# 1) Count sorting only (MEME off — the default). Check stage3 afterwards.
docker run --rm \
    -v /path/to/inputs:/input:ro \
    -v /path/to/outputs:/output \
    aptaselect \
    snakemake --cores 1

# 2) Later, enable MEME (set run_meme: true in config.yaml) and re-run.
#    Earlier results are kept; only the motif step is added.
```