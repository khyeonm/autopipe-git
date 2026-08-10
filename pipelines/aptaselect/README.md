# AptaSelect

AptaSelect identifies high-frequency aptamer candidate sequences from paired-end FASTQ files produced by SELEX experiments.

## Pipeline Stages

1. **Join** – Reverse-complements Read 2, finds optimal overlap with Read 1, and merges into a single sequence
2. **Selection Filtering** – Extracts the region between configurable selection primer patterns
3. **1st Sort Filtering** – Validates presence of 1st sort primer patterns within the extracted sequence
4. **2nd Sort Filtering** – Validates 2nd sort primer patterns with a required between-length constraint
5. **Aggregation & Ranking** – Deduplicates and ranks sequences by frequency at each stage
6. **Motif Analysis (MEME)** *(optional, appended)* – Extracts the variable core between two flanking motifs from the top-ranked sequences and runs MEME to find the shared motif

## Required Inputs

- `test_1.fq` — Read 1 paired-end FASTQ file
- `test_2.fq` — Read 2 paired-end FASTQ file

Gzipped FASTQ files (`.fq.gz`) are also supported.

## Outputs

- `stage1_joined_ranked.tsv` — Frequency-ranked joined sequences
- `stage2_selection_ranked.tsv` — Frequency-ranked sequences after selection filtering
- `stage3_sort1_ranked.tsv` — Frequency-ranked sequences after 1st sort filtering
- `stage4_sort2_ranked.tsv` — Frequency-ranked sequences after 2nd sort filtering
- `summary.txt` — Statistics for each stage
- `variable_cores.fasta` — *(motif step)* the variable cores extracted between the two flanking motifs; produced whenever the motif/gap values are provided, even if MEME is not run
- `meme_out/` — *(motif step)* MEME motif-discovery results (produced only when `run_meme: true`)

## How to Run

```bash
# Build the Docker image
docker build -t aptaselect /path/to/aptaselect/

# Run the pipeline
docker run --rm \
  -v /path/to/input:/input:ro \
  -v /path/to/output:/output \
  aptaselect \
  snakemake --cores 4 -s /pipeline/Snakefile
```

## Configuration

Edit `config.yaml` to adjust parameters:

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

## Motif Analysis (MEME) Step

This optional step runs **after** the count sorting. It takes the top count-ranked
sequences, extracts the **variable core** between two known flanking motifs, writes
those cores to `variable_cores.fasta`, and (optionally) runs **MEME** on them to
discover the shared motif. Only the variable core is fed to MEME — primers and other
fixed regions are excluded.

### Two-step workflow

1. **Run the count sorting first** with the shipped config (motif step off). Inspect
   the ranked TSVs.
2. **Later, run MEME** by filling in the three required values and setting
   `run_meme: true`, then re-running. Snakemake reuses the completed count-sorting
   results — they are not recomputed or deleted, and MEME results are written to a
   separate `meme_out/` folder.

### Required values (no defaults — set per experiment)

| Parameter | Description |
|---|---|
| `core_flank5` | 5' motif immediately **before** the variable core |
| `core_flank3` | 3' motif immediately **after** the variable core |
| `core_gap` | exact length (bp) of the variable core between the two flanks |

If any of these is blank, the motif step does not run. If `run_meme: true` is set
while any is still blank, the pipeline **refuses to start**.

### Adjustable motif parameters

| Parameter | Default | Description |
|---|---|---|
| `run_meme` | `false` | Master switch for the MEME step |
| `top_seq_source` | `stage4_sort2_ranked.tsv` | Which ranked table to take the top sequences from |
| `top_n` | `100` | Number of top-ranked sequences to extract cores from |
| `flank_max_mismatches` | `0` | Mismatches allowed when locating each flank (0 = exact; >0 can shift boundaries near homopolymers) |
| `meme_nmotifs` | `5` | Number of motifs for MEME to find |
| `meme_mod` | `zoops` | MEME site distribution (`oops` / `zoops` / `anr`) |
| `meme_minw` / `meme_maxw` | `6` / `20` | Motif width range (max is auto-capped to the core length) |
| `meme_revcomp` | `false` | Also search the reverse-complement strand |

### Notes on the MEME install

- MEME is installed into its **own conda environment** (`/opt/conda/envs/meme`) so it
  does not clash with the workflow engine, and is invoked by full path.
- It runs **single-process** (no `-p`), so it never starts MPI daemons.
- **ghostscript** is included so sequence logos are emitted as **PNG**, not only EPS.
- If you prefer the MEME website, just take `variable_cores.fasta` and upload it there.