# AptaSelect

Identifies high-frequency aptamer candidate sequences from paired-end FASTQ files produced by SELEX experiments.

## Pipeline Stages

1. **Join** — Reverse-complements Read 2, finds optimal overlap, merges paired reads using quality scores
2. **Selection Filtering** — Locates selection adapter patterns and extracts the enclosed region
3. **1st Sort Filtering** — Validates the extracted sequence contains 1st sort patterns with correct spacing (40 bp)
4. **2nd Sort Filtering** — Validates the extracted sequence contains 2nd sort patterns with correct spacing (20 bp)
5. **Aggregation & Ranking** — Counts identical sequences and outputs the top N candidates

## Required Inputs

- `{sample}_1.fq` — Read 1 FASTQ file
- `{sample}_2.fq` — Read 2 FASTQ file

Place input files in the input directory. Sample names are configured in `config.yaml`.

## Expected Outputs

- `{sample}_top_candidates.tsv` — Ranked aptamer candidates with counts
- `{sample}_pipeline_stats.txt` — Statistics from all pipeline stages

## Configuration

Edit `config.yaml` to adjust parameters including overlap settings, adapter patterns, sort pattern lengths, mismatch tolerance, and the number of top candidates to report.

## Running

```bash
# Build the Docker image
docker build -t autopipe-aptaselect /path/to/aptaselect/

# Run the pipeline
docker run --rm \
  -v /path/to/input:/input:ro \
  -v /path/to/output:/output \
  autopipe-aptaselect \
  snakemake --cores 4 -s /pipeline/Snakefile --directory /output
```