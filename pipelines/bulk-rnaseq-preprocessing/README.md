# Bulk RNA-seq Preprocessing Pipeline

Preprocesses bulk RNA-seq data from raw FASTQ files using the **nf-core/rnaseq** pipeline (v3.22.2) with the **STAR-Salmon** workflow, producing aligned BAMs, gene/transcript-level quantification, and a comprehensive MultiQC quality report.

## Required Inputs

Place all input files in a single directory:

1. **samplesheet.csv** — nf-core/rnaseq sample sheet describing your FASTQ files. Format:

   ```
   sample,fastq_1,fastq_2,strandedness
   SAMPLE1,/input/SAMPLE1_R1.fastq.gz,/input/SAMPLE1_R2.fastq.gz,auto
   SAMPLE2,/input/SAMPLE2_R1.fastq.gz,/input/SAMPLE2_R2.fastq.gz,auto
   ```

   - Use `/input/` prefix for all FASTQ paths (they are mounted at runtime).
   - For single-end data, leave `fastq_2` empty.
   - Set `strandedness` to `auto`, `forward`, `reverse`, or `unstranded`.

2. **FASTQ files** — Raw sequencing reads (`.fastq.gz`).

3. **(Optional) Custom references** — genome FASTA, GTF, or pre-built STAR/Salmon indices. Configure paths in `config.yaml`.

## Expected Outputs

All outputs are written to `<output_dir>/nfcore/`:

| Directory | Contents |
|---|---|
| `multiqc/star_salmon/` | MultiQC HTML report summarizing all QC metrics |
| `star_salmon/` | STAR-aligned BAM files and Salmon quantification |
| `fastqc/` | Per-sample FastQC reports |
| `trimgalore/` | Trimmed FASTQ files and trimming logs |

## Configuration

Edit `config.yaml` to customise the run:

- **genome**: Reference genome key (default: `GRCh38`). Set to `null` if providing custom FASTA/GTF.
- **igenomes_base**: Set to `null` (default) to disable iGenomes.
- **fasta / gtf / star_index / salmon_index**: Paths to custom reference files.
- **trimmer**: `trimgalore` (default) or `fastp`.
- **max_cpus / max_memory / max_time**: Resource limits per Nextflow process.
- **extra_params**: Any additional nf-core/rnaseq CLI flags.

## How to Run

### Build the Docker image

```bash
docker build -t autopipe-bulk-rnaseq-preprocessing pipelines/bulk-rnaseq-preprocessing/
```

### Execute the pipeline

```bash
docker run --rm \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v /path/to/input:/input:ro \
  -v /path/to/output:/output \
  -e HOST_INPUT_DIR=/path/to/input \
  -e HOST_OUTPUT_DIR=/path/to/output \
  -e HOST_PIPELINE_DIR=/path/to/pipeline \
  autopipe-bulk-rnaseq-preprocessing \
  snakemake --cores 16 -p
```

> **Note:** The Docker socket mount (`-v /var/run/docker.sock`) is required because nf-core/rnaseq pulls its own tool containers at runtime.