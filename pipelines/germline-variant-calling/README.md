# germline-variant-calling

Calls germline variants from BAM files using GATK HaplotypeCaller in GVCF mode. Supports multiple samples in parallel, producing per-sample `.g.vcf.gz` files ready for joint genotyping.

## Inputs

| File | Description |
|------|-------------|
| `{sample}.bam` + `.bai` | Sorted, indexed BAM files |
| `reference.fa` + `.fai` + `.dict` | Reference genome FASTA with index and sequence dictionary |

## Outputs

| File | Description |
|------|-------------|
| `{sample}.g.vcf.gz` | Per-sample GVCF with raw variant calls |
| `{sample}.g.vcf.gz.tbi` | Tabix index for the GVCF |

## Configuration (`config.yaml`)

```yaml
samples:
  - NA12891_S1      # list sample names (without .bam extension)

reference: "/input/hg19_human.fa"   # path to reference FASTA inside container

java_opts: "-Xmx8g"          # Java heap size
threads: 4                   # threads for native pair-HMM
haplotypecaller_extra: ""    # extra HaplotypeCaller flags
```

## Running

```bash
# Build
docker build -t germline-variant-calling:1.0.0 .

# Run
docker run --rm \
  -v /path/to/data:/input:ro \
  -v /path/to/output:/output \
  germline-variant-calling:1.0.0 \
  snakemake --cores 4
```