# 10X snMultiome Analysis Pipeline (Gastric WNT)

Reproduces the snMultiome (joint scRNA-seq + scATAC-seq) analysis from **"Epithelial WNT secretion drives niche escape of developing gastric cancer"** using Seurat and Signac.

## What it does

This pipeline processes 10X Genomics snMultiome data from wild-type (RZ) and Kras-mutant (RZK) mouse gastric tissue through a 7-step workflow: QC filtering, MACS2 peak calling, RNA/ATAC normalization, SCT-based integration, UMAP clustering with cell-type annotation, differential gene expression (MAST), and motif enrichment analysis (JASPAR2020).

## Required Inputs

Place input data in a directory with the following structure:

```
input_dir/
├── Mouse_RZ_WEN/                    # Wild-type sample
│   ├── filtered_feature_bc_matrix.h5
│   ├── atac_fragments.tsv.gz
│   ├── atac_fragments.tsv.gz.tbi
│   └── per_barcode_metrics.csv
└── Mouse_RZK_WEN/                   # Mutant sample
    ├── filtered_feature_bc_matrix.h5
    ├── atac_fragments.tsv.gz
    ├── atac_fragments.tsv.gz.tbi
    └── per_barcode_metrics.csv
```

Each sample requires:
- **filtered_feature_bc_matrix.h5**: 10X Cell Ranger ARC output (Gene Expression + Peaks)
- **atac_fragments.tsv.gz** + **.tbi**: ATAC fragment file with tabix index
- **per_barcode_metrics.csv**: Per-barcode QC metrics from Cell Ranger ARC

## Expected Outputs

```
output_dir/
├── rds/                             # Seurat objects at each stage
│   ├── wild_qc.rds
│   ├── mutant_qc.rds
│   ├── wild_macs2.rds
│   ├── mutant_macs2.rds
│   ├── wild_normalized.rds
│   ├── mutant_normalized.rds
│   └── integrated_clustered.rds
├── plots/                           # Publication figures
│   ├── qc_violin_plots.pdf
│   ├── fig3A_UMAP.pdf
│   ├── fig3B_dotplot.pdf
│   ├── fig3C_dittobarplot.pdf
│   ├── fig3D_Wnt7b_feature.pdf
│   ├── figS3A_featureplots.pdf
│   ├── figS3B_wnt_family_violin.pdf
│   ├── fig3I_motifs_Wnt_vs_Lgr5_mutant.pdf
│   └── fig3I_motifs_Wnt_vs_Lgr5_wildtype.pdf
├── tables/                          # Analysis results
│   ├── DEG_scMAST.csv
│   ├── KRAS_signaling_overlap.csv
│   ├── DA_peaks_Wnt_vs_Lgr5_mutant.csv
│   ├── DA_peaks_Wnt_vs_Lgr5_wildtype.csv
│   ├── enriched_motifs_mutant.csv
│   └── enriched_motifs_wildtype.csv
└── logs/                            # Step-level logs
```

## How to Run

### Build the Docker image

```bash
docker build -t snmultiome-gastric-wnt .
```

### Run the pipeline

```bash
docker run --rm \
  -v /path/to/input:/input:ro \
  -v /path/to/output:/output \
  snmultiome-gastric-wnt \
  snakemake --snakefile /pipeline/Snakefile --cores 20 -p
```

## Configuration

All parameters are in `config.yaml`. Key settings:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `qc.nCount_ATAC_min` | 1000 | Min ATAC counts per cell |
| `qc.nCount_RNA_min` | 1000 | Min RNA counts per cell |
| `qc.percent_mt_max` | 15 | Max mitochondrial % |
| `clustering_resolution_initial` | 0.3 | Initial clustering (6 clusters) |
| `clustering_resolution_refine` | 0.6 | Refined clustering (splits Pit) |
| `de_test` | MAST | DE test method |
| `workers` | 20 | Parallel workers |
| `future_ram_gb` | 80 | Max RAM for R futures |

## Pipeline Steps

1. **Load & QC** — Read 10X .h5, compute QC metrics, filter cells
2. **MACS2 Peaks** — Call peaks, remove blacklist regions, quantify
3. **Normalize** — SCTransform + LogNormalize (RNA); TF-IDF + SVD (ATAC)
4. **Integrate & Cluster** — SCT anchor integration, PCA, UMAP, clustering
5. **Visualization** — UMAP, DotPlot, dittoBarPlot, FeaturePlots
6. **Differential Expression** — MAST DEG + KRAS signaling overlap
7. **Motif Enrichment** — DA peaks (Wnt7+ vs Lgr5+) + JASPAR2020 motifs

## Cell Types Identified

| Cluster | Cell Type | Key Markers |
|---------|-----------|-------------|
| 0 | Neck | Muc6, Cftr |
| 1 | Lgr5+ | Lgr5 |
| 2 | SPEM | Glipr1 |
| 3 | Wnt7+ | Wnt7b |
| 4 | Proliferating | Mki67, Foxm1 |
| 5 | Pre-Pit | Gkn1, Tff1 |
| 6 | Pit | Muc5ac |