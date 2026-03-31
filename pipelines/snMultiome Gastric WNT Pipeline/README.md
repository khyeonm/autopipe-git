# snMultiome Gastric WNT Pipeline

A reproducible 10X Genomics snMultiome (joint snRNA-seq + scATAC-seq) analysis pipeline for the study **"Epithelial WNT secretion drives niche escape of developing gastric cancer"**. Compares wild-type (RZ) vs Kras-mutant (RZK) mouse gastric epithelial cells using Seurat and Signac.

## Required Inputs

Place input data under a single directory with the following structure:

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

These are standard CellRanger-ARC outputs for each sample.

## Expected Outputs

| File | Description |
|------|-------------|
| `01_wild_qc.rds`, `01_mutant_qc.rds` | QC-filtered Seurat objects |
| `02_wild_peaks.rds`, `02_mutant_peaks.rds` | Objects with MACS2 peak assay |
| `03_wild_normalized.rds`, `03_mutant_normalized.rds` | Normalized objects (SCT + TF-IDF) |
| `04_integrated_clustered.rds` | Final integrated object with cell-type labels |
| `05_figures/fig3A_UMAP.pdf` | UMAP colored by cluster × genotype |
| `05_figures/fig3B_dotplot.pdf` | DotPlot of marker genes |
| `05_figures/fig3C_dittobar.pdf` | Cell-type proportion barplot |
| `05_figures/fig3D_Wnt7b_feature.pdf` | Wnt7b FeaturePlot split by genotype |
| `05_figures/fig3H_coverage_Wnt7b.pdf` | CoveragePlot with peak-gene linkage |
| `05_figures/fig3I_motif_enrichment_*.pdf` | TF motif enrichment bar plots |
| `05_figures/figS3A_feature_plots.pdf` | Individual marker FeaturePlots |
| `05_figures/figS3B_wnt_violin.pdf` | Wnt family VlnPlot |
| `06_DEG_MAST.csv` | Differential expression results (MAST) |
| `06_DEG_KRAS_overlap.csv` | DEGs overlapping KRAS_SIGNALING_UP pathway |
| `07_DA_peaks_wnt_vs_lgr5_*.csv` | Differential accessibility peaks |
| `07_motif_enrichment_*.csv` | Enriched TF motifs from JASPAR2020 |

## Pipeline Steps

1. **Load & QC** — Read 10X multiome h5, compute QC metrics (nucleosome signal, TSS enrichment, %mt), filter cells
2. **MACS2 Peak Calling** — Call peaks de novo, remove blacklisted regions, quantify new peak matrix
3. **Normalization** — SCTransform + LogNormalize (RNA); TF-IDF + SVD (ATAC + MACS2)
4. **Integration & Clustering** — SCT integration, PCA, UMAP, Louvain clustering (res 0.3 + 0.6), cell-type annotation
5. **Visualization** — UMAP, DotPlot, dittoBarPlot, FeaturePlot, CoveragePlot, VlnPlot
6. **Differential Expression** — MAST test (Mutant vs Wild-type), KRAS pathway overlap via MSigDB/gprofiler2
7. **Motif Enrichment** — DA peaks (Wnt7+ vs Lgr5+), JASPAR2020 motif scanning, enrichment visualization

## QC Thresholds

| Metric | Threshold |
|--------|-----------|
| nCount_ATAC | 1,000 – 100,000 |
| nCount_RNA | 1,000 – 25,000 |
| Nucleosome signal | < 2 |
| TSS enrichment | > 1 |
| Mitochondrial % | < 15% |

## How to Run

```bash
# Build the Docker image
docker build -t snmultiome-gastric-wnt .

# Run the pipeline
docker run --rm \
  -v /path/to/input:/input:ro \
  -v /path/to/output:/output \
  snmultiome-gastric-wnt \
  conda run --no-capture-output -n pipeline \
  snakemake --snakefile /pipeline/Snakefile --cores 20
```

## Configuration

All parameters are in `config.yaml`. Key settings include QC thresholds, clustering resolutions (0.3 initial, 0.6 refined), marker gene lists, DE test method (MAST), DA test method (LR), and JASPAR motif collection. See the config file for full documentation.

## Cell Types Identified

| Cluster | Cell Type | Key Markers |
|---------|-----------|-------------|
| 0 | Neck | Muc6, Cftr |
| 1 | Lgr5+ | Lgr5 |
| 2 | SPEM | Glipr1 |
| 3 | Wnt7+ | Wnt7b |
| 4 | Proliferating | Mki67, Foxm1, Top2a |
| 5 | Pre-Pit | Gkn1, Tff1 |
| 6 | Pit | Muc5ac, Gkn2 |