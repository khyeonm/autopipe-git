# snMultiome Gastric WNT Pipeline

A 7-step Snakemake pipeline for 10X Genomics snMultiome (joint snRNA-seq + scATAC-seq) analysis, reproducing the computational results from *"Epithelial WNT secretion drives niche escape of developing gastric cancer"*. Processes wild-type (RZ) and KRAS-mutant (RZK) mouse gastric tissue through QC, peak calling, normalisation, integration, clustering, differential analysis, and motif enrichment.

## Required Inputs

Place your Cell Ranger ARC outputs in a directory with the following structure:

```
input_data/
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

## Expected Outputs

```
output/
├── 01_qc/                           # QC'd Seurat objects + violin plots
├── 02_peaks/                        # MACS2 peak-called objects
├── 03_normalised/                   # SCTransform + TF-IDF normalised objects
├── 04_integrated/                   # Integrated + clustered object (RZ_RZK_integrated.rds)
├── 05_figures/                      # Publication figures (PDF)
│   ├── fig3A_UMAP.pdf
│   ├── fig3B_dotplot.pdf
│   ├── fig3C_dittobarplot.pdf
│   ├── fig3D_Wnt7b_featureplot.pdf
│   ├── figS3A_marker_featureplots.pdf
│   └── figS3B_wnt_family_violin.pdf
├── 06_differential/                 # DE/DA results (CSV)
│   ├── DEG_MAST_results.csv
│   ├── KRAS_signaling_overlap.csv
│   ├── DA_peaks_Wnt7_vs_Lgr5_mutant.csv
│   ├── DA_peaks_Wnt7_vs_Lgr5_wildtype.csv
│   ├── Wnt7b_expression_stats.csv
│   └── Wnt7b_accessibility_stats.csv
├── 07_motifs/                       # Motif enrichment results + bar plots
│   ├── motif_enrichment_Wnt_vs_Lgr5_mutant.csv
│   ├── motif_enrichment_Wnt_vs_Lgr5_wildtype.csv
│   ├── fig3I_motif_barplot_mutant.pdf
│   └── fig3I_motif_barplot_wildtype.pdf
└── logs/                            # Per-step log files
```

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_load_and_qc.R` | Load 10X h5, compute QC metrics (nucleosome signal, TSS enrichment, %mito), filter cells |
| 2 | `02_macs2_peaks.R` | Call peaks with MACS2, remove blacklist regions, create `macs2` ChromatinAssay |
| 3 | `03_normalize.R` | SCTransform (RNA) + log-normalise (RNA) + TF-IDF/SVD (ATAC & macs2) |
| 4 | `04_integrate_cluster.R` | SCT integration, PCA, UMAP, clustering at res 0.3 & 0.6, cell-type annotation |
| 5 | `05_visualization.R` | UMAP, DotPlot, dittoBarPlot, Wnt7b FeaturePlot, marker FeaturePlots, Wnt family VlnPlot |
| 6 | `06_differential.R` | DEG (MAST), KRAS signaling overlap, DA peaks (LR), Wnt7b expression/accessibility stats |
| 7 | `07_motif_enrichment.R` | JASPAR2020 TF motif enrichment on DA peaks, bar plots |

## QC Thresholds

| Metric | Filter |
|--------|--------|
| nCount_ATAC | 1,000 – 100,000 |
| nCount_RNA | 1,000 – 25,000 |
| Nucleosome signal | < 2 |
| TSS enrichment | > 1 |
| Percent mitochondrial | < 15% |

## How to Run

### Build the Docker image

```bash
docker build -t snmultiome-gastric-wnt .
```

### Run the pipeline

```bash
docker run --rm \
  -v /path/to/your/input_data:/input:ro \
  -v /path/to/output:/output \
  snmultiome-gastric-wnt \
  snakemake --cores 20 --snakefile /pipeline/Snakefile -d /pipeline
```

## Configuration

All parameters are in `config.yaml`. Key options:

- **QC thresholds**: `qc.max_nCount_ATAC`, `qc.min_nCount_RNA`, etc.
- **Clustering resolution**: `cluster_resolution_initial` (0.3), `cluster_resolution_refined` (0.6)
- **DE test**: `de_test` (MAST), `de_padj_cutoff` (0.01), `de_log2fc_cutoff` (1.0)
- **DA test**: `da_test` (LR), `da_padj_cutoff` (0.005)
- **Motif enrichment**: `jaspar_collection` (CORE), `motif_padj_cutoff` (0.05)
- **Parallelism**: `threads` (20), `max_memory_gb` (80)
- **Cluster labels**: mapped in `cluster_labels` (0→Neck, 1→Lgr5+, etc.)

## Cell Type Annotations

| Cluster | Cell Type | Key Markers |
|---------|-----------|-------------|
| 0 | Neck | Muc6, Cftr |
| 1 | Lgr5+ | Lgr5 |
| 2 | SPEM | Glipr1 |
| 3 | Wnt7+ | Wnt7b |
| 4 | Proliferating | Mki67, Foxm1, Top2a, Smc2 |
| 5 | Pre-Pit | Gkn1, Tff1 |
| 6 | Pit | Muc5ac, Gkn2 |