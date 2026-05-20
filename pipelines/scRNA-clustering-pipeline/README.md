# scRNA-clustering-pipeline

A **Scanpy-based single-cell RNA-seq analysis pipeline** for mouse scRNA-seq data (h5ad format).

## What it does

Performs a complete scRNA-seq analysis workflow including:
1. **Quality Control** — filters cells by gene count, total counts, and mitochondrial percentage
2. **Normalization** — library size normalization + log1p transformation
3. **Highly Variable Gene Selection** — Seurat v3 flavor
4. **PCA** — Principal Component Analysis on HVGs
5. **Leiden Clustering** — kNN graph + Leiden algorithm
6. **UMAP Embedding** — 2D visualization colored by cluster

## Required Input

| File | Format | Description |
|------|--------|-------------|
| `10k.h5ad` | h5ad | Mouse scRNA-seq raw count matrix (AnnData format) |

- Mitochondrial genes must be prefixed with `mt-` (mouse convention)

## Expected Outputs

| File | Description |
|------|-------------|
| `qc_metrics.csv` | Per-cell QC metrics (n_genes, n_counts, pct_counts_mt) |
| `qc_plots.pdf` | Violin plots, scatter plots, histograms |
| `filtered.h5ad` | QC-filtered AnnData object |
| `normalized.h5ad` | Normalized + HVG-selected AnnData |
| `pca.h5ad` | PCA-computed AnnData |
| `pca_plots.pdf` | Elbow plot, PC1 vs PC2 scatter |
| `final.h5ad` | Final AnnData with clusters + UMAP coordinates |
| `clustering_plots.pdf` | UMAP colored by cluster + n_genes |
| `summary_report.txt` | Human-readable summary of analysis |

## How to Run

```bash
# Build Docker image
docker build -t scRNA-clustering-pipeline .

# Run with Docker
docker run --rm \
    -v /path/to/input:/input:ro \
    -v /path/to/output:/output \
    scRNA-clustering-pipeline \
    snakemake --cores 8
```

## Configuration

All parameters are in `config.yaml`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `qc.min_genes` | 200 | Minimum genes per cell |
| `qc.max_genes` | 6000 | Maximum genes per cell |
| `qc.max_mito_percent` | 20 | Max mitochondrial % |
| `qc.mito_prefix` | `mt-` | Mitochondrial gene prefix |
| `normalize.n_hvg` | 2000 | Number of HVGs |
| `normalize.hvg_flavor` | `seurat_v3` | HVG selection method |
| `pca.n_components` | 30 | PCA dimensions |
| `clustering.n_neighbors` | 15 | kNN neighbors |
| `clustering.resolution` | 0.5 | Leiden resolution |
| `umap.min_dist` | 0.3 | UMAP min distance |

## License

MIT