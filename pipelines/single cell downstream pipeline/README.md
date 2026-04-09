# Single Cell Downstream Pipeline

Single cell RNA-seq downstream analysis pipeline using Scanpy. Takes raw count h5ad data and performs full downstream analysis from QC to clustering.

## Analysis Steps

1. **QC & Filtering** — filter cells by min genes, min/max counts, mitochondrial gene %
2. **Normalization** — total-count normalization (target_sum=1e4) + log1p transformation
3. **Highly Variable Genes** — select top 2000 HVGs using Seurat v3 method
4. **PCA** — principal component analysis (50 components)
5. **Neighbors** — k-nearest neighbor graph (k=15)
6. **UMAP** — 2D embedding for visualization
7. **Leiden Clustering** — community detection (resolution=1.0)
8. **Visualization** — UMAP plots colored by cluster and QC metrics

## Required Input

- `.h5ad` file with raw counts (AnnData format)

## Expected Outputs

- `umap_leiden.png` — UMAP plot colored by Leiden clusters
- `umap_qc_metrics.png` — UMAP colored by QC metrics (n_genes, total_counts, mito%)
- `qc_violin.png` — Violin plots of QC metrics
- `hvg_plot.png` — Highly variable genes plot
- `pca_elbow.png` — PCA variance ratio (elbow) plot
- `processed.h5ad` — Processed AnnData object
- `analysis.log` — Full analysis log

## Configuration (config.yaml)

| Parameter | Default | Description |
|-----------|---------|-------------|
| min_genes | 200 | Min genes per cell |
| min_cells | 3 | Min cells per gene |
| max_mito_pct | 20 | Max mitochondrial gene % |
| min_counts | 500 | Min total counts per cell |
| max_counts | 50000 | Max total counts per cell |
| n_top_genes | 2000 | Number of HVGs |
| n_pcs | 50 | Number of PCs |
| n_neighbors | 15 | k for kNN graph |
| resolution | 1.0 | Leiden clustering resolution |

## How to Run

```bash
# Build Docker image
docker build -t autopipe-single-cell-downstream-pipeline .

# Run pipeline
docker run --rm \
  -v /path/to/input:/input:ro \
  -v /path/to/output:/output \
  autopipe-single-cell-downstream-pipeline \
  snakemake --cores 4 -p
```