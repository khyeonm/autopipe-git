# Scanpy Analysis Pipeline

Single-cell RNA-seq downstream analysis pipeline using Scanpy. This pipeline performs comprehensive quality control, normalization, dimensionality reduction, clustering, and visualization of single-cell transcriptomic data.

## Features

- **Quality Control**: Cell and gene filtering based on gene counts and mitochondrial percentage
- **Normalization**: Log-normalization with configurable target sum
- **Feature Selection**: Highly variable gene detection using Seurat V3 method
- **Dimensionality Reduction**: PCA and UMAP embeddings
- **Clustering**: Leiden clustering with configurable resolution
- **Differential Expression**: Wilcoxon-based marker gene detection per cluster
- **Visualization**: UMAP plots, marker gene heatmaps, PCA variance plots
- **Summary Statistics**: Cluster composition and QC metrics

## Input

- **Required**: An `.h5ad` file containing raw count matrix
  - Should have cells as observations and genes as variables
  - Gene names should be properly annotated (e.g., ENSEMBL or gene symbols)
  - Place file at `/input/data.h5ad` (or configure via `config.yaml`)

## Output

- `annotated.h5ad`: Fully processed AnnData object with all analysis results
- `plots/umap_clusters.png/pdf`: UMAP visualization colored by clusters
- `plots/marker_heatmap.png/pdf`: Heatmap of top marker genes per cluster
- `plots/pca_variance.png/pdf`: PCA variance explained plot
- `plots/cluster_stats.csv`: Cluster composition and statistics
- `summary_qc_report.csv`: Quality control summary metrics

## Configuration

Edit `config.yaml` to customize parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_genes` | 200 | Minimum genes per cell |
| `max_genes` | 7500 | Maximum genes per cell |
| `max_mito_percent` | 20.0 | Maximum mitochondrial gene percentage |
| `n_hvg` | 2000 | Number of highly variable genes |
| `n_pcs` | 50 | Number of PCA components |
| `n_neighbors` | 15 | KNN graph neighbors |
| `resolution` | 0.5 | Leiden clustering resolution |
| `de_n_top` | 50 | Top marker genes per cluster |

## Usage

```bash
# Build Docker image
docker build -t scanpy-analysis .

# Run pipeline
docker run --rm \
  -v /path/to/input:/input:ro \
  -v /path/to/output:/output \
  scanpy-analysis snakemake --cores 8
```

## Analysis Steps

1. **QC**: Filter cells/genes, calculate mitochondrial percentage
2. **Normalization**: Log-normalize counts
3. **HVG Selection**: Identify top variable genes
4. **PCA**: Principal component analysis
5. **KNN Graph**: Build neighborhood graph
6. **Clustering**: Leiden algorithm
7. **UMAP**: 2D embedding
8. **DE Analysis**: Marker gene detection
9. **Annotation**: Cluster statistics
10. **Visualization**: Generate plots

## Requirements

- Docker
- 8+ CPU cores recommended
- Memory: 16GB minimum, 32GB+ recommended for large datasets

## License

MIT