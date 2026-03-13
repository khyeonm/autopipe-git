# scLENS Pipeline

scLENS (Single-cell Low-dimension Embedding using effective Noise Subtraction) is a dimensionality reduction tool for scRNA-seq data that automatically detects biologically meaningful signals using random matrix theory (RMT), without requiring manual parameter tuning.

## Required Inputs

- **CSV file**: scRNA-seq count matrix where rows = cells, columns = genes. First row must contain gene names/IDs, first column must contain cell IDs.

## Expected Outputs

| File | Description |
|------|-------------|
| `pca.csv` | PCA embedding (noise-filtered dimensions) |
| `umap.csv` | UMAP 2D coordinates |
| `test_data.h5ad` | Full result in AnnData format (compatible with Scanpy) |
| `umap_dist.png` | UMAP scatter plot colored by cluster |
| `sclens.log` | Run log |

## Configuration (`config.yaml`)

```yaml
input_csv: "/input/Z8eq_l.csv"   # path to input CSV (required)
threads: 8                         # number of CPU threads
device: "cpu"                      # "cpu" or "gpu"
```

## Citation

Kim, H., Chang, W., Chae, S.J. et al. scLENS: data-driven signal detection for unbiased scRNA-seq data analysis. *Nat Commun* 15, 3575 (2024). https://doi.org/10.1038/s41467-024-47884-3