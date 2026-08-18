# ssam-spatial-celltyping

Cell segmentation-free spatial cell-type inference using [SSAM](https://github.com/pnucolab/ssam) (Spot-based Spatial cell-type Analysis by Multidimensional mRNA density estimation). Given single-molecule mRNA coordinates, SSAM applies KDE to build a gene-expression density map, finds representative expression vectors at local maxima, clusters them *de novo*, and maps a cell type to every pixel of the tissue.

This pipeline follows the canonical SSAM workflow (see the project's [tl;dr guide](https://github.com/pnucolab/ssam/blob/develop/docs/tldr.rst)).

## Inputs

| File | Format | Description |
|------|--------|-------------|
| spot table CSV | CSV | Single-molecule mRNA locations. **3 columns: gene, x, y** (column names auto-detected, e.g. `feature_name`, `x_location`/`y_location`) |

### CSV example
```
gene,x,y
Slc17a7,123.4,456.7
Gad1,234.5,123.4
...
```

## Outputs

| File | Description |
|------|-------------|
| `celltype_map.png` | Spatial *de novo* cell-type assignment map (background masked) |
| `kde_map.png` | KDE gene-expression density map |
| `celltype_abundance.csv` | Per-cluster pixel count and fraction |

## Configuration (`config.yaml`)

Defaults reproduce the canonical SSAM reference workflow.

```yaml
csv:               "/input/s3_spot_table.csv"  # required — spot table path
bandwidth:         2.5    # KDE bandwidth in µm (ssam default: 2.5)
sampling_distance: 1.0    # KDE grid spacing in µm/pixel (ssam default: 1.0)
norm_threshold:    ""     # empty → ssam's bandwidth-derived default; set a number to override
search_size:       3      # find_localmax neighborhood size (ssam default: 3)
resolution:        0.6    # Leiden clustering resolution (canonical: 0.6)
min_norm:          0.05   # filter_celltypemaps background masking (canonical: 0.05)
threads:           4
```

## Running

```bash
# Build
docker build -t ssam-spatial-celltyping:1.0.2 .

# Run
docker run --rm \
  -v /path/to/data:/input:ro \
  -v /path/to/output:/output \
  ssam-spatial-celltyping:1.0.2 \
  snakemake --cores 4
```

## Reference

Park et al. "Cell segmentation-free inference of cell types from in situ transcriptomics data." *Nature Communications* 12, 3545 (2021).
