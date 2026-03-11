# ssam-spatial-celltyping

Cell segmentation-free spatial cell-type inference using [SSAM](https://github.com/pnucolab/ssam) (Spot-based Spatial cell-type Analysis by Multidimensional mRNA density estimation). Given single-molecule mRNA coordinates and reference cell-type signatures, SSAM applies KDE to build a gene expression density map and assigns cell types to every pixel of the tissue image.

## Inputs

| File | Format | Description |
|------|--------|-------------|
| `coordinates.csv` | CSV | mRNA molecule locations. **3 columns: gene, x, y** (column names flexible) |
| `signatures.csv` | CSV | Cell-type reference signatures. **Rows = cell types, Columns = genes** |

### coordinates.csv example
```
gene,x,y
Slc17a7,123.4,456.7
Gad1,234.5,123.4
...
```

### signatures.csv example
```
,Slc17a7,Gad1,Pvalb,...
ExcitatoryNeuron,12.3,0.1,0.0,...
InhibitoryNeuron,0.2,8.5,6.1,...
...
```

## Outputs

| File | Description |
|------|-------------|
| `celltype_map.png` | Spatial cell-type assignment map |
| `kde_map.png` | KDE gene expression density map |
| `ssam_result.h5ad` | Full AnnData result for downstream analysis |
| `celltype_abundance.csv` | Per-cell-type pixel count and fraction |

## Configuration (`config.yaml`)

```yaml
coordinates: "/input/coordinates.csv"   # required
signatures:  "/input/signatures.csv"    # required
bandwidth:   2.5     # KDE bandwidth in µm (larger = smoother)
threshold:   0.2     # expression threshold (0–1)
map_width:   2000    # output map width in pixels
threads:     4
```

## Running

```bash
# Build
docker build -t ssam-spatial-celltyping:1.0.0 .

# Run
docker run --rm \
  -v /path/to/data:/input:ro \
  -v /path/to/output:/output \
  ssam-spatial-celltyping:1.0.0 \
  snakemake --cores 4
```

## Reference

Park et al. "Cell segmentation-free inference of cell types from in situ transcriptomics data." *Nature Communications* 12, 3545 (2021).