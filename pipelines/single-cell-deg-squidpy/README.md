# single-cell-deg-squidpy

This AutoPipe workflow performs single-cell differential expression analysis from an AnnData `.h5ad` file. It uses Scanpy for marker testing and includes Squidpy spatial neighborhood analysis when spatial coordinates are available in `adata.obsm['spatial']`.

## Inputs

- `/input/sample.h5ad`: AnnData object containing a count/expression matrix and cell metadata in `adata.obs`.
- `groupby`: an `adata.obs` column that defines groups for differential expression, such as `cell_type`, `leiden`, `condition`, or `sample_group`.

## Outputs

- `/output/deg_results.tsv`: full differential expression table with group, gene, score, p-value, adjusted p-value, log fold change, and expression fractions when available.
- `/output/top_markers_by_group.tsv`: top ranked markers per group.
- `/output/analysis_summary.json`: run metadata, dimensions, settings, groups, and generated output list.
- `/output/qc_overview.png`: basic group-size and library-size QC plot.
- `/output/deg_dotplot.png`: dot plot of top marker genes by group.
- `/output/squidpy_neighborhood_enrichment.tsv`: Squidpy cluster neighborhood enrichment matrix when spatial coordinates are present.
- `/output/squidpy_spatial_autocorr.tsv`: Squidpy Moran's I spatial autocorrelation results for top variable genes when spatial coordinates are present.
- `/output/single_cell_deg.log`: command log.

## Configuration

Edit `config.yaml` before running:

- `input_h5ad`: path to the mounted AnnData file, usually `/input/<file>.h5ad`.
- `groupby`: required observation column for DEG groups.
- `method`: Scanpy DEG method, commonly `wilcoxon`, `t-test`, `t-test_overestim_var`, or `logreg`.
- `reference`: `rest` or a specific category from the `groupby` column.
- `use_raw`: whether to use `adata.raw` for DEG.
- `layer`: optional layer name for DEG.
- `normalize_total` and `log1p`: enable standard normalization before testing.
- `spatial_n_neighs`, `spatial_radius`, and `spatial_coord_type`: Squidpy spatial graph options.

## Run With Docker

```bash
docker build -t single-cell-deg-squidpy .
docker run --rm \
  -v /path/to/input:/input:ro \
  -v /path/to/output:/output \
  single-cell-deg-squidpy \
  snakemake --cores 4 --snakefile /pipeline/Snakefile
```

The pipeline never modifies input data; `/input` is mounted read-only and all results are written to `/output`.