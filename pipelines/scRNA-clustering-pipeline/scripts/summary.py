#!/usr/bin/env python3
"""
Generate a summary report of the scRNA-seq analysis results.
"""

import scanpy as sc
import pandas as pd

# Read config
config = snakemake.config
output_dir = config["output_dir"]

# Load final data
adata = sc.read(snakemake.input[0])
print(f"Loaded final dataset: {adata.n_obs} cells x {adata.n_vars} genes")

# Load QC metrics
qc_metrics = pd.read_csv(snakemake.input[1])

# Compute statistics
n_clusters = len(adata.obs["leiden"].unique())
cluster_sizes = adata.obs["leiden"].value_counts().sort_index()

# Build report
lines = []
lines.append("=" * 60)
lines.append("scRNA-seq Analysis Summary Report")
lines.append("=" * 60)
lines.append("")
lines.append("DATASET OVERVIEW")
lines.append(f"  Total cells after QC: {adata.n_obs}")
lines.append(f"  Total genes: {adata.n_vars}")
lines.append(f"  Highly variable genes: {adata.var['highly_variable'].sum()}")
lines.append("")
lines.append("QC METRICS")
lines.append(f"  Median genes per cell: {qc_metrics['n_genes'].median():.0f}")
lines.append(f"  Median counts per cell: {qc_metrics['n_counts'].median():.0f}")
lines.append(f"  Median mitochondrial %: {qc_metrics['pct_counts_mt'].median():.1f}%")
lines.append("")
lines.append("PCA")
lines.append(f"  Number of components: {adata.obsm['X_pca'].shape[1]}")
lines.append("")
lines.append(f"LEIDEN CLUSTERING ({n_clusters} clusters)")
for cluster_id, count in cluster_sizes.items():
    pct = count / adata.n_obs * 100
    lines.append(f"  Cluster {cluster_id}: {count} cells ({pct:.1f}%)")
lines.append("")
lines.append("=" * 60)

report = "\n".join(lines)

with open(snakemake.output[0], "w") as f:
    f.write(report)

print(report)
print("Summary report saved")