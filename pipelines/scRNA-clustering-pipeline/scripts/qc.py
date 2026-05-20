#!/usr/bin/env python3
"""
Quality Control for scRNA-seq data.
- Compute QC metrics (n_genes, n_counts, mito %)
- Filter cells based on thresholds
- Generate QC plots (violin plot, scatter plot, histogram)
"""

import sys
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Read config and paths
config = snakemake.config
input_file = config["input_file"]
qc_config = config["qc"]

min_genes = qc_config["min_genes"]
max_genes = qc_config["max_genes"]
max_mito_percent = qc_config["max_mito_percent"]
mito_prefix = qc_config["mito_prefix"]

# Load data
adata = sc.read(input_file, cache=True)
print(f"Loaded dataset: {adata.n_obs} cells x {adata.n_vars} genes")

# Make observation names unique if needed
if not adata.obs_names.is_unique:
    adata.obs_names_make_unique()

# Compute basic QC metrics manually
adata.obs["n_counts"] = np.array(adata.X.sum(axis=1)).flatten()
adata.obs["n_genes"] = np.array((adata.X > 0).sum(axis=1)).flatten()

# Compute mitochondrial gene percentage
mito_genes = adata.var_names.str.startswith(mito_prefix)
adata.var["mt"] = mito_genes
adata.raw = adata  # store raw counts
mito_counts = adata[:, mito_genes].X.sum(axis=1)
total_counts = adata.obs["n_counts"]
if hasattr(mito_counts, "A1"):
    mito_counts = mito_counts.A1
adata.obs["pct_counts_mt"] = (mito_counts / total_counts) * 100

# Save QC metrics
qc_metrics = adata.obs[["n_genes", "n_counts", "pct_counts_mt"]].copy()
qc_metrics.to_csv(snakemake.output[0], index=False)
print(f"QC metrics saved: {qc_metrics.shape[0]} cells, {qc_metrics.shape[1]} metrics")

# Generate QC plots
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Violin plot using seaborn
import seaborn as sns
qc_long = qc_metrics.melt(value_vars=["n_genes", "pct_counts_mt"], var_name="metric", value_name="value")
sns.violinplot(data=qc_long, x="metric", y="value", ax=axes[0], inner="quartile")
axes[0].set_title("Violin Plot")
axes[0].tick_params(axis='x', rotation=45)

# Scatter plot: n_genes vs pct_counts_mt
axes[1].scatter(adata.obs["n_genes"], adata.obs["pct_counts_mt"], s=1, alpha=0.3, c="steelblue")
axes[1].set_xlabel("Number of Genes")
axes[1].set_ylabel("% Mitochondrial")
axes[1].set_title("n_genes vs %MT")

# Histogram of n_counts
axes[2].hist(adata.obs["n_counts"], bins=100, color="steelblue", edgecolor="white")
axes[2].set_xlabel("Total Counts")
axes[2].set_ylabel("Count")
axes[2].set_title("Distribution of n_counts")

plt.tight_layout()
fig.savefig(snakemake.output[1], dpi=150, bbox_inches="tight")
plt.close(fig)
print("QC plots saved")

# Filter cells
adata_subset = adata[
    (adata.obs["n_genes"] >= min_genes) &
    (adata.obs["n_genes"] <= max_genes) &
    (adata.obs["pct_counts_mt"] <= max_mito_percent),
    :
].copy()

print(f"Before QC: {adata.n_obs} cells")
print(f"After QC: {adata_subset.n_obs} cells")
print(f"Removed {adata.n_obs - adata_subset.n_obs} cells ({100*(1 - adata_subset.n_obs/adata.n_obs):.1f}%)")

# Filter genes with 0 counts
sc.pp.filter_genes(adata_subset, min_counts=1)
print(f"After gene filtering: {adata_subset.n_obs} cells x {adata_subset.n_vars} genes")

# Save filtered data
adata_subset.write(snakemake.output[2], compression="gzip")
print("Filtered h5ad saved")