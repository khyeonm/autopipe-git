#!/usr/bin/env python3
"""
PCA analysis for scRNA-seq data.
- Run PCA on highly variable genes
- Generate PCA plots (elbow plot, PCA scatter)
"""

import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Read config
config = snakemake.config
pca_config = config["pca"]

# Load normalized data
adata = sc.read(snakemake.input[0])
print(f"Loaded normalized dataset: {adata.n_obs} cells x {adata.n_vars} genes")

# Run PCA on highly variable genes
n_components = pca_config["n_components"]
sc.tl.pca(adata, svd_solver="arpack", n_comps=n_components)
print(f"PCA computed: {adata.obsm['X_pca'].shape[1]} components")

# Generate PCA plots
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Elbow plot: variance explained per PC
variance_ratio = adata.uns["pca"]["variance_ratio"]
axes[0].plot(range(1, len(variance_ratio) + 1), variance_ratio, marker='o', markersize=3, linestyle='-', color='steelblue')
axes[0].set_xlabel("Principal Component")
axes[0].set_ylabel("Variance Explained")
axes[0].set_title("PCA Elbow Plot")
axes[0].grid(True, alpha=0.3)

# PC1 vs PC2 scatter
axes[1].scatter(adata.obsm["X_pca"][:, 0], adata.obsm["X_pca"][:, 1], s=1, alpha=0.3, c="steelblue")
axes[1].set_xlabel(f"PC1 ({variance_ratio[0]*100:.1f}%)")
axes[1].set_ylabel(f"PC2 ({variance_ratio[1]*100:.1f}%)")
axes[1].set_title("PC1 vs PC2")
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(snakemake.output[1], dpi=150, bbox_inches="tight")
plt.close(fig)
print("PCA plots saved")

# Save PCA data
adata.write(snakemake.output[0], compression="gzip")
print("PCA h5ad saved")