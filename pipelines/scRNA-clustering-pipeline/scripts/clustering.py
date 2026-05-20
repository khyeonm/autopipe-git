#!/usr/bin/env python3
"""
Leiden clustering and UMAP embedding for scRNA-seq data.
- Build kNN graph
- Leiden clustering
- UMAP embedding
- Generate UMAP plots colored by cluster
"""

import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# Read config
config = snakemake.config
clustering_config = config["clustering"]
umap_config = config["umap"]

# Load PCA data
adata = sc.read(snakemake.input[0])
print(f"Loaded PCA dataset: {adata.n_obs} cells x {adata.n_vars} genes")

# Build neighborhood graph
n_neighbors = clustering_config["n_neighbors"]
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=30)
print(f"Neighborhood graph built with k={n_neighbors}")

# Leiden clustering
resolution = clustering_config["resolution"]
random_state = clustering_config["random_state"]
sc.tl.leiden(adata, resolution=resolution, random_state=random_state)
n_clusters = len(adata.obs["leiden"].unique())
print(f"Leiden clustering: {n_clusters} clusters (resolution={resolution})")

# UMAP embedding (uses pre-computed neighbors from sc.pp.neighbors)
umap_min_dist = umap_config["min_dist"]
umap_random_state = umap_config["random_state"]
sc.tl.umap(adata, min_dist=umap_min_dist, random_state=umap_random_state)
print("UMAP embedding computed")

# Generate UMAP plots using matplotlib directly
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# UMAP colored by cluster
cluster_colors = sns.color_palette("tab20", n_clusters)
cluster_dict = {str(c): cluster_colors[int(c) % n_clusters] for c in adata.obs["leiden"].cat.categories}
for cluster_id in adata.obs["leiden"].cat.categories:
    mask = adata.obs["leiden"] == cluster_id
    axes[0].scatter(
        adata.obsm["X_umap"][mask, 0],
        adata.obsm["X_umap"][mask, 1],
        s=5, alpha=0.6, c=cluster_dict[cluster_id], label=f"Cluster {cluster_id}", edgecolors='none'
    )
axes[0].set_xlabel("UMAP 1")
axes[0].set_ylabel("UMAP 2")
axes[0].set_title(f"UMAP (Leiden, {n_clusters} clusters)")
axes[0].legend(markerscale=5, fontsize=8, loc="best")
axes[0].set_axis_off()

# UMAP colored by n_genes
scatter_obj = axes[1].scatter(
    adata.obsm["X_umap"][:, 0],
    adata.obsm["X_umap"][:, 1],
    s=5, c=adata.obs["n_genes"], cmap="viridis", alpha=0.6, edgecolors='none'
)
plt.colorbar(scatter_obj, ax=axes[1], label="Number of Genes")
axes[1].set_xlabel("UMAP 1")
axes[1].set_ylabel("UMAP 2")
axes[1].set_title("UMAP (n_genes)")
axes[1].set_axis_off()

plt.tight_layout()
fig.savefig(snakemake.output[1], dpi=150, bbox_inches="tight")
plt.close(fig)
print("Clustering plots saved")

# Save final data
adata.write(snakemake.output[0], compression="gzip")
print("Final h5ad saved")