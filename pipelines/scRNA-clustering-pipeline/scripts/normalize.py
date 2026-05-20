#!/usr/bin/env python3
"""
Normalization and Highly Variable Gene selection for scRNA-seq data.
- Library size normalization
- Log1p transformation
- Highly variable gene selection (Seurat v3 flavor)
"""

import scanpy as sc

# Read config
config = snakemake.config
norm_config = config["normalize"]

# Load filtered data
adata = sc.read(snakemake.input[0])
print(f"Loaded filtered dataset: {adata.n_obs} cells x {adata.n_vars} genes")

# Store raw counts
adata.raw = adata

# Library size normalization (counts per cell divided by total counts, multiplied by 1e4)
sc.pp.normalize_total(adata, target_sum=1e4)

# Log1p transformation
sc.pp.log1p(adata)

# Highly variable gene selection
n_hvg = norm_config["n_hvg"]
hvg_flavor = norm_config["hvg_flavor"]

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=n_hvg,
    flavor=hvg_flavor,
    subset=False
)

print(f"Selected {adata.var['highly_variable'].sum()} highly variable genes")

# Save normalized data
adata.write(snakemake.output[0], compression="gzip")
print("Normalized h5ad saved")