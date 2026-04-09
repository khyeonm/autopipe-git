#!/usr/bin/env python
"""Single cell downstream analysis pipeline.

Steps:
1. QC & Filtering (min genes, min cells, mito%, total counts)
2. Normalization & Log transformation
3. Highly Variable Genes (HVGs) selection
4. PCA
5. Neighbors graph
6. UMAP embedding
7. Leiden clustering
8. Save UMAP plot and processed h5ad
"""

import sys
import yaml
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    # Load config
    with open("/pipeline/config.yaml", "r") as f:
        config = yaml.safe_load(f)

    input_path = config["input_h5ad"]
    output_dir = config["output_dir"]
    
    # QC parameters
    min_genes = config.get("min_genes", 200)
    min_cells = config.get("min_cells", 3)
    max_mito_pct = config.get("max_mito_pct", 20)
    min_counts = config.get("min_counts", 500)
    max_counts = config.get("max_counts", 50000)
    
    # Analysis parameters
    n_top_genes = config.get("n_top_genes", 2000)
    n_pcs = config.get("n_pcs", 50)
    n_neighbors = config.get("n_neighbors", 15)
    resolution = config.get("resolution", 1.0)
    random_state = config.get("random_state", 42)

    # --------------------------------------------------
    # 1. Load data
    # --------------------------------------------------
    print(f"Loading data from {input_path} ...")
    adata = sc.read_h5ad(input_path)
    print(f"  Raw data shape: {adata.shape}")

    # Make observation and variable names unique
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    # --------------------------------------------------
    # 2. QC & Filtering
    # --------------------------------------------------
    print("Running QC & Filtering ...")
    
    # Basic filtering
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # Mitochondrial gene QC
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    
    # Filter by mito% and total counts
    adata = adata[adata.obs["pct_counts_mt"] < max_mito_pct, :].copy()
    adata = adata[adata.obs["total_counts"] > min_counts, :].copy()
    adata = adata[adata.obs["total_counts"] < max_counts, :].copy()
    
    print(f"  After QC: {adata.shape}")

    # Save QC violin plot
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    sc.pl.violin(adata, keys=["n_genes_by_counts"], ax=axes[0], show=False)
    sc.pl.violin(adata, keys=["total_counts"], ax=axes[1], show=False)
    sc.pl.violin(adata, keys=["pct_counts_mt"], ax=axes[2], show=False)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/qc_violin.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved qc_violin.png")

    # --------------------------------------------------
    # 3. Normalization & Log transformation
    # --------------------------------------------------
    print("Normalizing and log-transforming ...")
    # Store raw counts
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # --------------------------------------------------
    # 4. Highly Variable Genes
    # --------------------------------------------------
    print(f"Selecting top {n_top_genes} highly variable genes ...")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat_v3", layer="counts")
    
    # Save HVG plot
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(f"{output_dir}/hvg_plot.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved hvg_plot.png")

    # --------------------------------------------------
    # 5. PCA
    # --------------------------------------------------
    print(f"Running PCA (n_comps={n_pcs}) ...")
    adata.raw = adata.copy()
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, random_state=random_state)
    
    # Elbow plot
    sc.pl.pca_variance_ratio(adata, n_pcs=n_pcs, show=False)
    plt.savefig(f"{output_dir}/pca_elbow.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved pca_elbow.png")

    # --------------------------------------------------
    # 6. Neighbors
    # --------------------------------------------------
    print(f"Computing neighbor graph (n_neighbors={n_neighbors}) ...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)

    # --------------------------------------------------
    # 7. UMAP
    # --------------------------------------------------
    print("Computing UMAP embedding ...")
    sc.tl.umap(adata, random_state=random_state)

    # --------------------------------------------------
    # 8. Leiden Clustering
    # --------------------------------------------------
    print(f"Running Leiden clustering (resolution={resolution}) ...")
    sc.tl.leiden(adata, resolution=resolution, random_state=random_state)
    n_clusters = adata.obs["leiden"].nunique()
    print(f"  Found {n_clusters} clusters")

    # --------------------------------------------------
    # 9. Save UMAP plot
    # --------------------------------------------------
    print("Saving UMAP plot ...")
    sc.pl.umap(adata, color=["leiden"], legend_loc="on data", 
               title=f"Leiden clustering (res={resolution})",
               show=False, frameon=False)
    plt.savefig(f"{output_dir}/umap_leiden.png", dpi=200, bbox_inches="tight")
    plt.close()
    print("  Saved umap_leiden.png")

    # Also save UMAP colored by QC metrics
    sc.pl.umap(adata, color=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
               show=False, frameon=False, ncols=3)
    plt.savefig(f"{output_dir}/umap_qc_metrics.png", dpi=200, bbox_inches="tight")
    plt.close()
    print("  Saved umap_qc_metrics.png")

    # --------------------------------------------------
    # 10. Save processed data
    # --------------------------------------------------
    output_h5ad = f"{output_dir}/processed.h5ad"
    print(f"Saving processed data to {output_h5ad} ...")
    adata.write_h5ad(output_h5ad)
    print("Done!")

if __name__ == "__main__":
    main()