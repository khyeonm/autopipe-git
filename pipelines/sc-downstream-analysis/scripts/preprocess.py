#!/usr/bin/env python3
"""Step 1: Preprocessing — normalize, log-transform, HVG selection, scaling, PCA."""

import argparse
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--target-sum", type=float, default=1e4)
    args = parser.parse_args()

    print(f"[preprocess] Loading {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"[preprocess] Shape: {adata.shape}")

    # Basic QC filter
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)

    # Highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes, subset=True)
    print(f"[preprocess] HVGs selected: {adata.shape[1]}")

    # Scale
    sc.pp.scale(adata, max_value=10)

    # PCA
    n_pcs = min(args.n_pcs, adata.shape[0] - 1, adata.shape[1] - 1)
    sc.tl.pca(adata, n_comps=n_pcs)
    print(f"[preprocess] PCA done with {n_pcs} components")

    adata.write_h5ad(args.output)
    print(f"[preprocess] Saved to {args.output}")

if __name__ == "__main__":
    main()