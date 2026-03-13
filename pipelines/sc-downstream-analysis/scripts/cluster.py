#!/usr/bin/env python3
"""Step 2: Leiden clustering."""

import argparse
import scanpy as sc
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--resolution", type=float, default=0.5)
    parser.add_argument("--n-neighbors", type=int, default=15)
    args = parser.parse_args()

    print(f"[cluster] Loading {args.input}")
    adata = sc.read_h5ad(args.input)

    # Neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep="X_pca")

    # Leiden clustering
    sc.tl.leiden(adata, resolution=args.resolution)
    n_clusters = adata.obs["leiden"].nunique()
    print(f"[cluster] Found {n_clusters} Leiden clusters (resolution={args.resolution})")

    # Save cluster assignments
    cluster_df = pd.DataFrame({
        "cell_barcode": adata.obs_names,
        "leiden_cluster": adata.obs["leiden"].values
    })
    cluster_df.to_csv(args.output_csv, index=False)
    print(f"[cluster] Cluster CSV saved to {args.output_csv}")

    adata.write_h5ad(args.output_h5ad)
    print(f"[cluster] h5ad saved to {args.output_h5ad}")

if __name__ == "__main__":
    main()