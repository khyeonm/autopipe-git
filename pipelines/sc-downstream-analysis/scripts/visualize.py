#!/usr/bin/env python3
"""Step 3: UMAP visualization."""

import argparse
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--output-png", required=True)
    parser.add_argument("--min-dist", type=float, default=0.3)
    parser.add_argument("--spread", type=float, default=1.0)
    args = parser.parse_args()

    print(f"[visualize] Loading {args.input}")
    adata = sc.read_h5ad(args.input)

    # Compute UMAP
    sc.tl.umap(adata, min_dist=args.min_dist, spread=args.spread)
    print("[visualize] UMAP computed")

    # Plot
    n_clusters = adata.obs["leiden"].nunique()
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    sc.pl.umap(adata, color="leiden", title=f"Leiden Clusters (n={n_clusters})",
               legend_loc="on data", legend_fontsize=8, ax=axes[0], show=False)

    sc.pl.umap(adata, color="leiden", title="Leiden Clusters (legend)",
               legend_loc="right margin", ax=axes[1], show=False)

    plt.tight_layout()
    plt.savefig(args.output_png, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[visualize] UMAP PNG saved to {args.output_png}")

    adata.write_h5ad(args.output_h5ad)
    print(f"[visualize] Final h5ad saved to {args.output_h5ad}")

if __name__ == "__main__":
    main()