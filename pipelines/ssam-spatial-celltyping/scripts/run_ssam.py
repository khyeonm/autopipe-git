#!/usr/bin/env python3
"""
SSAM pipeline script — uses pnucolab/ssam API (github.com/pnucolab/ssam)
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Priority-ordered candidate column names for gene, x, y
GENE_CANDIDATES = ["gene", "target", "gene_name", "genename", "feature_name"]
X_CANDIDATES    = ["x", "xcoord", "x_coord", "x_location", "xc", "rotated_x", "x_um"]
Y_CANDIDATES    = ["y", "ycoord", "y_coord", "y_location", "yc", "rotated_y", "y_um"]


def pick_column(df_cols, candidates):
    """Return the first candidate (case-insensitive) found in df_cols, or None."""
    lower_map = {c.lower(): c for c in df_cols}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def load_from_csv(coord_path):
    print(f"[SSAM] Loading spot table from {coord_path} ...")
    df = pd.read_csv(coord_path)
    print(f"[SSAM] Columns: {list(df.columns)}")
    print(f"[SSAM] Shape: {df.shape}")

    gene_col = pick_column(df.columns, GENE_CANDIDATES)
    x_col    = pick_column(df.columns, X_CANDIDATES)
    y_col    = pick_column(df.columns, Y_CANDIDATES)

    missing = [name for name, col in [("gene", gene_col), ("x", x_col), ("y", y_col)] if col is None]
    if missing:
        sys.exit(
            f"[SSAM] ERROR: Could not find columns for {missing}.\n"
            f"  Available columns: {list(df.columns)}\n"
            f"  Tried gene={GENE_CANDIDATES}, x={X_CANDIDATES}, y={Y_CANDIDATES}"
        )

    print(f"[SSAM] Using columns — gene='{gene_col}', x='{x_col}', y='{y_col}'")

    df = df[[gene_col, x_col, y_col]].copy()
    df.columns = ["gene", "x", "y"]

    df["x"] = pd.to_numeric(df["x"], errors="coerce")
    df["y"] = pd.to_numeric(df["y"], errors="coerce")
    df = df.dropna()

    print(f"[SSAM] Transcripts: {len(df):,} rows, {df['gene'].nunique()} genes")
    return df


def find_localmax_safe(analysis, threshold):
    """
    Call find_localmax with the correct keyword argument.
    The pnucolab/ssam API accepts either 'min_norm' or 'norm_thres'
    depending on the version. We try both.
    """
    import inspect
    sig = inspect.signature(analysis.find_localmax)
    params = list(sig.parameters.keys())
    print(f"[SSAM] find_localmax signature params: {params}")

    if "min_norm" in params:
        analysis.find_localmax(min_norm=threshold)
    elif "norm_thres" in params:
        analysis.find_localmax(norm_thres=threshold)
    elif "threshold" in params:
        analysis.find_localmax(threshold=threshold)
    else:
        print(f"[SSAM] Warning: unknown find_localmax params {params}, trying positional.")
        analysis.find_localmax(threshold)


def check_store_ready(store_path):
    """
    Check which stages are already completed in ssam_store.
    Returns a dict of {stage: bool}.
    """
    def has_data(subdir):
        p = os.path.join(store_path, subdir)
        if not os.path.isdir(p):
            return False
        files = [f for f in os.listdir(p) if not f.startswith('.')]
        return len(files) > 0

    return {
        "kde":        has_data("vf") and has_data("kde_computed"),
        "localmax":   has_data("local_maxs"),
        "normalized": has_data("normalized_vectors"),
        "scaled":     has_data("scaled_vectors"),
    }


def main():
    parser = argparse.ArgumentParser(description="SSAM spatial analysis (pnucolab/ssam)")
    parser.add_argument("--csv",         required=True,  help="Spot table CSV")
    parser.add_argument("--output-dir",  required=True)
    parser.add_argument("--bandwidth",   type=float, default=2.5)
    parser.add_argument("--threshold",   type=float, default=0.5)
    parser.add_argument("--map-width",   type=int,   default=1000)
    parser.add_argument("--threads",     type=int,   default=4)
    parser.add_argument("--resolution",  type=float, default=0.6)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    store_path = os.path.join(args.output_dir, "ssam_store")
    ready = check_store_ready(store_path)

    print(f"[SSAM] Store status — KDE: {ready['kde']}, LocalMax: {ready['localmax']}, "
          f"Normalized: {ready['normalized']}, Scaled: {ready['scaled']}")

    # ── 1. Load data (needed for KDE; skip if KDE already done) ──────────────
    import ssam
    ds       = ssam.SSAMDataset(store_path)
    analysis = ssam.SSAMAnalysis(ds, ncores=args.threads, verbose=True)

    if not ready["kde"]:
        df = load_from_csv(args.csv)

        um_per_pixel = 0.1
        df["x"] = (df["x"] - df["x"].min()) * um_per_pixel + 10
        df["y"] = (df["y"] - df["y"].min()) * um_per_pixel + 10

        width  = df["x"].max() - df["x"].min() + 10
        height = df["y"].max() - df["y"].min() + 10
        sampling_distance = 1.0

        print(f"[SSAM] Tissue size (µm): {width:.1f} × {height:.1f}")
        print(f"[SSAM] Map: {args.map_width} px wide, sampling={sampling_distance:.4f} µm/px")

        df_indexed = df.set_index("gene")

        print(f"[SSAM] Initialising SSAMDataset at {store_path} ...")
        print("[SSAM] Running KDE ...")
        analysis.run_kde(
            df_indexed,
            width=width,
            height=height,
            bandwidth=args.bandwidth,
            sampling_distance=sampling_distance,
        )
    else:
        print("[SSAM] ✓ KDE already done — skipping.")

    # ── 4. Find local maxima ──────────────────────────────────────────────────
    if not ready["localmax"]:
        print(f"[SSAM] Finding local maxima (threshold={args.threshold}) ...")
        find_localmax_safe(analysis, args.threshold)
    else:
        print(f"[SSAM] ✓ Local maxima already found — skipping. (threshold={args.threshold})")

    # ── 5. Normalise & scale vectors ──────────────────────────────────────────
    if not ready["normalized"]:
        print("[SSAM] Normalising vectors ...")
        analysis.normalize_vectors()
    else:
        print("[SSAM] ✓ Normalization already done — skipping.")

    if not ready["scaled"]:
        print("[SSAM] Scaling vectors ...")
        analysis.scale_vectors()
    else:
        print("[SSAM] ✓ Scaling already done — skipping.")

    # ── 6. Cluster & map cell types (unsupervised) ────────────────────────────
    print(f"[SSAM] Clustering vectors (resolution={args.resolution}) ...")
    analysis.cluster_vectors(resolution=args.resolution, metric="correlation")

    print("[SSAM] Mapping cell types ...")
    analysis.map_celltypes()

    # ── 7. Save KDE map ───────────────────────────────────────────────────────
    print("[SSAM] Saving KDE map ...")
    vf_norm = np.array(ds.vf_norm).squeeze()
    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(np.log1p(vf_norm).T, cmap="hot", origin="lower")
    plt.colorbar(im, ax=ax, label="log1p(KDE norm)")
    ax.set_title("KDE Density Map")
    ax.set_xlabel("x (pixels)")
    ax.set_ylabel("y (pixels)")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "kde_map.png"), dpi=150)
    plt.close()

    # ── 8. Save cell-type map ─────────────────────────────────────────────────
    print("[SSAM] Saving cell-type map ...")
    ct_map = np.array(ds.celltype_maps).squeeze()
    n_types = int(ct_map.max()) + 1
    cmap = plt.get_cmap("tab20", n_types + 1)
    fig, ax = plt.subplots(figsize=(12, 12))
    im = ax.imshow(ct_map.T, cmap=cmap, vmin=-0.5, vmax=n_types - 0.5, origin="lower")
    cbar = plt.colorbar(im, ax=ax, ticks=range(n_types))
    cbar.ax.set_yticklabels([f"Cluster {i}" for i in range(n_types)], fontsize=7)
    ax.set_title("Cell-type Map (SSAM)")
    ax.set_xlabel("x (pixels)")
    ax.set_ylabel("y (pixels)")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "celltype_map.png"), dpi=150)
    plt.close()

    # ── 9. Abundance table ────────────────────────────────────────────────────
    print("[SSAM] Computing abundance ...")
    valid  = ct_map[ct_map >= 0]
    counts = np.bincount(valid, minlength=n_types)
    total  = counts.sum()
    abund  = pd.DataFrame({
        "cluster":     [f"Cluster {i}" for i in range(n_types)],
        "pixel_count": counts,
        "fraction":    counts / total if total > 0 else counts,
    })
    abund.to_csv(os.path.join(args.output_dir, "celltype_abundance.csv"), index=False)
    print(abund.to_string(index=False))
    print(f"\n[SSAM] Done! Results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()

