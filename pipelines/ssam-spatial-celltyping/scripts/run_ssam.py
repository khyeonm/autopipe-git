#!/usr/bin/env python3
"""
SSAM pipeline script — follows the official pnucolab/ssam workflow
(github.com/pnucolab/ssam, docs/tldr.rst).
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


def check_store_ready(store_path):
    """Check which stages are already completed in ssam_store."""
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
    parser.add_argument("--bandwidth",   type=float, default=2.5,
                        help="KDE bandwidth in micrometers (ssam default: 2.5).")
    parser.add_argument("--sampling-distance", type=float, default=1.0,
                        help="KDE grid spacing in micrometers/pixel (ssam default: 1.0).")
    parser.add_argument("--norm-threshold", type=float, default=None,
                        help="Optional norm-threshold override for local-maxima detection. "
                             "If omitted, SSAM's canonical bandwidth-derived default is used "
                             "(norm_threshold = 2 / (sqrt(2*pi)*bandwidth)**ndim).")
    parser.add_argument("--search-size", type=int, default=3,
                        help="find_localmax search neighborhood size in pixels (ssam default: 3).")
    parser.add_argument("--resolution",  type=float, default=0.6,
                        help="Leiden clustering resolution (canonical: 0.6).")
    parser.add_argument("--min-norm",    type=float, default=0.1,
                        help="filter_celltypemaps min_norm (ssam default: 0.1).")
    parser.add_argument("--outlier-min-r", type=float, default=0.3,
                        help="cluster_vectors MedoidCorrelation min_r for outlier removal "
                             "(ssam default 0.8; negative disables outlier removal).")
    parser.add_argument("--filter-min-r", type=float, default=0.3,
                        help="filter_celltypemaps min_r (ssam default 0.6).")
    parser.add_argument("--no-scale", action="store_true",
                        help="Skip scale_vectors; reuse the normalized vectors in the scaled slots "
                             "(scaled_vectors := normalized_vectors, vf_scaled := vf_normalized). "
                             "Deviates from canonical SSAM: no per-gene standardization.")
    parser.add_argument("--threads",     type=int, default=4)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    store_path = os.path.join(args.output_dir, "ssam_store")
    ready = check_store_ready(store_path)

    print(f"[SSAM] Store status — KDE: {ready['kde']}, LocalMax: {ready['localmax']}, "
          f"Normalized: {ready['normalized']}, Scaled: {ready['scaled']}")

    import ssam
    ds       = ssam.SSAMDataset(store_path)
    analysis = ssam.SSAMAnalysis(ds, ncores=args.threads, verbose=True)

    # ── 1. KDE (real micrometer coordinates; grid spacing = sampling_distance) ──
    if not ready["kde"]:
        df = load_from_csv(args.csv)

        # Canonical ssam: normalize coordinates and keep the REAL scale (bandwidth
        # is in the same units as x/y), tissue size = max coordinate.
        df["x"] = df["x"] - df["x"].min()
        df["y"] = df["y"] - df["y"].min()
        width  = float(df["x"].max())
        height = float(df["y"].max())
        sampling_distance = args.sampling_distance

        print(f"[SSAM] Tissue size (µm): {width:.1f} × {height:.1f}")
        print(f"[SSAM] sampling={sampling_distance} µm/px, bandwidth={args.bandwidth} µm "
              f"-> grid ~{int(width/sampling_distance)} x {int(height/sampling_distance)} px")

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

    # ── 2. Find local maxima (canonical bandwidth-derived threshold; optional override) ──
    if not ready["localmax"]:
        # Canonical ssam relies on the bandwidth-derived thresholds that run_kde
        # sets (norm_threshold = 2 / (sqrt(2*pi)*bandwidth)**ndim). Only override
        # when the user explicitly provides one. NOTE: in ssam 1.1.0 find_localmax's
        # only positional arg is `search_size`, NOT a threshold, so any threshold
        # MUST go through set_thresholds() -- never positionally into find_localmax().
        if args.norm_threshold is not None:
            analysis.set_thresholds(norm_threshold=args.norm_threshold)
            print(f"[SSAM] Finding local maxima (norm_threshold override={args.norm_threshold}) ...")
        else:
            print("[SSAM] Finding local maxima (canonical bandwidth-derived threshold) ...")
        analysis.find_localmax(search_size=args.search_size)
    else:
        print(f"[SSAM] ✓ Local maxima already found — skipping.")

    # ── 3. Normalise & scale vectors ─────────────────────────────────────────────
    if not ready["normalized"]:
        print("[SSAM] Normalising vectors ...")
        analysis.normalize_vectors()
    else:
        print("[SSAM] ✓ Normalization already done — skipping.")

    if not ready["scaled"]:
        if args.no_scale:
            # Skip scale_vectors and let the downstream cluster/map steps read the
            # already-computed normalized vectors instead (they read scaled_vectors /
            # vf_scaled). No per-gene standardization is applied.
            print("[SSAM] --no-scale: reusing normalized vectors as scaled (skipping scale_vectors) ...")
            ds.scaled_vectors = ds.normalized_vectors
            ds.vf_scaled = ds.vf_normalized
        else:
            print("[SSAM] Scaling vectors ...")
            analysis.scale_vectors()
    else:
        print("[SSAM] ✓ Scaling already done — skipping.")

    # ── 4. Cluster & map cell types (de novo), then filter (official workflow) ───
    if args.outlier_min_r < 0:
        print(f"[SSAM] Clustering vectors (resolution={args.resolution}, outlier removal disabled) ...")
        analysis.cluster_vectors(resolution=args.resolution, metric="correlation",
                                 outlier_detection_method=None)
    else:
        print(f"[SSAM] Clustering vectors (resolution={args.resolution}, "
              f"outlier min_r={args.outlier_min_r}) ...")
        analysis.cluster_vectors(resolution=args.resolution, metric="correlation",
                                 outlier_detection_kwargs={"min_r": args.outlier_min_r})

    if len(np.asarray(ds.centroids)) == 0:
        sys.exit(
            "[SSAM] ERROR: clustering produced 0 clusters, so every downstream output would be "
            "empty. This normally means the outlier filter removed all vectors -- lower "
            "--outlier-min-r (currently %s), or set it negative to disable outlier removal."
            % args.outlier_min_r
        )

    print("[SSAM] Mapping cell types ...")
    analysis.map_celltypes()

    print(f"[SSAM] Filtering cell-type map (min_norm={args.min_norm}, min_r={args.filter_min_r}) ...")
    analysis.filter_celltypemaps(min_norm=args.min_norm, min_r=args.filter_min_r)

    # ── 5. Save KDE map ──────────────────────────────────────────────────────────
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

    # ── 6. Save cell-type map (use the FILTERED map) ─────────────────────────────
    print("[SSAM] Saving cell-type map ...")
    ct_map = np.array(ds.filtered_celltype_maps).squeeze()
    n_types = int(ct_map.max()) + 1
    # Mask background/unassigned pixels (-1) so they render neutral, not as cluster 0.
    ct_masked = np.ma.masked_less(ct_map, 0)
    cmap = plt.get_cmap("tab20", n_types).copy()
    cmap.set_bad("#ffffff")
    fig, ax = plt.subplots(figsize=(12, 12))
    im = ax.imshow(ct_masked.T, cmap=cmap, vmin=-0.5, vmax=n_types - 0.5, origin="lower")
    cbar = plt.colorbar(im, ax=ax, ticks=range(n_types))
    cbar.ax.set_yticklabels([f"Cluster {i}" for i in range(n_types)], fontsize=7)
    ax.set_title("Cell-type Map (SSAM)")
    ax.set_xlabel("x (pixels)")
    ax.set_ylabel("y (pixels)")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "celltype_map.png"), dpi=150)
    plt.close()

    # ── 7. Abundance table ───────────────────────────────────────────────────────
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
