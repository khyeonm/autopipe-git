#!/usr/bin/env python3
"""Single-cell DEG analysis with optional Squidpy spatial statistics."""

from __future__ import annotations

import argparse
import json
import os
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq


def str_to_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run single-cell DEG and optional Squidpy spatial analysis.")
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--groupby", required=True)
    parser.add_argument("--reference", default="rest")
    parser.add_argument("--method", default="wilcoxon")
    parser.add_argument("--layer", default="")
    parser.add_argument("--use-raw", default="false")
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--normalize-total", default="true")
    parser.add_argument("--target-sum", type=float, default=10000)
    parser.add_argument("--log1p", default="true")
    parser.add_argument("--n-top-genes", type=int, default=25)
    parser.add_argument("--spatial-library-key", default="")
    parser.add_argument("--spatial-coord-type", default="generic")
    parser.add_argument("--spatial-radius", default="")
    parser.add_argument("--spatial-n-neighs", type=int, default=6)
    parser.add_argument("--spatial-autocorr-genes", type=int, default=100)
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args()


def ensure_output_files(output_dir: Path) -> dict[str, Path]:
    paths = {
        "deg": output_dir / "deg_results.tsv",
        "top": output_dir / "top_markers_by_group.tsv",
        "summary": output_dir / "analysis_summary.json",
        "qc": output_dir / "qc_overview.png",
        "dotplot": output_dir / "deg_dotplot.png",
        "nhood": output_dir / "squidpy_neighborhood_enrichment.tsv",
        "spatial_autocorr": output_dir / "squidpy_spatial_autocorr.tsv",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    return paths


def write_empty_table(path: Path, columns: list[str], reason: str) -> None:
    df = pd.DataFrame(columns=columns)
    df.attrs["reason"] = reason
    df.to_csv(path, sep="\t", index=False)


def basic_qc_plot(adata: sc.AnnData, groupby: str, output: Path) -> None:
    obs = adata.obs.copy()
    if "total_counts" not in obs.columns:
        obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
    order = obs[groupby].astype(str).value_counts().index.tolist()
    sns.countplot(data=obs, y=groupby, order=order, ax=axes[0], color="#4c78a8")
    axes[0].set_title("Cells per group")
    axes[0].set_xlabel("Cells")
    axes[0].set_ylabel(groupby)

    sns.histplot(obs["total_counts"], bins=60, ax=axes[1], color="#f58518")
    axes[1].set_title("Library size")
    axes[1].set_xlabel("Total counts")
    axes[1].set_ylabel("Cells")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def preprocess(adata: sc.AnnData, args: argparse.Namespace) -> sc.AnnData:
    adata = adata.copy()
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    if str_to_bool(args.normalize_total):
        sc.pp.normalize_total(adata, target_sum=args.target_sum)
    if str_to_bool(args.log1p):
        sc.pp.log1p(adata)
    return adata


def run_deg(adata: sc.AnnData, args: argparse.Namespace, paths: dict[str, Path]) -> pd.DataFrame:
    rank_kwargs = {
        "groupby": args.groupby,
        "method": args.method,
        "reference": args.reference,
        "use_raw": str_to_bool(args.use_raw),
    }
    if args.layer:
        rank_kwargs["layer"] = args.layer

    sc.tl.rank_genes_groups(adata, **rank_kwargs)
    deg = sc.get.rank_genes_groups_df(adata, group=None)

    expected = ["group", "names", "scores", "logfoldchanges", "pvals", "pvals_adj"]
    for column in expected:
        if column not in deg.columns:
            deg[column] = np.nan

    deg = deg.rename(columns={"names": "gene"})
    deg = deg.sort_values(["group", "pvals_adj", "pvals"], na_position="last")
    deg.to_csv(paths["deg"], sep="\t", index=False)

    top = deg.groupby("group", group_keys=False).head(args.n_top_genes).copy()
    top.to_csv(paths["top"], sep="\t", index=False)
    return top


def marker_dotplot(adata: sc.AnnData, groupby: str, top: pd.DataFrame, n_top: int, output: Path) -> None:
    genes = top["gene"].dropna().astype(str).drop_duplicates().head(max(n_top, 1)).tolist()
    if not genes:
        write_placeholder_png(output, "No marker genes available for dotplot")
        return

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pl.dotplot(
            adata,
            var_names=genes,
            groupby=groupby,
            standard_scale="var",
            show=False,
            save=None,
        )
    fig = plt.gcf()
    fig.set_size_inches(max(8, len(genes) * 0.35), 5)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_placeholder_png(path: Path, message: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def run_squidpy(adata: sc.AnnData, args: argparse.Namespace, top: pd.DataFrame, paths: dict[str, Path]) -> dict[str, str]:
    status = {}
    if "spatial" not in adata.obsm:
        write_empty_table(paths["nhood"], ["group", "neighbor_group", "zscore", "count"], "adata.obsm['spatial'] not found")
        write_empty_table(paths["spatial_autocorr"], ["gene", "I", "pval_norm", "pval_norm_fdr_bh"], "adata.obsm['spatial'] not found")
        status["squidpy"] = "skipped: adata.obsm['spatial'] not found"
        return status

    coord_kwargs = {"coord_type": args.spatial_coord_type, "n_neighs": args.spatial_n_neighs}
    if args.spatial_radius:
        coord_kwargs["radius"] = float(args.spatial_radius)
    if args.spatial_library_key:
        coord_kwargs["library_key"] = args.spatial_library_key

    sq.gr.spatial_neighbors(adata, **coord_kwargs)
    sq.gr.nhood_enrichment(adata, cluster_key=args.groupby)
    zscore = adata.uns[f"{args.groupby}_nhood_enrichment"]["zscore"]
    counts = adata.uns[f"{args.groupby}_nhood_enrichment"]["count"]
    categories = adata.obs[args.groupby].cat.categories.astype(str).tolist()

    rows = []
    for i, group in enumerate(categories):
        for j, neighbor_group in enumerate(categories):
            rows.append({
                "group": group,
                "neighbor_group": neighbor_group,
                "zscore": zscore[i, j],
                "count": counts[i, j],
            })
    pd.DataFrame(rows).to_csv(paths["nhood"], sep="\t", index=False)

    genes = top["gene"].dropna().astype(str).drop_duplicates().head(args.spatial_autocorr_genes).tolist()
    genes = [gene for gene in genes if gene in adata.var_names]
    if genes:
        sq.gr.spatial_autocorr(adata, mode="moran", genes=genes)
        adata.uns["moranI"].reset_index(names="gene").to_csv(paths["spatial_autocorr"], sep="\t", index=False)
    else:
        write_empty_table(paths["spatial_autocorr"], ["gene", "I", "pval_norm", "pval_norm_fdr_bh"], "No DEG genes available in adata.var_names")

    status["squidpy"] = "completed"
    return status


def main() -> None:
    args = parse_args()
    sc.settings.n_jobs = args.threads
    paths = ensure_output_files(Path(args.output_dir))

    adata = sc.read_h5ad(args.input_h5ad)
    if args.groupby not in adata.obs.columns:
        raise ValueError(f"groupby column '{args.groupby}' was not found in adata.obs")

    adata.obs[args.groupby] = adata.obs[args.groupby].astype("category")
    if adata.obs[args.groupby].cat.categories.size < 2:
        raise ValueError(f"groupby column '{args.groupby}' must contain at least two groups")

    basic_qc_plot(adata, args.groupby, paths["qc"])
    adata = preprocess(adata, args)
    top = run_deg(adata, args, paths)
    marker_dotplot(adata, args.groupby, top, args.n_top_genes, paths["dotplot"])
    spatial_status = run_squidpy(adata, args, top, paths)

    summary = {
        "input_h5ad": args.input_h5ad,
        "n_cells": int(adata.n_obs),
        "n_genes_after_filtering": int(adata.n_vars),
        "groupby": args.groupby,
        "groups": adata.obs[args.groupby].cat.categories.astype(str).tolist(),
        "method": args.method,
        "reference": args.reference,
        "use_raw": str_to_bool(args.use_raw),
        "layer": args.layer or None,
        "normalize_total": str_to_bool(args.normalize_total),
        "log1p": str_to_bool(args.log1p),
        "outputs": {key: str(value) for key, value in paths.items()},
        **spatial_status,
    }
    paths["summary"].write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()