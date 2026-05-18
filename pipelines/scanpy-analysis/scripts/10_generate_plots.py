#!/usr/bin/env python3
"""Step 10: Generate Visualization Plots"""

import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
import logging

# Get parameters from environment
INPUT_PATH = os.environ.get('INPUT_PATH')
UMAP_PNG = os.environ.get('UMAP_PNG')
UMAP_PDF = os.environ.get('UMAP_PDF')
HEATMAP_PNG = os.environ.get('HEATMAP_PNG')
HEATMAP_PDF = os.environ.get('HEATMAP_PDF')
PCA_PNG = os.environ.get('PCA_PNG')
PCA_PDF = os.environ.get('PCA_PDF')
CLUSTER_STATS = os.environ.get('CLUSTER_STATS')
LOG_PATH = os.environ.get('LOG_PATH')

os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
for f in [UMAP_PNG, UMAP_PDF, HEATMAP_PNG, HEATMAP_PDF, PCA_PNG, PCA_PDF]:
    os.makedirs(os.path.dirname(f), exist_ok=True)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(LOG_PATH), logging.StreamHandler()])
logger = logging.getLogger(__name__)

def main():
    logger.info(f"Loading data from {INPUT_PATH}")
    adata = sc.read(INPUT_PATH)
    logger.info(f"Data loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    
    sns.set_style("whitegrid")
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['savefig.dpi'] = 150
    
    # 1. UMAP plot colored by clusters
    logger.info("Generating UMAP plot")
    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.umap(adata, color='leiden', ax=ax, title='UMAP - Cell Clusters',
               show=False, palette='tab20', size=50)
    plt.tight_layout()
    plt.savefig(UMAP_PNG, bbox_inches='tight')
    plt.savefig(UMAP_PDF, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved UMAP plot to {UMAP_PNG}")
    
    # 2. Marker gene heatmap
    logger.info("Generating marker gene heatmap")
    de_results = sc.get.rank_genes_groups_df(adata, group=None, key='rank_genes')
    
    # Find the cluster column name (could be 'group', 'groups', or 'names')
    cluster_col = None
    for col in ['group', 'groups', 'names']:
        if col in de_results.columns:
            cluster_col = col
            break
    
    if cluster_col is None:
        cluster_col = de_results.columns[0]
    
    logger.info(f"Using DE cluster column: {cluster_col}")
    
    clusters = sorted(adata.obs['leiden'].unique(), key=lambda x: int(x) if str(x).isdigit() else x)
    marker_genes = []
    
    for cluster in clusters:
        cluster_markers = de_results[de_results[cluster_col] == cluster].head(10)
        marker_genes.extend(cluster_markers['names'].values)
    
    seen = set()
    unique_genes = []
    for gene in marker_genes:
        if gene not in seen:
            seen.add(gene)
            unique_genes.append(gene)
    
    if len(unique_genes) > 0 and all(gene in adata.var_names for gene in unique_genes):
        adata_marker = adata[:, unique_genes].copy()
        n_cells_display = min(2000, adata.n_obs)
        if adata.n_obs > n_cells_display:
            rng = np.random.default_rng(42)
            idx = rng.choice(adata.n_obs, n_cells_display, replace=False)
            adata_marker = adata_marker[idx, :]
        
        sc.pp.scale(adata_marker, max_value=5)
        df = adata_marker.to_df()
        
        fig, ax = plt.subplots(figsize=(15, max(10, len(unique_genes) * 0.3)))
        sns.heatmap(df.T, cmap='viridis', cbar_kws={'label': 'Expression (scaled)'},
                    yticklabels=True, xticklabels=True, ax=ax)
        ax.set_xlabel('Marker Genes')
        ax.set_ylabel('Cells')
        ax.set_title('Marker Gene Expression Heatmap')
        plt.tight_layout()
        plt.savefig(HEATMAP_PNG, bbox_inches='tight', dpi=150)
        plt.savefig(HEATMAP_PDF, bbox_inches='tight')
        plt.close()
    logger.info(f"Saved heatmap to {HEATMAP_PNG}")
    
    # 3. PCA variance explained plot
    logger.info("Generating PCA variance plot")
    fig, ax1 = plt.subplots(figsize=(10, 6))
    n_pcs = adata.obsm['X_pca'].shape[1]
    variance_ratio = adata.uns['pca']['variance_ratio']
    
    ax1.bar(range(1, n_pcs + 1), variance_ratio, alpha=0.7, color='steelblue')
    ax1.set_xlabel('Principal Component', fontsize=12)
    ax1.set_ylabel('Variance Explained', fontsize=12)
    ax1.set_title('PCA: Variance Explained by Each Component', fontsize=14)
    
    cumsum = np.cumsum(variance_ratio)
    ax2 = ax1.twinx()
    ax2.plot(range(1, n_pcs + 1), cumsum, 'r-', linewidth=2, label='Cumulative')
    ax2.set_ylabel('Cumulative Variance', color='red', fontsize=12)
    ax2.tick_params(axis='y', labelcolor='red')
    
    plt.tight_layout()
    plt.savefig(PCA_PNG, bbox_inches='tight')
    plt.savefig(PCA_PDF, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved PCA plot to {PCA_PNG}")
    
    # 4. Cluster statistics table
    logger.info("Generating cluster statistics")
    cluster_stats_df = adata.obs['leiden'].value_counts().reset_index()
    cluster_stats_df.columns = ['cluster', 'n_cells']
    cluster_stats_df['pct_cells'] = (cluster_stats_df['n_cells'] / adata.n_obs * 100).round(2)
    cluster_stats_df = cluster_stats_df.sort_values('cluster', key=lambda x: pd.to_numeric(x))
    cluster_stats_df['median_genes'] = adata.obs.groupby('leiden')['n_genes_by_counts'].median().values
    cluster_stats_df['median_counts'] = adata.obs.groupby('leiden')['total_counts'].median().values
    
    cluster_stats_df.to_csv(CLUSTER_STATS, index=False)
    logger.info(f"Saved cluster statistics to {CLUSTER_STATS}")
    
    logger.info("All plots generated successfully")

if __name__ == '__main__':
    main()