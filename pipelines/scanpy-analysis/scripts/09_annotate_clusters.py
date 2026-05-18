#!/usr/bin/env python3
"""Step 9: Cell Type Annotation and Summary"""
import scanpy as sc
import pandas as pd
import numpy as np
import logging
import os

INPUT_PATH = os.environ.get('INPUT_PATH')
OUTPUT_PATH = os.environ.get('OUTPUT_PATH')
LOG_PATH = os.environ.get('LOG_PATH')

os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(LOG_PATH), logging.StreamHandler()])
logger = logging.getLogger(__name__)

def main():
    logger.info(f"Loading data from {INPUT_PATH}")
    adata = sc.read(INPUT_PATH)
    logger.info(f"Data loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Get cluster statistics
    cluster_stats = adata.obs['leiden'].value_counts().reset_index()
    cluster_stats.columns = ['cluster', 'n_cells']
    cluster_stats['pct_cells'] = (cluster_stats['n_cells'] / adata.n_obs * 100).round(2)
    
    logger.info("Cluster statistics:")
    logger.info(cluster_stats.to_string(index=False))
    
    # Get top marker genes for each cluster
    de_results = sc.get.rank_genes_groups_df(adata, group=None, key='rank_genes')
    logger.info(f"DE results columns: {de_results.columns.tolist()}")
    
    # Find the cluster column name (could be 'groups' or 'names' depending on scanpy version)
    cluster_col = 'groups' if 'groups' in de_results.columns else (
        'names' if 'names' in de_results.columns else de_results.columns[0]
    )
    logger.info(f"Using cluster column: {cluster_col}")
    
    # Create annotation summary
    annotation_summary = []
    for cluster in sorted(adata.obs['leiden'].unique()):
        cluster_markers = de_results[de_results[cluster_col] == cluster].head(10)
        top_markers = ', '.join(cluster_markers['names'].values)
        annotation_summary.append({
            'cluster': cluster,
            'n_cells': int(cluster_stats[cluster_stats['cluster'] == cluster]['n_cells'].values[0]),
            'pct_cells': float(cluster_stats[cluster_stats['cluster'] == cluster]['pct_cells'].values[0]),
            'top_markers': top_markers
        })
    
    annotation_df = pd.DataFrame(annotation_summary)
    logger.info("\nTop 10 markers per cluster:")
    logger.info(annotation_df.to_string(index=False))
    
    # Store cluster info in adata
    adata.obs['cluster_size'] = adata.obs['leiden'].map(
        dict(zip(cluster_stats['cluster'], cluster_stats['n_cells']))
    )
    
    # Save data
    logger.info(f"Saving to {OUTPUT_PATH}")
    adata.write(OUTPUT_PATH)
    logger.info("Annotation complete")

if __name__ == '__main__':
    main()