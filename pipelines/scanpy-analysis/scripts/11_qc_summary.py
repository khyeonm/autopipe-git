#!/usr/bin/env python3
"""Step 11: Generate QC Summary Report"""

import scanpy as sc
import pandas as pd
import logging
import os

# Get parameters from environment
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
    
    # Generate summary statistics
    summary = {
        'n_cells': [adata.n_obs],
        'n_genes': [adata.n_vars],
        'median_genes_per_cell': [adata.obs['n_genes_by_counts'].median()],
        'median_counts_per_cell': [adata.obs['total_counts'].median()],
        'mean_mito_percent': [adata.obs['pct_counts_mito'].mean()],
        'median_mito_percent': [adata.obs['pct_counts_mito'].median()],
        'max_mito_percent': [adata.obs['pct_counts_mito'].max()],
    }
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(OUTPUT_PATH, index=False)
    
    logger.info(f"QC summary saved to {OUTPUT_PATH}")
    logger.info(f"Summary:\n{summary_df}")

if __name__ == '__main__':
    main()