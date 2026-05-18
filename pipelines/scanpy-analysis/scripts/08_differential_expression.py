#!/usr/bin/env python3
"""Step 8: Differential Expression Analysis"""
import scanpy as sc
import pandas as pd
import logging
import os

INPUT_PATH = os.environ.get('INPUT_PATH')
OUTPUT_PATH = os.environ.get('OUTPUT_PATH')
LOG_PATH = os.environ.get('LOG_PATH')
N_TOP = int(os.environ.get('DE_N_TOP', 50))
MIN_FRAC = float(os.environ.get('DE_MIN_FRAC', 0.25))
LOGFC_THRESH = float(os.environ.get('DE_LOGFC_THRESH', 0.25))

os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(LOG_PATH), logging.StreamHandler()])
logger = logging.getLogger(__name__)

def main():
    logger.info(f"Loading data from {INPUT_PATH}")
    adata = sc.read(INPUT_PATH)
    logger.info(f"Data loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    
    logger.info(f"Running DE analysis: n_top={N_TOP}, min_frac={MIN_FRAC}, logfc_thresh={LOGFC_THRESH}")
    sc.tl.rank_genes_groups(
        adata, groupby='leiden', method='wilcoxon', key_added='rank_genes',
        n_genes=N_TOP + 10, use_raw=False
    )
    
    de_results = sc.get.rank_genes_groups_df(adata, group=None, key='rank_genes')
    logger.info(f"DE results: {len(de_results)} total comparisons")
    
    logger.info(f"Saving to {OUTPUT_PATH}")
    adata.write(OUTPUT_PATH)
    logger.info("Differential expression analysis complete")

if __name__ == '__main__':
    main()