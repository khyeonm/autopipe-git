#!/usr/bin/env python3
"""Step 7: UMAP Embedding"""
import scanpy as sc
import logging
import os

INPUT_PATH = os.environ.get('INPUT_PATH')
OUTPUT_PATH = os.environ.get('OUTPUT_PATH')
LOG_PATH = os.environ.get('LOG_PATH')
UMAP_MIN_DIST = float(os.environ.get('UMAP_MIN_DIST', 0.3))

os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(LOG_PATH), logging.StreamHandler()])
logger = logging.getLogger(__name__)

def main():
    logger.info(f"Loading data from {INPUT_PATH}")
    adata = sc.read(INPUT_PATH)
    logger.info(f"Data loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    
    logger.info(f"Running UMAP with min_dist={UMAP_MIN_DIST}")
    sc.tl.umap(adata, min_dist=UMAP_MIN_DIST)
    
    logger.info(f"Saving to {OUTPUT_PATH}")
    adata.write(OUTPUT_PATH)
    logger.info("UMAP embedding complete")

if __name__ == '__main__':
    main()