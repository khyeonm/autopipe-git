#!/usr/bin/env python3
"""Step 6: Leiden Clustering"""
import scanpy as sc
import logging
import os

INPUT_PATH = os.environ.get('INPUT_PATH')
OUTPUT_PATH = os.environ.get('OUTPUT_PATH')
LOG_PATH = os.environ.get('LOG_PATH')
RESOLUTION = float(os.environ.get('RESOLUTION', 0.5))

os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(LOG_PATH), logging.StreamHandler()])
logger = logging.getLogger(__name__)

def main():
    logger.info(f"Loading data from {INPUT_PATH}")
    adata = sc.read(INPUT_PATH)
    logger.info(f"Data loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    
    logger.info(f"Running Leiden clustering with resolution={RESOLUTION}")
    sc.tl.leiden(adata, resolution=RESOLUTION, key_added='leiden')
    
    n_clusters = adata.obs['leiden'].nunique()
    logger.info(f"Found {n_clusters} clusters")
    logger.info(f"Cluster sizes: {adata.obs['leiden'].value_counts().to_dict()}")
    
    logger.info(f"Saving to {OUTPUT_PATH}")
    adata.write(OUTPUT_PATH)
    logger.info("Clustering complete")

if __name__ == '__main__':
    main()