#!/usr/bin/env python3
"""Step 5: Build KNN Graph"""
import scanpy as sc
import logging
import os

INPUT_PATH = os.environ.get('INPUT_PATH')
OUTPUT_PATH = os.environ.get('OUTPUT_PATH')
LOG_PATH = os.environ.get('LOG_PATH')
N_NEIGHBORS = int(os.environ.get('N_NEIGHBORS', 15))
N_PCS = int(os.environ.get('N_PCS', 50))

os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(LOG_PATH), logging.StreamHandler()])
logger = logging.getLogger(__name__)

def main():
    logger.info(f"Loading data from {INPUT_PATH}")
    adata = sc.read(INPUT_PATH)
    logger.info(f"Data loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    
    logger.info(f"Building KNN graph with n_neighbors={N_NEIGHBORS}, n_pcs={N_PCS}")
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)
    
    logger.info(f"Saving to {OUTPUT_PATH}")
    adata.write(OUTPUT_PATH)
    logger.info("KNN graph construction complete")

if __name__ == '__main__':
    main()