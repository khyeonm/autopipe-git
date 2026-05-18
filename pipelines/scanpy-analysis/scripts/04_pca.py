#!/usr/bin/env python3
"""Step 4: Principal Component Analysis"""
import scanpy as sc
import logging
import os

INPUT_PATH = os.environ.get('INPUT_PATH')
OUTPUT_PATH = os.environ.get('OUTPUT_PATH')
LOG_PATH = os.environ.get('LOG_PATH')
N_PCS = int(os.environ.get('N_PCS', 50))

os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(LOG_PATH), logging.StreamHandler()])
logger = logging.getLogger(__name__)

def main():
    logger.info(f"Loading data from {INPUT_PATH}")
    adata = sc.read(INPUT_PATH)
    logger.info(f"Data loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    
    logger.info("Scaling data")
    sc.pp.scale(adata, max_value=10)
    
    logger.info(f"Running PCA with n_comps={N_PCS}")
    sc.tl.pca(adata, n_comps=N_PCS, svd_solver='arpack')
    
    logger.info(f"Saving to {OUTPUT_PATH}")
    adata.write(OUTPUT_PATH)
    logger.info("PCA complete")

if __name__ == '__main__':
    main()