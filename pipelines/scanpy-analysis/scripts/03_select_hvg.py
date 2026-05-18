#!/usr/bin/env python3
"""Step 3: Select Highly Variable Genes"""
import scanpy as sc
import logging
import os

INPUT_PATH = os.environ.get('INPUT_PATH')
OUTPUT_PATH = os.environ.get('OUTPUT_PATH')
LOG_PATH = os.environ.get('LOG_PATH')
N_HVG = int(os.environ.get('N_HVG', 2000))

os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(LOG_PATH), logging.StreamHandler()])
logger = logging.getLogger(__name__)

def main():
    logger.info(f"Loading data from {INPUT_PATH}")
    adata = sc.read(INPUT_PATH)
    logger.info(f"Data loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    
    logger.info(f"Selecting top {N_HVG} highly variable genes (flavor='seurat')")
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor='seurat', subset=True, inplace=True)
    
    logger.info(f"After HVG selection: {adata.n_obs} cells, {adata.n_vars} genes")
    logger.info(f"Saving to {OUTPUT_PATH}")
    adata.write(OUTPUT_PATH)
    logger.info("HVG selection complete")

if __name__ == '__main__':
    main()