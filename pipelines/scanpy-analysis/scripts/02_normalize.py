#!/usr/bin/env python3
"""
Step 2: Normalization - Log-normalize the count data
"""

import scanpy as sc
import logging
import os

# Get parameters from environment
INPUT_PATH = os.environ.get('INPUT_PATH')
OUTPUT_PATH = os.environ.get('OUTPUT_PATH')
LOG_PATH = os.environ.get('LOG_PATH')
NORMALIZE_TARGET = int(os.environ.get('NORMALIZE_TARGET', 10000))

# Setup logging
os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_PATH),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def main():
    logger.info(f"Loading data from {INPUT_PATH}")
    adata = sc.read(INPUT_PATH)
    logger.info(f"Data loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Normalize the data
    logger.info(f"Normalizing with target_sum={NORMALIZE_TARGET}")
    sc.pp.normalize_total(adata, target_sum=NORMALIZE_TARGET)
    sc.pp.log1p(adata)
    
    logger.info(f"Saving normalized data to {OUTPUT_PATH}")
    adata.write(OUTPUT_PATH)
    logger.info("Normalization complete")


if __name__ == '__main__':
    main()