#!/usr/bin/env python3
"""
Step 1: Quality Control
- Filter cells by gene count and mitochondrial percentage
- Calculate QC metrics
"""

import scanpy as sc
import numpy as np
import pandas as pd
import logging
import os
import sys

# Get parameters from environment
INPUT_PATH = os.environ.get('INPUT_PATH')
OUTPUT_PATH = os.environ.get('OUTPUT_PATH')
LOG_PATH = os.environ.get('LOG_PATH')
MIN_CELLS = int(os.environ.get('MIN_CELLS', 3))
MIN_GENES = int(os.environ.get('MIN_GENES', 200))
MAX_GENES = int(os.environ.get('MAX_GENES', 7500))
MAX_MITO_PERCENT = float(os.environ.get('MAX_MITO_PERCENT', 20.0))

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
    logger.info(f"Initial data: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Calculate mitochondrial gene percentage
    mito_genes = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
    adata.var['mito'] = mito_genes
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mito'],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    logger.info("QC metrics calculated")
    logger.info(f"  Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
    logger.info(f"  Median counts/cell: {adata.obs['total_counts'].median():.0f}")
    logger.info(f"  Median mito%: {adata.obs['pct_counts_mito'].median():.2f}%")
    
    # Filter cells
    logger.info(f"Filtering: min_genes={MIN_GENES}, max_genes={MAX_GENES}, max_mito%={MAX_MITO_PERCENT}")
    
    initial_cells = adata.n_obs
    initial_genes = adata.n_vars
    
    # Filter based on gene counts
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    sc.pp.filter_cells(adata, max_genes=MAX_GENES)
    
    # Filter based on mitochondrial percentage
    adata = adata[adata.obs['pct_counts_mito'] < MAX_MITO_PERCENT, :]
    
    # Filter genes detected in too few cells
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    
    logger.info(f"After QC: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Save filtered data
    logger.info(f"Saving to {OUTPUT_PATH}")
    adata.write(OUTPUT_PATH)
    logger.info("QC filtering complete")


if __name__ == '__main__':
    main()