#!/usr/bin/env Rscript
# =============================================================================
#  03_normalize.R — RNA and ATAC normalization
#
#  RNA processing (dual normalization as in the original code):
#    1. SCTransform: variance-stabilized normalization → stored in SCT assay
#       (used for integration anchoring)
#    2. LogNormalize on RNA assay (used for downstream DE analysis with MAST)
#    3. FindVariableFeatures(vst, nfeatures=2000) on RNA assay
#
#  ATAC processing (applied to both CellRanger peaks AND MACS2 peaks):
#    1. FindTopFeatures(min.cutoff='q0') — keep all peaks
#    2. RunTFIDF — term frequency–inverse document frequency normalization
#    3. RunSVD — singular value decomposition (LSI) for dimensionality reduction
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(yaml)
  library(future)
})

config <- yaml::read_yaml("/pipeline/config.yaml")
plan("multicore", workers = config$workers)
options(future.globals.maxSize = config$max_ram_gb * 1024^3)
set.seed(config$seed)

cat("=== Step 03: Normalization ===\n")

# =============================================================================
#  RNA normalization
# =============================================================================
# SCTransform performs regularized negative binomial regression to remove
# technical variation while preserving biological heterogeneity.
# We also run standard LogNormalize for DE analysis compatibility.
process_expression <- function(obj) {
  cat("  SCTransform...\n")
  obj <- SCTransform(obj)

  cat("  LogNormalize on RNA assay...\n")
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = config$processing$nfeatures_vst
  )

  return(obj)
}

# =============================================================================
#  ATAC normalization (TF-IDF + LSI)
# =============================================================================
# TF-IDF: Normalizes peak accessibility by cell depth and peak frequency.
# RunSVD computes the LSI (Latent Semantic Indexing) embedding.
# We run this on both the CellRanger ATAC peaks and the MACS2-called peaks.
process_accessibility <- function(obj) {
  cat("  ATAC assay: TF-IDF + SVD...\n")
  DefaultAssay(obj) <- "ATAC"
  obj <- FindTopFeatures(obj, min.cutoff = config$processing$atac_min_cutoff)
  obj <- RunTFIDF(obj)
  obj <- RunSVD(obj)

  cat("  MACS2 assay: TF-IDF + SVD...\n")
  DefaultAssay(obj) <- "macs2"
  obj <- FindTopFeatures(obj, min.cutoff = config$processing$atac_min_cutoff)
  obj <- RunTFIDF(obj)
  obj <- RunSVD(obj)

  return(obj)
}

# =============================================================================
#  Process both samples
# =============================================================================
cat("\n--- Wild-type: Normalization ---\n")
wild <- readRDS("/output/02_wild_peaks.rds")
wild <- process_expression(wild)
wild <- process_accessibility(wild)
saveRDS(wild, "/output/03_wild_normalized.rds")
rm(wild); gc()

cat("\n--- Mutant: Normalization ---\n")
mutant <- readRDS("/output/02_mutant_peaks.rds")
mutant <- process_expression(mutant)
mutant <- process_accessibility(mutant)
saveRDS(mutant, "/output/03_mutant_normalized.rds")
rm(mutant); gc()

cat("\n=== Step 03 complete ===\n")