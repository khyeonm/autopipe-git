#!/usr/bin/env Rscript
# ===========================================================================
#  03_normalize.R — RNA + ATAC normalization & dimensional reduction
#  Replicates: process_expression() + process_accessibility() from original
#
#  RNA processing (process_expression):
#    1. SCTransform (normalizes, scales, finds variable features → stored in SCT)
#    2. NormalizeData on RNA assay (LogNormalize for DE analysis)
#    3. FindVariableFeatures (vst, nfeatures=2000) on RNA assay
#
#  ATAC processing (process_accessibility):
#    For both "ATAC" and "macs2" assays:
#    1. FindTopFeatures (min.cutoff = 'q0' — use all features)
#    2. RunTFIDF
#    3. RunSVD (LSI dimensional reduction)
# ===========================================================================

suppressPackageStartupMessages({
    library(yaml)
    library(future)
    library(Seurat)
    library(Signac)
})

# ── Load configuration ────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
config_path <- args[which(args == "--config") + 1]
cfg <- read_yaml(config_path)

set.seed(cfg$seed)
plan("multicore", workers = cfg$workers)
options(future.globals.maxSize = cfg$future_ram_gb * 1024^3)
output_dir <- cfg$output_dir

cat("=== Step 3: Normalize RNA & ATAC ===\n")

# ── Load MACS2 objects ────────────────────────────────────────────────────
wild   <- readRDS(file.path(output_dir, "rds", "wild_macs2.rds"))
mutant <- readRDS(file.path(output_dir, "rds", "mutant_macs2.rds"))

# ── Function: process_expression ──────────────────────────────────────────
# SCTransform for integration; LogNormalize + VST for downstream DE
process_expression <- function(obj) {
    cat("  SCTransform...\n")
    obj <- SCTransform(obj)

    cat("  LogNormalize RNA + FindVariableFeatures...\n")
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                nfeatures = cfg$sct_nfeatures)
    return(obj)
}

# ── Function: process_accessibility ───────────────────────────────────────
# TF-IDF + SVD on both the original ATAC peaks and the MACS2-called peaks
process_accessibility <- function(obj) {
    # Original ATAC assay
    cat("  TF-IDF + SVD on ATAC assay...\n")
    DefaultAssay(obj) <- "ATAC"
    obj <- FindTopFeatures(obj, min.cutoff = "q0")
    obj <- RunTFIDF(obj)
    obj <- RunSVD(obj)

    # MACS2-called peaks assay
    cat("  TF-IDF + SVD on macs2 assay...\n")
    DefaultAssay(obj) <- "macs2"
    obj <- FindTopFeatures(obj, min.cutoff = "q0")
    obj <- RunTFIDF(obj)
    obj <- RunSVD(obj)

    return(obj)
}

# ── Process wild-type ─────────────────────────────────────────────────────
cat("Normalizing wild-type...\n")
wild <- process_expression(wild)
wild <- process_accessibility(wild)

# ── Process mutant ────────────────────────────────────────────────────────
cat("Normalizing mutant...\n")
mutant <- process_expression(mutant)
mutant <- process_accessibility(mutant)

# ── Save ──────────────────────────────────────────────────────────────────
saveRDS(wild,   file.path(output_dir, "rds", "wild_normalized.rds"))
saveRDS(mutant, file.path(output_dir, "rds", "mutant_normalized.rds"))
cat("=== Step 3 complete ===\n")