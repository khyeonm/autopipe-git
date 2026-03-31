# ==============================================================================
# 03_normalize.R — SCTransform (RNA) + TF-IDF/SVD (ATAC & macs2)
# ==============================================================================
# RNA:   SCTransform → stored in SCT assay (for integration)
#        Also log-normalises the RNA assay (for DE analysis)
# ATAC:  FindTopFeatures(q0) → RunTFIDF → RunSVD  (both ATAC and macs2 assays)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(future)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
cfg  <- read_yaml(args[1])

plan("multicore", workers = cfg$threads)
options(future.globals.maxSize = cfg$max_memory_gb * 1024^3)
set.seed(1234)

outdir <- cfg$output_dir

# ---- Helper: RNA normalisation ----------------------------------------------
process_expression <- function(obj, nfeatures) {
  cat("   SCTransform...\n")
  obj <- SCTransform(obj)

  # Also log-normalise the RNA assay for downstream DE
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  return(obj)
}

# ---- Helper: ATAC normalisation ---------------------------------------------
process_accessibility <- function(obj) {
  # Process the Cell Ranger peaks assay
  DefaultAssay(obj) <- "ATAC"
  obj <- FindTopFeatures(obj, min.cutoff = "q0")
  obj <- RunTFIDF(obj)
  obj <- RunSVD(obj)

  # Process the MACS2 peaks assay
  DefaultAssay(obj) <- "macs2"
  obj <- FindTopFeatures(obj, min.cutoff = "q0")
  obj <- RunTFIDF(obj)
  obj <- RunSVD(obj)
  return(obj)
}

# ---- Load peak-called objects -----------------------------------------------
cat(">> Loading peak-called objects\n")
wild   <- readRDS(file.path(outdir, "02_peaks", "wild_type_peaks.rds"))
mutant <- readRDS(file.path(outdir, "02_peaks", "mutant_peaks.rds"))

# ---- Normalise RNA ----------------------------------------------------------
cat(">> Normalising wild-type RNA\n")
wild <- process_expression(wild, cfg$n_variable_features)
cat(">> Normalising mutant RNA\n")
mutant <- process_expression(mutant, cfg$n_variable_features)

# ---- Normalise ATAC ---------------------------------------------------------
cat(">> Normalising wild-type ATAC\n")
wild <- process_accessibility(wild)
cat(">> Normalising mutant ATAC\n")
mutant <- process_accessibility(mutant)

# ---- Save -------------------------------------------------------------------
dir.create(file.path(outdir, "03_normalised"), recursive = TRUE, showWarnings = FALSE)
saveRDS(wild,   file.path(outdir, "03_normalised", "wild_type_normalised.rds"))
saveRDS(mutant, file.path(outdir, "03_normalised", "mutant_normalised.rds"))
cat(">> Step 03 complete.\n")