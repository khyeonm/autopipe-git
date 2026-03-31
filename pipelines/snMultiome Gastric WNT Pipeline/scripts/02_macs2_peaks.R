#!/usr/bin/env Rscript
# =============================================================================
#  02_macs2_peaks.R — MACS2 peak calling on scATAC-seq data
#
#  After QC filtering, this script calls peaks de novo using MACS2 (via
#  Signac::CallPeaks), removes blacklisted regions, and quantifies a new
#  peak-by-cell matrix stored as the "macs2" assay.
#
#  This replaces the default CellRanger peak calls with more sensitive,
#  sample-specific peak sets — critical for differential accessibility.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(yaml)
  library(future)
})

config <- yaml::read_yaml("/pipeline/config.yaml")
plan("multicore", workers = config$workers)
options(future.globals.maxSize = config$max_ram_gb * 1024^3)
set.seed(config$seed)

cat("=== Step 02: MACS2 Peak Calling ===\n")

# Load genome annotations (needed for CreateChromatinAssay)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

# =============================================================================
#  MACS2 peak calling + feature matrix quantification
# =============================================================================
# CallPeaks runs MACS2 on the fragment file attached to the ATAC assay.
# We then remove peaks on non-standard chromosomes and those overlapping
# the mm10 blacklist (ENCODE-defined problematic regions).
# FeatureMatrix re-counts fragments in the new peak set.
call_macs2_counts <- function(obj, species) {
  cat("  Calling MACS2 peaks...\n")
  peaks <- CallPeaks(obj, macs2.path = config$macs2_path, verbose = TRUE)

  # Select the appropriate blacklist for the species
  if (species == "mouse") {
    blacklist <- blacklist_mm10
  } else {
    blacklist <- blacklist_hg38_unified
  }

  # Keep only standard chromosomes (chr1-chr19, chrX, chrY for mouse)
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

  # Remove peaks overlapping blacklisted regions
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist, invert = TRUE)
  cat("  MACS2 peaks after filtering:", length(peaks), "\n")

  # Quantify fragment counts in the new peak set
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(obj),
    features  = peaks,
    cells     = colnames(obj)
  )

  return(macs2_counts)
}

# Add the MACS2 peak assay to the Seurat object
add_peak_to_object <- function(obj, species, frag_file) {
  obj[["macs2"]] <- CreateChromatinAssay(
    counts     = call_macs2_counts(obj, species),
    fragments  = frag_file,
    annotation = annotations
  )
  return(obj)
}

# =============================================================================
#  Process both samples
# =============================================================================
cat("\n--- Wild-type: MACS2 ---\n")
wild <- readRDS("/output/01_wild_qc.rds")
wild <- add_peak_to_object(wild, config$species, config$wild_type$fragments)
saveRDS(wild, "/output/02_wild_peaks.rds")
rm(wild); gc()

cat("\n--- Mutant: MACS2 ---\n")
mutant <- readRDS("/output/01_mutant_qc.rds")
mutant <- add_peak_to_object(mutant, config$species, config$mutant$fragments)
saveRDS(mutant, "/output/02_mutant_peaks.rds")
rm(mutant); gc()

cat("\n=== Step 02 complete ===\n")