# ==============================================================================
# 02_macs2_peaks.R — MACS2 peak calling & new ChromatinAssay
# ==============================================================================
# Calls peaks on the QC'd objects using MACS2, removes blacklist regions
# and non-standard chromosomes, then quantifies a new "macs2" assay.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(future)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
cfg  <- read_yaml(args[1])

plan("multicore", workers = cfg$threads)
options(future.globals.maxSize = cfg$max_memory_gb * 1024^3)
set.seed(1234)

outdir <- cfg$output_dir

# ---- Genome annotation (needed for ChromatinAssay) -------------------------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

# ---- Helper: MACS2 peak calling --------------------------------------------
call_macs2_counts <- function(obj, species, macs2_path) {
  cat("   Calling peaks with MACS2...\n")
  peaks <- CallPeaks(obj, macs2.path = macs2_path, verbose = TRUE)

  # Remove non-standard chromosomes
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

  # Remove blacklist regions
  if (species == "mouse") {
    peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
  } else {
    peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  }

  # Quantify counts in each peak
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(obj),
    features  = peaks,
    cells     = colnames(obj)
  )
  return(macs2_counts)
}

add_peak_to_object <- function(obj, species, frag_file, macs2_path, annot) {
  macs2_counts <- call_macs2_counts(obj, species, macs2_path)
  obj[["macs2"]] <- CreateChromatinAssay(
    counts     = macs2_counts,
    fragments  = frag_file,
    annotation = annot
  )
  return(obj)
}

# ---- Load QC'd objects ------------------------------------------------------
cat(">> Loading QC'd objects\n")
wild   <- readRDS(file.path(outdir, "01_qc", "wild_type_qc.rds"))
mutant <- readRDS(file.path(outdir, "01_qc", "mutant_qc.rds"))

# ---- Call MACS2 peaks -------------------------------------------------------
cat(">> Wild-type MACS2 peak calling\n")
wild <- add_peak_to_object(wild, cfg$species, cfg$wild_type$fragment_file,
                           cfg$macs2_path, annotations)

cat(">> Mutant MACS2 peak calling\n")
mutant <- add_peak_to_object(mutant, cfg$species, cfg$mutant$fragment_file,
                             cfg$macs2_path, annotations)

# ---- Save -------------------------------------------------------------------
dir.create(file.path(outdir, "02_peaks"), recursive = TRUE, showWarnings = FALSE)
saveRDS(wild,   file.path(outdir, "02_peaks", "wild_type_peaks.rds"))
saveRDS(mutant, file.path(outdir, "02_peaks", "mutant_peaks.rds"))
cat(">> Step 02 complete.\n")