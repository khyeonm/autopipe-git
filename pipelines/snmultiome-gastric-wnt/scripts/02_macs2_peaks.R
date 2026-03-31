#!/usr/bin/env Rscript
# ===========================================================================
#  02_macs2_peaks.R — MACS2 peak calling & ChromatinAssay creation
#  Replicates: call_macs2_counts() + add_peak_to_object() from original code
#  Inputs:  QC-filtered Seurat objects
#  Outputs: Objects with added "macs2" ChromatinAssay
# ===========================================================================

suppressPackageStartupMessages({
    library(yaml)
    library(Seurat)
    library(Signac)
    library(EnsDb.Mmusculus.v79)
    library(BSgenome.Mmusculus.UCSC.mm10)
})

# ── Load configuration ────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
config_path <- args[which(args == "--config") + 1]
cfg <- read_yaml(config_path)

set.seed(cfg$seed)
output_dir <- cfg$output_dir

cat("=== Step 2: MACS2 peak calling ===\n")

# ── Set genome annotation ─────────────────────────────────────────────────
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

# ── Load QC objects ────────────────────────────────────────────────────────
wild   <- readRDS(file.path(output_dir, "rds", "wild_qc.rds"))
mutant <- readRDS(file.path(output_dir, "rds", "mutant_qc.rds"))

# ── Function: call_macs2_counts ───────────────────────────────────────────
# Calls peaks with MACS2, removes blacklist regions (mm10), then
# quantifies fragment counts in the called peaks.
call_macs2_counts <- function(obj) {
    cat("  Calling peaks with MACS2...\n")
    peaks <- CallPeaks(obj, macs2.path = cfg$macs2_path, verbose = TRUE)

    # Remove peaks on nonstandard chromosomes
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

    # Remove peaks overlapping mm10 blacklist regions
    peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

    cat(sprintf("  MACS2 peaks after filtering: %d\n", length(peaks)))

    # Quantify counts in each peak per cell
    macs2_counts <- FeatureMatrix(
        fragments = Fragments(obj),
        features  = peaks,
        cells     = colnames(obj)
    )
    return(macs2_counts)
}

# ── Function: add_peak_to_object ──────────────────────────────────────────
# Creates a new ChromatinAssay named "macs2" from MACS2 peak counts
add_peak_to_object <- function(obj, frag_file) {
    obj[["macs2"]] <- CreateChromatinAssay(
        counts     = call_macs2_counts(obj),
        fragments  = frag_file,
        annotation = annotations
    )
    return(obj)
}

# ── Process wild-type ─────────────────────────────────────────────────────
cat("Processing wild-type MACS2 peaks...\n")
wild <- add_peak_to_object(wild, cfg$wild_type$fragments)
cat(sprintf("  Wild-type macs2 peaks: %d\n", nrow(wild[["macs2"]])))

# ── Process mutant ────────────────────────────────────────────────────────
cat("Processing mutant MACS2 peaks...\n")
mutant <- add_peak_to_object(mutant, cfg$mutant$fragments)
cat(sprintf("  Mutant macs2 peaks: %d\n", nrow(mutant[["macs2"]])))

# ── Save ──────────────────────────────────────────────────────────────────
saveRDS(wild,   file.path(output_dir, "rds", "wild_macs2.rds"))
saveRDS(mutant, file.path(output_dir, "rds", "mutant_macs2.rds"))
cat("=== Step 2 complete ===\n")