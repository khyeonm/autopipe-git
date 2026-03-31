#!/usr/bin/env Rscript
# =============================================================================
#  01_load_and_qc.R — Load 10X snMultiome data & apply QC filters
#
#  This script loads paired snRNA-seq + scATAC-seq data from 10X Genomics
#  CellRanger-ARC output (filtered_feature_bc_matrix.h5 + atac_fragments.tsv.gz)
#  for two conditions: Wild-type (RZ) and Mutant (RZK, Kras-driven).
#
#  Workflow:
#    1. Set genome annotations (EnsDb.Mmusculus.v79 for mm10)
#    2. Read 10X multiome h5 → separate RNA and ATAC count matrices
#    3. Create Seurat objects with ChromatinAssay for ATAC
#    4. Compute QC metrics: %mt, nucleosome signal, TSS enrichment
#    5. Subset cells passing all QC thresholds from config.yaml
#    6. Save QC-filtered objects as RDS
# =============================================================================

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(dplyr)
  library(ggplot2)
  library(yaml)
  library(future)
})

# ---- Read config ----
config <- yaml::read_yaml("/pipeline/config.yaml")

# ---- Set up parallelism ----
plan("multicore", workers = config$workers)
options(future.globals.maxSize = config$max_ram_gb * 1024^3)
future.seed <- TRUE
set.seed(config$seed)

cat("=== Step 01: Load & QC ===\n")
cat("Species:", config$species, "\n")

# =============================================================================
#  1. Set genome annotations
# =============================================================================
# GetGRangesFromEnsDb extracts gene-level GRanges from the Ensembl database.
# We enforce UCSC-style chromosome names (chr1, chr2, ...) and tag them as mm10.
set_annotation <- function(species) {
  if (species == "mouse") {
    annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotation) <- "UCSC"
    genome(annotation) <- "mm10"
  } else {
    stop("Only 'mouse' species is supported in this pipeline.")
  }
  return(annotation)
}

annotations <- set_annotation(config$species)
cat("Annotations loaded:", length(annotations), "features\n")

# =============================================================================
#  2. Build a Seurat multiome object from 10X CellRanger-ARC output
# =============================================================================
# Read10X_h5 returns a named list: $`Gene Expression` and $Peaks.
# We create the RNA assay first, then add the ATAC ChromatinAssay.
# Only standard chromosomes are retained for ATAC peaks.
generate_object <- function(h5_file, frag_file, meta_file, species, annotations) {
  cat("  Loading h5:", h5_file, "\n")

  # Read per-barcode metadata from CellRanger-ARC
  metadata <- read.csv(file = meta_file, header = TRUE, row.names = 1)

  # Read the multiome h5 file (contains both Gene Expression and Peaks)
  input.10x <- Read10X_h5(h5_file)
  rna.counts <- input.10x$`Gene Expression`

  # Create Seurat object with RNA assay
  obj <- CreateSeuratObject(counts = rna.counts, assay = "RNA", meta.data = metadata)

  # Compute mitochondrial percentage
  # Mouse mt genes start with lowercase "mt-"
  if (species == "mouse") {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  } else {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  }

  # Prepare ATAC counts: keep only standard chromosomes
  atac.counts <- input.10x$Peaks
  grange.counts <- StringToGRanges(rownames(atac.counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac.counts <- atac.counts[as.vector(grange.use), ]

  # Create ChromatinAssay with fragment file and genome annotations
  obj[["ATAC"]] <- CreateChromatinAssay(
    counts     = atac.counts,
    sep        = c(":", "-"),
    genome     = genome(annotations),
    fragments  = frag_file,
    min.cells  = config$processing$chromatin_assay_min_cells,
    annotation = annotations
  )

  cat("  Object created:", ncol(obj), "cells\n")
  return(obj)
}

# =============================================================================
#  3. QC filtering
# =============================================================================
# Compute ATAC-specific QC metrics (nucleosome signal, TSS enrichment) then
# apply multi-modal thresholds from config.yaml.
subsetting_objects_QC <- function(obj) {
  DefaultAssay(obj) <- "ATAC"

  # NucleosomeSignal: ratio of mono-nucleosomal to nucleosome-free fragments
  obj <- NucleosomeSignal(obj)

  # TSSEnrichment: fold-enrichment of fragments at transcription start sites
  obj <- TSSEnrichment(obj)

  cat("  Pre-QC cells:", ncol(obj), "\n")

  # Apply QC thresholds (from the original RZ/RZK analysis code)
  obj <- subset(
    x = obj,
    subset = nCount_ATAC < config$qc$nCount_ATAC_max &
      nCount_RNA  < config$qc$nCount_RNA_max &
      nCount_ATAC > config$qc$nCount_ATAC_min &
      nCount_RNA  > config$qc$nCount_RNA_min &
      nucleosome_signal < config$qc$nucleosome_signal_max &
      TSS.enrichment    > config$qc$TSS_enrichment_min &
      percent.mt        < config$qc$percent_mt_max
  )

  cat("  Post-QC cells:", ncol(obj), "\n")
  return(obj)
}

# =============================================================================
#  4. Execute for both samples
# =============================================================================
# --- Wild-type (RZ) ---
cat("\n--- Wild-type (RZ) ---\n")
wild <- generate_object(
  config$wild_type$h5,
  config$wild_type$fragments,
  config$wild_type$metadata,
  config$species,
  annotations
)
wild@meta.data$kras_genotype <- "Wild_type"
wild <- subsetting_objects_QC(wild)

# --- Mutant (RZK) ---
cat("\n--- Mutant (RZK) ---\n")
mutant <- generate_object(
  config$mutant$h5,
  config$mutant$fragments,
  config$mutant$metadata,
  config$species,
  annotations
)
mutant@meta.data$kras_genotype <- "Mutant"
mutant <- subsetting_objects_QC(mutant)

# =============================================================================
#  5. Save QC-filtered objects
# =============================================================================
dir.create("/output", showWarnings = FALSE, recursive = TRUE)
saveRDS(wild,   "/output/01_wild_qc.rds")
saveRDS(mutant, "/output/01_mutant_qc.rds")

cat("\n=== Step 01 complete ===\n")
cat("Wild-type cells:", ncol(wild), "\n")
cat("Mutant cells:   ", ncol(mutant), "\n")