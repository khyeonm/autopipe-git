# ==============================================================================
# 01_load_and_qc.R — Load 10X snMultiome data, compute QC, filter
# ==============================================================================
# Reads raw Cell Ranger ARC output (filtered_feature_bc_matrix.h5 +
# atac_fragments.tsv.gz + per_barcode_metrics.csv) for wild-type (RZ)
# and mutant (RZK) samples. Computes QC metrics (nucleosome signal,
# TSS enrichment, percent mito) and filters cells.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(future)
  library(yaml)
  library(ggplot2)
})

# ---- Read config -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
cfg  <- read_yaml(args[1])

# ---- Parallelism ----------------------------------------------------------
plan("multicore", workers = cfg$threads)
options(future.globals.maxSize = cfg$max_memory_gb * 1024^3)
future.seed <- TRUE
set.seed(1234)

outdir <- cfg$output_dir

# ---- Genome annotation -----------------------------------------------------
cat(">> Setting genome annotation for", cfg$species, "\n")
if (cfg$species == "mouse") {
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "mm10"
} else {
  stop("Only 'mouse' species is currently supported in this pipeline.")
}

# ---- Helper: generate Seurat multiome object --------------------------------
generate_object <- function(h5_file, frag_file, meta_file, species, annot) {
  cat("   Loading h5:", h5_file, "\n")
  metadata   <- read.csv(file = meta_file, header = TRUE, row.names = 1)
  input.10x  <- Read10X_h5(h5_file)

  # --- RNA assay ---
  rna.counts <- input.10x$`Gene Expression`
  obj <- CreateSeuratObject(counts = rna.counts, assay = "RNA", meta.data = metadata)

  # Percent mitochondrial
  if (species == "mouse") {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  } else {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  }

  # --- ATAC assay (standard chromosomes only) ---
  atac.counts  <- input.10x$Peaks
  grange.counts <- StringToGRanges(rownames(atac.counts), sep = c(":", "-"))
  grange.use   <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac.counts  <- atac.counts[as.vector(grange.use), ]

  obj[["ATAC"]] <- CreateChromatinAssay(
    counts     = atac.counts,
    sep        = c(":", "-"),
    genome     = genome(annot),
    fragments  = frag_file,
    min.cells  = 10,
    annotation = annot
  )
  return(obj)
}

# ---- Helper: QC metrics + subsetting ----------------------------------------
qc_and_filter <- function(obj, cfg) {
  DefaultAssay(obj) <- "ATAC"
  obj <- NucleosomeSignal(obj)
  obj <- TSSEnrichment(obj)

  cat("   Cells before QC:", ncol(obj), "\n")
  obj <- subset(
    x = obj,
    subset = nCount_ATAC < cfg$qc$max_nCount_ATAC &
             nCount_RNA  < cfg$qc$max_nCount_RNA  &
             nCount_ATAC > cfg$qc$min_nCount_ATAC &
             nCount_RNA  > cfg$qc$min_nCount_RNA  &
             nucleosome_signal < cfg$qc$max_nucleosome_signal &
             TSS.enrichment    > cfg$qc$min_TSS_enrichment &
             percent.mt        < cfg$qc$max_percent_mt
  )
  cat("   Cells after QC: ", ncol(obj), "\n")
  return(obj)
}

# ---- Process wild-type (RZ) ------------------------------------------------
cat(">> Processing wild-type sample\n")
wild <- generate_object(
  h5_file   = cfg$wild_type$h5_file,
  frag_file = cfg$wild_type$fragment_file,
  meta_file = cfg$wild_type$meta_file,
  species   = cfg$species,
  annot     = annotations
)
wild@meta.data$kras_genotype <- "Wild_type"

# QC violin plot before filtering
DefaultAssay(wild) <- "ATAC"
wild <- NucleosomeSignal(wild)
wild <- TSSEnrichment(wild)
dir.create(file.path(outdir, "01_qc"), recursive = TRUE, showWarnings = FALSE)
pdf(file.path(outdir, "01_qc", "qc_violin_wild_type.pdf"), width = 14, height = 5)
print(VlnPlot(wild,
              features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
              ncol = 4, pt.size = 0))
dev.off()

wild <- qc_and_filter(wild, cfg)

# ---- Process mutant (RZK) --------------------------------------------------
cat(">> Processing mutant sample\n")
mutant <- generate_object(
  h5_file   = cfg$mutant$h5_file,
  frag_file = cfg$mutant$fragment_file,
  meta_file = cfg$mutant$meta_file,
  species   = cfg$species,
  annot     = annotations
)
mutant@meta.data$kras_genotype <- "Mutant"

DefaultAssay(mutant) <- "ATAC"
mutant <- NucleosomeSignal(mutant)
mutant <- TSSEnrichment(mutant)
pdf(file.path(outdir, "01_qc", "qc_violin_mutant.pdf"), width = 14, height = 5)
print(VlnPlot(mutant,
              features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
              ncol = 4, pt.size = 0))
dev.off()

mutant <- qc_and_filter(mutant, cfg)

# ---- Save ------------------------------------------------------------------
saveRDS(wild,   file.path(outdir, "01_qc", "wild_type_qc.rds"))
saveRDS(mutant, file.path(outdir, "01_qc", "mutant_qc.rds"))
cat(">> Step 01 complete.\n")