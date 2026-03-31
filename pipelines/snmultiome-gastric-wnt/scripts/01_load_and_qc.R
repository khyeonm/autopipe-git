#!/usr/bin/env Rscript
# ===========================================================================
#  01_load_and_qc.R — Load 10X snMultiome data & apply QC filters
#  Inputs:  10X .h5 matrices, ATAC fragment files, barcode metrics
#  Outputs: QC-filtered Seurat objects (wild_qc.rds, mutant_qc.rds)
# ===========================================================================

suppressPackageStartupMessages({
    library(yaml)
    library(Seurat)
    library(Signac)
    library(EnsDb.Mmusculus.v79)
    library(ggplot2)
})

# ── Load configuration ────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
config_path <- args[which(args == "--config") + 1]
cfg <- read_yaml(config_path)

set.seed(cfg$seed)
output_dir <- cfg$output_dir

# Create output directories
dir.create(file.path(output_dir, "rds"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"),  recursive = TRUE, showWarnings = FALSE)

cat("=== Step 1: Load data & QC ===\n")

# ── Set genome annotation (mouse mm10 via EnsDb.Mmusculus.v79) ────────────
# Exactly as in the original code: set_annotation("mouse")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"
cat("Genome annotation set: mm10 (EnsDb.Mmusculus.v79)\n")

# ── Function: generate_object ─────────────────────────────────────────────
# Replicates the original generate_object() function exactly:
#   - Reads 10X .h5 (Gene Expression + Peaks)
#   - Creates Seurat object with RNA assay
#   - Computes percent.mt (pattern: ^mt- for mouse)
#   - Filters ATAC peaks to standard chromosomes
#   - Creates ChromatinAssay with min.cells=10 and the EnsDb annotation
generate_object <- function(h5_file, frag_file, meta_file) {
    cat(sprintf("  Loading: %s\n", h5_file))

    # Read barcode metrics as metadata
    metadata <- read.csv(file = meta_file, header = TRUE, row.names = 1)

    # Read multiome .h5 (contains Gene Expression + Peaks)
    input.10x <- Read10X_h5(h5_file)
    rna.counts <- input.10x$`Gene Expression`

    # Create Seurat object with RNA counts
    obj <- CreateSeuratObject(counts = rna.counts, assay = "RNA", meta.data = metadata)

    # Calculate mitochondrial gene percentage (mouse: ^mt-)
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

    # Prepare ATAC counts — keep only standard chromosomes
    atac.counts <- input.10x$Peaks
    grange.counts <- StringToGRanges(rownames(atac.counts), sep = c(":", "-"))
    grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac.counts <- atac.counts[as.vector(grange.use), ]

    # Create ChromatinAssay with genome annotation
    obj[["ATAC"]] <- CreateChromatinAssay(
        counts     = atac.counts,
        sep        = c(":", "-"),
        genome     = genome(annotations),
        fragments  = frag_file,
        min.cells  = 10,
        annotation = annotations
    )

    return(obj)
}

# ── Function: qc_function ─────────────────────────────────────────────────
# Computes NucleosomeSignal and TSSEnrichment on the ATAC assay
qc_function <- function(obj) {
    DefaultAssay(obj) <- "ATAC"
    obj <- NucleosomeSignal(obj)
    obj <- TSSEnrichment(obj)
    return(obj)
}

# ── Function: subsetting_objects_QC ───────────────────────────────────────
# Applies the exact QC thresholds from the original code:
#   nCount_ATAC:       1,000 – 100,000
#   nCount_RNA:        1,000 – 25,000
#   nucleosome_signal: < 2
#   TSS.enrichment:    > 1
#   percent.mt:        < 15
subsetting_objects_QC <- function(obj) {
    obj <- qc_function(obj)
    obj <- subset(
        x = obj,
        subset = nCount_ATAC < cfg$qc$nCount_ATAC_max &
                 nCount_RNA  < cfg$qc$nCount_RNA_max  &
                 nCount_ATAC > cfg$qc$nCount_ATAC_min &
                 nCount_RNA  > cfg$qc$nCount_RNA_min  &
                 nucleosome_signal < cfg$qc$nucleosome_signal_max &
                 TSS.enrichment    > cfg$qc$TSS_enrichment_min    &
                 percent.mt        < cfg$qc$percent_mt_max
    )
    return(obj)
}

# ── Process wild-type sample ──────────────────────────────────────────────
cat("Processing wild-type sample...\n")
wild <- generate_object(
    h5_file   = cfg$wild_type$h5,
    frag_file = cfg$wild_type$fragments,
    meta_file = cfg$wild_type$metadata
)
wild@meta.data$kras_genotype <- cfg$wild_type$label
cat(sprintf("  Wild-type cells before QC: %d\n", ncol(wild)))

wild <- subsetting_objects_QC(wild)
cat(sprintf("  Wild-type cells after QC:  %d\n", ncol(wild)))

# ── Process mutant sample ─────────────────────────────────────────────────
cat("Processing mutant sample...\n")
mutant <- generate_object(
    h5_file   = cfg$mutant$h5,
    frag_file = cfg$mutant$fragments,
    meta_file = cfg$mutant$metadata
)
mutant@meta.data$kras_genotype <- cfg$mutant$label
cat(sprintf("  Mutant cells before QC: %d\n", ncol(mutant)))

mutant <- subsetting_objects_QC(mutant)
cat(sprintf("  Mutant cells after QC:  %d\n", ncol(mutant)))

# ── Generate QC violin plots ─────────────────────────────────────────────
pdf(file.path(output_dir, "plots", "qc_violin_plots.pdf"), width = 16, height = 8)
p_wt <- VlnPlot(wild,   features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 4, pt.size = 0) + ggtitle("Wild-type QC")
p_mt <- VlnPlot(mutant, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 4, pt.size = 0) + ggtitle("Mutant QC")
print(p_wt)
print(p_mt)
dev.off()
cat("QC violin plots saved.\n")

# ── Save QC-filtered objects ──────────────────────────────────────────────
saveRDS(wild,   file.path(output_dir, "rds", "wild_qc.rds"))
saveRDS(mutant, file.path(output_dir, "rds", "mutant_qc.rds"))
cat("=== Step 1 complete ===\n")