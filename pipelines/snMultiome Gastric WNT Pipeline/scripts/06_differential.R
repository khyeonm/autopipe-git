#!/usr/bin/env Rscript
# =============================================================================
#  06_differential.R — Differential expression (Mutant vs Wild-type)
#
#  Uses MAST (Model-based Analysis of Single-cell Transcriptomics) to find
#  differentially expressed genes between Kras-mutant and wild-type cells.
#
#  Additionally performs KRAS pathway overlap analysis:
#    1. Retrieve MSigDB Hallmark KRAS_SIGNALING_UP gene set (human orthologs)
#    2. Convert mouse DEGs to human orthologs via gprofiler2::gconvert
#    3. Report intersection of upregulated DEGs with KRAS signaling genes
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(msigdbr)
  library(gprofiler2)
  library(yaml)
  library(future)
})

config <- yaml::read_yaml("/pipeline/config.yaml")
plan("multicore", workers = config$workers)
options(future.globals.maxSize = config$max_ram_gb * 1024^3)
set.seed(config$seed)

cat("=== Step 06: Differential Expression ===\n")

# Load integrated object
RZ.RZK <- readRDS("/output/04_integrated_clustered.rds")
DefaultAssay(RZ.RZK) <- "RNA"

# =============================================================================
#  1. Find DEGs: Mutant vs Wild-type (all cells, MAST test)
# =============================================================================
cat("  Running FindMarkers (MAST): Mutant vs Wild_type...\n")
DEG <- FindMarkers(
  RZ.RZK,
  group.by  = "kras_genotype",
  ident.1   = "Mutant",
  test.use  = config$de$test_use
)

# Save full DEG table
write.csv(DEG, "/output/06_DEG_MAST.csv")
cat("  Total DEGs tested:", nrow(DEG), "\n")

# Filter significant upregulated DEGs
sig_up <- DEG[DEG$p_val_adj < config$de$padj_threshold &
              DEG$avg_log2FC > config$de$log2fc_threshold, ]
cat("  Significant upregulated DEGs (padj <", config$de$padj_threshold,
    ", log2FC >", config$de$log2fc_threshold, "):", nrow(sig_up), "\n")

# =============================================================================
#  2. KRAS pathway overlap analysis
# =============================================================================
cat("  KRAS pathway overlap analysis...\n")

# Get Hallmark KRAS_SIGNALING_UP gene set (human genes)
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

kras_signaling <- hallmark[hallmark$gs_name == "HALLMARK_KRAS_SIGNALING_UP", ]
kras_genes <- kras_signaling$gene_symbol
cat("  Hallmark KRAS_SIGNALING_UP genes:", length(kras_genes), "\n")

# Convert mouse DEG gene names to human orthologs
mouse_deg_genes <- rownames(sig_up)
human_orthologs <- gconvert(
  query    = mouse_deg_genes,
  organism = "hsapiens",
  target   = "ENSG"
)

cat("  Mouse DEGs with human orthologs:", length(human_orthologs$name), "\n")

# Find intersection with KRAS signaling
overlap_genes <- intersect(human_orthologs$name, kras_genes)
cat("  Overlap with KRAS_SIGNALING_UP:", length(overlap_genes), "genes\n")

# Save overlap results
overlap_df <- data.frame(
  kras_overlap_gene   = overlap_genes,
  stringsAsFactors    = FALSE
)
write.csv(overlap_df, "/output/06_DEG_KRAS_overlap.csv", row.names = FALSE)

cat("\n=== Step 06 complete ===\n")