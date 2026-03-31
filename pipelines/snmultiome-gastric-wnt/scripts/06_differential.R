#!/usr/bin/env Rscript
# ===========================================================================
#  06_differential.R — Differential gene expression (Mutant vs Wild-type)
#  Replicates: FindMarkers(group.by="kras_genotype", test.use="MAST")
#              + MSigDB Hallmark KRAS_SIGNALING_UP overlap analysis
#
#  Outputs:
#    - DEG_scMAST.csv: full DEG table
#    - KRAS_signaling_overlap.csv: overlap between significant DEGs and
#      Hallmark KRAS signaling genes
# ===========================================================================

suppressPackageStartupMessages({
    library(yaml)
    library(future)
    library(Seurat)
    library(Signac)
    library(MAST)
    library(msigdbr)
    library(gprofiler2)
})

# ── Load configuration ────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
config_path <- args[which(args == "--config") + 1]
cfg <- read_yaml(config_path)

set.seed(cfg$seed)
plan("multicore", workers = cfg$workers)
options(future.globals.maxSize = cfg$future_ram_gb * 1024^3)
output_dir <- cfg$output_dir

dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

cat("=== Step 6: Differential gene expression ===\n")

# ── Load integrated object ────────────────────────────────────────────────
obj <- readRDS(file.path(output_dir, "rds", "integrated_clustered.rds"))
DefaultAssay(obj) <- "RNA"

# ── DEG: Mutant vs Wild-type (MAST test) ──────────────────────────────────
# Exactly as in the original code: FindMarkers with group.by="kras_genotype"
cat("Running FindMarkers (MAST): Mutant vs Wild_type...\n")
DEG <- FindMarkers(obj,
    group.by = "kras_genotype",
    ident.1  = "Mutant",
    test.use = cfg$de_test
)

# Save full DEG table
deg_path <- file.path(output_dir, "tables", "DEG_scMAST.csv")
write.csv(DEG, deg_path)
cat(sprintf("DEG table saved: %d genes tested\n", nrow(DEG)))

# ── Filter significant upregulated DEGs ───────────────────────────────────
# Original filter: p_val_adj < 0.01 & avg_log2FC > 1
sig_up <- DEG[DEG$p_val_adj < cfg$de_pval_cutoff &
              DEG$avg_log2FC > cfg$de_log2fc_cutoff, ]
cat(sprintf("Significant upregulated DEGs: %d\n", nrow(sig_up)))

# ── MSigDB Hallmark KRAS_SIGNALING_UP overlap ────────────────────────────
# Original code: get human Hallmark, convert mouse DEG gene names to human
# via gconvert, then intersect
cat("Fetching Hallmark KRAS_SIGNALING_UP gene set...\n")
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark <- hallmark[, c("gs_name", "gene_symbol")]

kras_signaling <- hallmark$gene_symbol[
    hallmark$gs_name == "HALLMARK_KRAS_SIGNALING_UP"
]
cat(sprintf("KRAS_SIGNALING_UP genes: %d\n", length(kras_signaling)))

# Convert mouse gene symbols to human orthologs
cat("Converting mouse gene symbols to human orthologs (gProfiler)...\n")
mouse_genes <- rownames(sig_up)
converted <- gconvert(query = mouse_genes, organism = "hsapiens", target = "ENSG")

cat(sprintf("Converted genes: %d\n", length(converted$name)))

# Find intersection
overlap_genes <- intersect(converted$name, kras_signaling)
cat(sprintf("Overlap with KRAS_SIGNALING_UP: %d genes\n", length(overlap_genes)))

# Save overlap results
overlap_df <- data.frame(
    overlap_gene      = overlap_genes,
    in_KRAS_signaling = TRUE,
    stringsAsFactors  = FALSE
)
overlap_path <- file.path(output_dir, "tables", "KRAS_signaling_overlap.csv")
write.csv(overlap_df, overlap_path, row.names = FALSE)

# ── Summary statistics ────────────────────────────────────────────────────
cat(sprintf("\n--- DEG Summary ---\n"))
cat(sprintf("Total genes tested:              %d\n", nrow(DEG)))
cat(sprintf("Significant upregulated (padj<%.2f, log2FC>%.1f): %d\n",
            cfg$de_pval_cutoff, cfg$de_log2fc_cutoff, nrow(sig_up)))
cat(sprintf("KRAS_SIGNALING_UP overlap:       %d / %d\n",
            length(overlap_genes), length(kras_signaling)))

cat("=== Step 6 complete ===\n")