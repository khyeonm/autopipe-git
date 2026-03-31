#!/usr/bin/env Rscript
# =============================================================================
#  07_motif_enrichment.R — Differential accessibility + TF motif enrichment
#
#  Performs differential accessibility (DA) analysis between Wnt7+ and Lgr5+
#  populations in both mutant and wild-type conditions, then identifies
#  enriched TF motifs in the DA peak sets using JASPAR2020.
#
#  DA analysis uses logistic regression (LR) with ATAC fragment counts
#  as a latent variable to control for sequencing depth differences.
#
#  Motif analysis workflow:
#    1. AddMotifs: scan all peaks for JASPAR2020 CORE vertebrate motifs
#    2. FindMotifs: hypergeometric test for motif enrichment in DA peaks
#    3. Visualize top enriched motifs as bar plots
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(motifmatchr)
  library(JASPAR2020)
  library(TFBSTools)
  library(chromVAR)
  library(yaml)
  library(future)
})

config <- yaml::read_yaml("/pipeline/config.yaml")
# Use single worker for DA analysis (as in the original code: plan("multicore", workers = 1))
plan("multicore", workers = 1)
options(future.globals.maxSize = config$max_ram_gb * 1024^3)
set.seed(config$seed)

cat("=== Step 07: Differential Accessibility & Motif Enrichment ===\n")

fig_dir <- "/output/05_figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Load integrated object
RZ.RZK <- readRDS("/output/04_integrated_clustered.rds")
linkage_genome <- BSgenome.Mmusculus.UCSC.mm10

# =============================================================================
#  1. Differential accessibility: Wnt7+ vs Lgr5+ (Mutant)
# =============================================================================
cat("  DA peaks: Wnt7+ vs Lgr5+ in Mutant...\n")
DefaultAssay(RZ.RZK) <- "ATAC"

da_peaks_mutant <- FindMarkers(
  RZ.RZK,
  ident.1    = "3_Mutant",
  ident.2    = "1_Mutant",
  group.by   = "cl_genotype",
  min.pct    = config$da$min_pct,
  logfc.threshold = config$da$logfc_threshold,
  test.use   = config$da$test_use,
  latent.vars = config$da$latent_vars
)

# Annotate DA peaks with closest gene
closest_genes <- ClosestFeature(RZ.RZK, regions = rownames(da_peaks_mutant))
rownames(closest_genes) <- closest_genes$query_region
da_peaks_mutant$gene         <- closest_genes$gene_name
da_peaks_mutant$gene_biotype <- closest_genes$gene_biotype
da_peaks_mutant$distance     <- closest_genes$distance

# Filter significant DA peaks
sig_da_mut <- da_peaks_mutant[da_peaks_mutant$p_val_adj < 0.05, ]
cat("  Significant DA peaks (Mutant):", nrow(sig_da_mut), "\n")

write.csv(da_peaks_mutant, "/output/07_DA_peaks_wnt_vs_lgr5_mutant.csv")

# =============================================================================
#  2. Differential accessibility: Wnt7+ vs Lgr5+ (Wild-type)
# =============================================================================
cat("  DA peaks: Wnt7+ vs Lgr5+ in Wild-type...\n")

da_peaks_wildtype <- FindMarkers(
  RZ.RZK,
  ident.1    = "3_Wild_type",
  ident.2    = "1_Wild_type",
  group.by   = "cl_genotype",
  min.pct    = config$da$min_pct,
  logfc.threshold = config$da$logfc_threshold,
  test.use   = config$da$test_use,
  latent.vars = config$da$latent_vars
)

closest_genes <- ClosestFeature(RZ.RZK, regions = rownames(da_peaks_wildtype))
rownames(closest_genes) <- closest_genes$query_region
da_peaks_wildtype$gene         <- closest_genes$gene_name
da_peaks_wildtype$gene_biotype <- closest_genes$gene_biotype
da_peaks_wildtype$distance     <- closest_genes$distance

sig_da_wt <- da_peaks_wildtype[da_peaks_wildtype$p_val_adj < 0.05, ]
cat("  Significant DA peaks (Wild-type):", nrow(sig_da_wt), "\n")

write.csv(da_peaks_wildtype, "/output/07_DA_peaks_wnt_vs_lgr5_wildtype.csv")

# =============================================================================
#  3. Motif enrichment analysis
# =============================================================================
cat("  Adding JASPAR2020 motifs...\n")

# Get JASPAR2020 position frequency matrices (CORE vertebrate collection)
pfm <- getMatrixSet(
  x    = JASPAR2020,
  opts = list(
    collection  = config$motif$jaspar_collection,
    tax_group   = config$motif$jaspar_tax_group,
    all_versions = FALSE
  )
)

# Scan all ATAC peaks for motif occurrences
RZ.RZK <- AddMotifs(object = RZ.RZK, genome = linkage_genome, pfm = pfm)

# ---- Motif enrichment: Mutant (Wnt7+ vs Lgr5+) ----
cat("  Motif enrichment: Mutant...\n")
sig_peaks_mut <- rownames(
  da_peaks_mutant[da_peaks_mutant$p_val_adj < config$da$padj_threshold, ]
)
motif_mut <- FindMotifs(object = RZ.RZK, features = sig_peaks_mut)
motif_mut$logP <- -log10(motif_mut$p.adjust)

sig_motif_mut <- motif_mut[motif_mut$p.adjust < config$motif$enriched_padj, ]
cat("  Significant motifs (Mutant):", nrow(sig_motif_mut), "\n")

write.csv(motif_mut, "/output/07_motif_enrichment_mutant.csv", row.names = FALSE)

# Plot top 15 enriched motifs (Mutant)
data_mut <- motif_mut[1:min(15, nrow(motif_mut)), ]
y_order_mut <- rev(data_mut$motif.name)

p_mot_mut <- ggplot(data_mut, aes(x = logP, y = motif.name, fill = percent.observed)) +
  geom_bar(stat = "identity", width = 0.2) +
  scale_fill_gradientn(colours = c("royalblue", "rosybrown2", "red"), limits = c(0, 60)) +
  scale_y_discrete(limits = y_order_mut) +
  xlim(0, 30) +
  ggtitle("Wnt7b+ vs Lgr5+\n(RZK / Mutant)") +
  xlab("-log10(p.adjust)") + ylab("TF motif") +
  theme_linedraw() +
  theme(
    plot.title   = element_text(hjust = 0.5, color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(color = "black", size = 12),
    axis.title.y = element_text(color = "black", size = 12)
  )

ggsave(p_mot_mut,
       filename = file.path(fig_dir, "fig3I_motif_enrichment_mutant.pdf"),
       width = 6, height = 6, dpi = 300)

# ---- Motif enrichment: Wild-type (Wnt7+ vs Lgr5+) ----
cat("  Motif enrichment: Wild-type...\n")
sig_peaks_wt <- rownames(
  da_peaks_wildtype[da_peaks_wildtype$p_val_adj < config$da$padj_threshold, ]
)
motif_wt <- FindMotifs(object = RZ.RZK, features = sig_peaks_wt)
motif_wt$logP <- -log10(motif_wt$p.adjust)

sig_motif_wt <- motif_wt[motif_wt$p.adjust < config$motif$enriched_padj, ]
cat("  Significant motifs (Wild-type):", nrow(sig_motif_wt), "\n")

write.csv(motif_wt, "/output/07_motif_enrichment_wildtype.csv", row.names = FALSE)

# Plot top 15 enriched motifs (Wild-type)
data_wt <- motif_wt[1:min(15, nrow(motif_wt)), ]
y_order_wt <- rev(data_wt$motif.name)

p_mot_wt <- ggplot(data_wt, aes(x = logP, y = motif.name, fill = percent.observed)) +
  geom_bar(stat = "identity", width = 0.2) +
  scale_fill_gradientn(colours = c("royalblue", "rosybrown2", "red"), limits = c(0, 60)) +
  scale_y_discrete(limits = y_order_wt) +
  xlim(0, 30) +
  ggtitle("Wnt7b+ vs Lgr5+\n(RZ / Wild-type)") +
  xlab("-log10(p.adjust)") + ylab("TF motif") +
  theme_linedraw() +
  theme(
    plot.title   = element_text(hjust = 0.5, color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(color = "black", size = 12),
    axis.title.y = element_text(color = "black", size = 12)
  )

ggsave(p_mot_wt,
       filename = file.path(fig_dir, "fig3I_motif_enrichment_wildtype.pdf"),
       width = 6, height = 6, dpi = 300)

cat("\n=== Step 07 complete ===\n")