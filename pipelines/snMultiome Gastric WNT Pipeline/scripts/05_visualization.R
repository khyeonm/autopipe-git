#!/usr/bin/env Rscript
# =============================================================================
#  05_visualization.R — Generate publication figures (Fig 3A–H, S3A–B)
#
#  Produces the following plots from the integrated snMultiome object:
#    Fig 3A: UMAP colored by cl_genotype
#    Fig 3B: DotPlot of marker genes across cell types
#    Fig 3C: dittoBarPlot of cell-type proportions by genotype
#    Fig 3D: Wnt7b FeaturePlot split by genotype
#    Fig 3H: CoveragePlot of Wnt7b locus linked to expression
#    Fig S3A: Individual FeaturePlots for all marker genes
#    Fig S3B: VlnPlot of all Wnt family genes by genotype
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(dittoSeq)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(yaml)
  library(future)
})

config <- yaml::read_yaml("/pipeline/config.yaml")
plan("multicore", workers = config$workers)
options(future.globals.maxSize = config$max_ram_gb * 1024^3)
set.seed(config$seed)

cat("=== Step 05: Visualization ===\n")

# Create output directory for figures
fig_dir <- "/output/05_figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Load integrated object
RZ.RZK <- readRDS("/output/04_integrated_clustered.rds")

# ---- Color definitions (matching the original code) ----
my_col <- c(
  "0_Wild_type" = "#8936EF", "0_Mutant" = "#8936EF",
  "1_Wild_type" = "#F2CA19", "1_Mutant" = "#F2CA19",
  "2_Wild_type" = "#FF00BD", "2_Mutant" = "#FF00BD",
  "3_Wild_type" = "#E11845", "3_Mutant" = "#E11845",
  "4_Wild_type" = "#0057E9", "4_Mutant" = "#0057E9",
  "5_Wild_type" = "#87E911", "5_Mutant" = "#87E911",
  "6_Wild_type" = "#018300", "6_Mutant" = "#018300"
)

cluster_col <- unlist(config$cluster_colors)
mylevel <- c("Lgr5+", "SPEM", "Neck", "Proliferating", "Pre-Pit", "Pit", "Wnt7+")

# =============================================================================
#  Custom FeaturePlot function (from the original 2_RZRZK_figure.R)
# =============================================================================
# Extracts UMAP coordinates and expression from Seurat's FeaturePlot data,
# then re-renders with ggplot2 using YlOrRd color scale and point styling.
plot_featureplot <- function(obj, gene) {
  DefaultAssay(obj) <- "RNA"
  p11 <- FeaturePlot(obj, features = gene, pt.size = 1,
                     min.cutoff = 0.3, max.cutoff = 2)
  ggplot(data = p11$data, aes_string(x = "UMAP_1", y = "UMAP_2", fill = gene)) +
    geom_point(shape = 21, stroke = 0.3, color = "black", alpha = 1, size = 2) +
    scale_fill_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd")) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}

# FeaturePlot with configurable min/max cutoffs and Spectral palette
plot_featureplot_v3 <- function(obj, gene, min_val, max_val) {
  DefaultAssay(obj) <- "RNA"
  # Use a reversed Spectral palette for high contrast
  ht_custom_col <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(30))
  p11 <- FeaturePlot(obj, features = gene, pt.size = 1,
                     min.cutoff = min_val, max.cutoff = max_val, order = TRUE)
  ggplot(data = p11$data, aes_string(x = "UMAP_1", y = "UMAP_2", fill = gene)) +
    geom_point(shape = 21, stroke = NA, color = "black", alpha = 1, size = 2) +
    scale_fill_gradientn(colours = ht_custom_col, name = "Expression\nLevel") +
    ggtitle(gene) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 16))
}

# =============================================================================
#  Fig 3A: UMAP colored by cluster × genotype
# =============================================================================
cat("  Fig 3A: UMAP by cl_genotype...\n")
p1 <- DimPlot(RZ.RZK, group.by = "cl_genotype", cols = my_col, pt.size = 1.3) +
  NoLegend()
ggsave(p1, filename = file.path(fig_dir, "fig3A_UMAP.pdf"),
       width = 8, height = 9)

# =============================================================================
#  Fig 3B: DotPlot of marker genes
# =============================================================================
cat("  Fig 3B: DotPlot of marker genes...\n")
DefaultAssay(RZ.RZK) <- "RNA"
Idents(RZ.RZK) <- RZ.RZK@meta.data$clusters
Idents(RZ.RZK) <- factor(Idents(RZ.RZK), levels = mylevel)

myGenes <- config$marker_genes

p2 <- DotPlot(RZ.RZK, features = myGenes, dot.scale = 10) +
  scale_colour_gradient2(
    low = "dodgerblue3", mid = "ghostwhite", high = "firebrick",
    limits = c(-2.5, 2.5)
  ) +
  xlab("Markers") + ylab("Clusters") +
  coord_flip() +
  theme(axis.text = element_text(size = 15), legend.position = "top")
ggsave(p2, filename = file.path(fig_dir, "fig3B_dotplot.pdf"),
       width = 6, height = 7)

# =============================================================================
#  Fig 3C: dittoBarPlot of cell-type proportions
# =============================================================================
cat("  Fig 3C: dittoBarPlot...\n")
Idents(RZ.RZK) <- factor(RZ.RZK@meta.data$clusters, levels = mylevel)
p3 <- dittoBarPlot(
  RZ.RZK, Idents(RZ.RZK),
  group.by = "kras_genotype",
  x.reorder = c(2, 1),
  var.labels.reorder = c(1, 6, 2, 5, 4, 3, 7),
  color.panel = cluster_col
)
ggsave(p3, filename = file.path(fig_dir, "fig3C_dittobar.pdf"),
       width = 6, height = 6)

# =============================================================================
#  Fig 3D: Wnt7b FeaturePlot split by genotype
# =============================================================================
cat("  Fig 3D: Wnt7b FeaturePlot...\n")
p1_wnt <- plot_featureplot(
  subset(RZ.RZK, subset = kras_genotype == "Wild_type"), "Wnt7b"
)
p2_wnt <- plot_featureplot(
  subset(RZ.RZK, subset = kras_genotype == "Mutant"), "Wnt7b"
)
p4 <- p1_wnt + p2_wnt
ggsave(p4, filename = file.path(fig_dir, "fig3D_Wnt7b_feature.pdf"),
       width = 12, height = 6, dpi = 300)

# =============================================================================
#  Fig S3A: Individual FeaturePlots for all markers
# =============================================================================
cat("  Fig S3A: Individual FeaturePlots...\n")
fp_genes <- config$feature_plot_genes
pdf(file.path(fig_dir, "figS3A_feature_plots.pdf"), width = 10, height = 10)
for (g in fp_genes) {
  p <- plot_featureplot_v3(RZ.RZK, g$gene, g$min, g$max)
  print(p)
}
dev.off()

# =============================================================================
#  Fig S3B: Wnt family VlnPlot by genotype
# =============================================================================
cat("  Fig S3B: Wnt family VlnPlot...\n")
DefaultAssay(RZ.RZK) <- "RNA"
wnt_genes <- config$wnt_genes
p_wnt_vln <- VlnPlot(
  object   = RZ.RZK,
  features = wnt_genes,
  group.by = "kras_genotype"
)
ggsave(p_wnt_vln, filename = file.path(fig_dir, "figS3B_wnt_violin.pdf"),
       width = 20, height = 15)

# =============================================================================
#  Fig 3H: CoveragePlot of Wnt7b locus with peak-gene linkage
# =============================================================================
cat("  Fig 3H: CoveragePlot (Wnt7b linkage)...\n")
linkage_genome <- BSgenome.Mmusculus.UCSC.mm10
target_gene <- config$linkage_target_gene

DefaultAssay(RZ.RZK) <- "macs2"

# RegionStats computes GC content for each peak (needed for LinkPeaks)
RZ.RZK <- RegionStats(RZ.RZK, genome = linkage_genome)

# LinkPeaks correlates peak accessibility with gene expression
RZ.RZK <- LinkPeaks(
  object           = RZ.RZK,
  peak.assay       = "macs2",
  expression.assay = "RNA",
  genes.use        = c(target_gene)
)

# Order cl_genotype levels for the coverage plot
lv.list <- c()
for (i in c(1, 2, 0, 4, 5, 6, 3)) {
  lv.list <- c(lv.list, sprintf("%s_Wild_type", i))
  lv.list <- c(lv.list, sprintf("%s_Mutant", i))
}
RZ.RZK$cl_genotype <- factor(RZ.RZK$cl_genotype, levels = lv.list)

p6 <- CoveragePlot(
  object              = RZ.RZK,
  region              = target_gene,
  features            = target_gene,
  expression.assay    = "RNA",
  group.by            = "cl_genotype",
  extend.upstream     = 0,
  extend.downstream   = 0
)
ggsave(p6, filename = file.path(fig_dir, "fig3H_coverage_Wnt7b.pdf"),
       width = 8, height = 8, dpi = 300)

# Save the object with linkage info for downstream steps
saveRDS(RZ.RZK, "/output/04_integrated_clustered.rds")

cat("\n=== Step 05 complete ===\n")
cat("Figures saved to:", fig_dir, "\n")