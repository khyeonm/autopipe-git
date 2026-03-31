# ==============================================================================
# 05_visualization.R — Generate publication figures (Fig 3A–D, S3A–B)
# ==============================================================================
# Reproduces the key visualisation panels from the paper:
#   Fig 3A: UMAP coloured by cl_genotype
#   Fig 3B: DotPlot of marker genes across cell types
#   Fig 3C: dittoBarPlot of cell-type proportions by genotype
#   Fig 3D: Wnt7b FeaturePlot split by genotype
#   Fig S3A: FeaturePlots for all marker genes
#   Fig S3B: VlnPlot of Wnt family genes by genotype
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(RColorBrewer)
  library(dittoSeq)
  library(yaml)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
cfg  <- read_yaml(args[1])
outdir <- cfg$output_dir

# ---- Load integrated object -------------------------------------------------
cat(">> Loading integrated object\n")
obj <- readRDS(file.path(outdir, "04_integrated", "RZ_RZK_integrated.rds"))

# ---- Colour palettes from the original scripts ------------------------------
my_col <- c(
  "0_Wild_type" = "#8936EF", "0_Mutant" = "#8936EF",
  "1_Wild_type" = "#F2CA19", "1_Mutant" = "#F2CA19",
  "2_Wild_type" = "#FF00BD", "2_Mutant" = "#FF00BD",
  "3_Wild_type" = "#E11845", "3_Mutant" = "#E11845",
  "4_Wild_type" = "#0057E9", "4_Mutant" = "#0057E9",
  "5_Wild_type" = "#87E911", "5_Mutant" = "#87E911",
  "6_Wild_type" = "#018300", "6_Mutant" = "#018300"
)
cluster_col <- c(
  "Neck" = "#8936EF", "Lgr5+" = "#F2CA19", "SPEM" = "#FF00BD",
  "Wnt7+" = "#E11845", "Proliferating" = "#0057E9",
  "Pre-Pit" = "#87E911", "Pit" = "#018300"
)
mylevel <- c("Lgr5+", "SPEM", "Neck", "Proliferating", "Pre-Pit", "Pit", "Wnt7+")

# ---- Output directory -------------------------------------------------------
figdir <- file.path(outdir, "05_figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# Fig 3A — UMAP coloured by cluster × genotype
# ===========================================================================
cat(">> Generating Fig 3A\n")
p3a <- DimPlot(obj, group.by = "cl_genotype", cols = my_col, pt.size = 1.3) + NoLegend()
ggsave(file.path(figdir, "fig3A_UMAP.pdf"), plot = p3a, width = 8, height = 9)

# ===========================================================================
# Fig 3B — DotPlot of marker genes
# ===========================================================================
cat(">> Generating Fig 3B\n")
DefaultAssay(obj) <- "RNA"
Idents(obj) <- factor(obj@meta.data$clusters, levels = mylevel)
marker_genes <- cfg$marker_genes

p3b <- DotPlot(obj, features = marker_genes, dot.scale = 10) +
  scale_colour_gradient2(low = "dodgerblue3", mid = "ghostwhite", high = "firebrick",
                         limits = c(-2.5, 2.5)) +
  xlab("Markers") + ylab("Clusters") + coord_flip() +
  theme(axis.text = element_text(size = 15), legend.position = "top")
ggsave(file.path(figdir, "fig3B_dotplot.pdf"), plot = p3b, width = 6, height = 7)

# ===========================================================================
# Fig 3C — dittoBarPlot (cell-type proportions by genotype)
# ===========================================================================
cat(">> Generating Fig 3C\n")
Idents(obj) <- factor(obj@meta.data$clusters, levels = mylevel)
p3c <- dittoBarPlot(obj, Idents(obj), group.by = "kras_genotype",
                    x.reorder = c(2, 1),
                    var.labels.reorder = c(1, 6, 2, 5, 4, 3, 7),
                    color.panel = cluster_col)
ggsave(file.path(figdir, "fig3C_dittobarplot.pdf"), plot = p3c, width = 6, height = 6)

# ===========================================================================
# Fig 3D — Wnt7b FeaturePlot (split by genotype)
# ===========================================================================
cat(">> Generating Fig 3D\n")

# Helper: custom FeaturePlot using modern aes() syntax
plot_featureplot <- function(obj, gene) {
  DefaultAssay(obj) <- "RNA"
  p <- FeaturePlot(obj, features = gene, pt.size = 1,
                   min.cutoff = 0.3, max.cutoff = 2)
  # Extract data and detect column names dynamically
  pdata <- p$data
  umap_cols <- grep("^[Uu][Mm][Aa][Pp]", colnames(pdata), value = TRUE)
  umap1 <- umap_cols[1]
  umap2 <- umap_cols[2]

  ggplot(data = pdata, aes(x = .data[[umap1]], y = .data[[umap2]], fill = .data[[gene]])) +
    geom_point(shape = 21, stroke = 0.3, color = "black", alpha = 1, size = 2) +
    scale_fill_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd")) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}

p_wt  <- plot_featureplot(subset(obj, subset = kras_genotype == "Wild_type"), "Wnt7b")
p_mut <- plot_featureplot(subset(obj, subset = kras_genotype == "Mutant"), "Wnt7b")
p3d   <- p_wt + p_mut
ggsave(file.path(figdir, "fig3D_Wnt7b_featureplot.pdf"), plot = p3d, width = 12, height = 6)

# ===========================================================================
# Fig S3A — FeaturePlots for all marker genes (Spectral colour scale)
# ===========================================================================
cat(">> Generating Fig S3A\n")
plot_featureplot_v3 <- function(obj, gene, min_val, max_val) {
  DefaultAssay(obj) <- "RNA"
  ht_custom_col <- rev(grDevices::hcl.colors(30, palette = "Spectral"))
  p <- FeaturePlot(obj, features = gene, pt.size = 1,
                   min.cutoff = min_val, max.cutoff = max_val, order = TRUE)
  pdata <- p$data
  umap_cols <- grep("^[Uu][Mm][Aa][Pp]", colnames(pdata), value = TRUE)
  umap1 <- umap_cols[1]
  umap2 <- umap_cols[2]

  ggplot(data = pdata, aes(x = .data[[umap1]], y = .data[[umap2]], fill = .data[[gene]])) +
    geom_point(shape = 21, stroke = NA, color = "black", alpha = 1, size = 2) +
    scale_fill_gradientn(colours = ht_custom_col, name = "Expression\nLevel") +
    ggtitle(gene) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 16))
}

# Min/max cutoffs extracted from the original 2_RZRZK_figure.R
marker_params <- list(
  list("Muc5ac", 0, 3),   list("Gkn2", 0, 3),      list("Tff1", 0.5, 6),
  list("Gkn1", 0.5, 5),   list("Smc2", 0.5, 3),     list("Top2a", 0.5, 4),
  list("Hmgb2", 1, 4),    list("Foxm1", 0, 2),       list("Mki67", 0.3, 3),
  list("Muc6", 0.3, 4.5), list("Cftr", 0.5, 3.5),    list("Cd44", 2, 4.2),
  list("Glipr1", 1.5, 4.5), list("Lgr5", 0.5, 3.5)
)

pdf(file.path(figdir, "figS3A_marker_featureplots.pdf"), width = 10, height = 10)
for (mp in marker_params) {
  p <- plot_featureplot_v3(obj, mp[[1]], mp[[2]], mp[[3]])
  print(p)
}
dev.off()

# ===========================================================================
# Fig S3B — Wnt family VlnPlot by genotype
# ===========================================================================
cat(">> Generating Fig S3B\n")
DefaultAssay(obj) <- "RNA"
wnt_genes <- cfg$wnt_family_genes

pdf(file.path(figdir, "figS3B_wnt_family_violin.pdf"), width = 16, height = 12)
print(VlnPlot(obj, features = wnt_genes, group.by = "kras_genotype"))
dev.off()

cat(">> Step 05 complete.\n")