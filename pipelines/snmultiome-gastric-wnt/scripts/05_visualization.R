#!/usr/bin/env Rscript
# ===========================================================================
#  05_visualization.R — Publication-quality figures (Fig 3A–D, S3A–B)
#  Replicates all visualization functions from the original code:
#    - UMAP colored by genotype + cluster (Fig 3A)
#    - DotPlot of marker genes (Fig 3B)
#    - dittoBarPlot of cluster proportions (Fig 3C)
#    - Wnt7b FeaturePlot split by genotype (Fig 3D)
#    - Individual marker FeaturePlots (Fig S3A)
#    - Wnt family VlnPlot (Fig S3B)
# ===========================================================================

suppressPackageStartupMessages({
    library(yaml)
    library(Seurat)
    library(Signac)
    library(ggplot2)
    library(RColorBrewer)
    library(dittoSeq)
    library(paletteer)
})

# ── Load configuration ────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
config_path <- args[which(args == "--config") + 1]
cfg <- read_yaml(config_path)

set.seed(cfg$seed)
output_dir <- cfg$output_dir
plot_dir   <- file.path(output_dir, "plots")

cat("=== Step 5: Visualization ===\n")

# ── Load integrated object ────────────────────────────────────────────────
obj <- readRDS(file.path(output_dir, "rds", "integrated_clustered.rds"))

# ── Color palettes (from original code) ───────────────────────────────────
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
    "Neck"          = "#8936EF",
    "Lgr5+"         = "#F2CA19",
    "SPEM"          = "#FF00BD",
    "Wnt7+"         = "#E11845",
    "Proliferating" = "#0057E9",
    "Pre-Pit"       = "#87E911",
    "Pit"           = "#018300"
)

mylevel <- c("Lgr5+", "SPEM", "Neck", "Proliferating", "Pre-Pit", "Pit", "Wnt7+")

# ── Helper: customized FeaturePlot (from plot_featureplot_v3) ─────────────
plot_featureplot_v3 <- function(obj, gene, min_cut, max_cut) {
    DefaultAssay(obj) <- "RNA"
    ht_custom_col <- rev(paletteer_c("grDevices::Spectral", 30))
    p <- FeaturePlot(obj, features = gene, pt.size = 1,
                     min.cutoff = min_cut, max.cutoff = max_cut, order = TRUE)
    ggplot(data = p$data, aes_string(x = "UMAP_1", y = "UMAP_2", fill = gene)) +
        geom_point(shape = 21, stroke = NA, color = "black", alpha = 1, size = 2) +
        scale_fill_gradientn(colours = ht_custom_col, name = "Expression\nLevel") +
        ggtitle(gene) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, size = 16))
}

# ── Helper: Wnt7b FeaturePlot with YlOrRd palette (from plot_featureplot) ─
plot_featureplot <- function(obj, gene) {
    DefaultAssay(obj) <- "RNA"
    p <- FeaturePlot(obj, features = gene, pt.size = 1,
                     min.cutoff = 0.3, max.cutoff = 2)
    ggplot(data = p$data, aes_string(x = "UMAP_1", y = "UMAP_2", fill = gene)) +
        geom_point(shape = 21, stroke = 0.3, color = "black", alpha = 1, size = 2) +
        scale_fill_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd")) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5))
}

# ══════════════════════════════════════════════════════════════════════════
#  Fig 3A — UMAP colored by cl_genotype
# ══════════════════════════════════════════════════════════════════════════
cat("  Fig 3A: UMAP by cluster × genotype...\n")
pdf(file.path(plot_dir, "fig3A_UMAP.pdf"), width = 8, height = 9)
p1 <- DimPlot(obj, group.by = "cl_genotype", cols = my_col, pt.size = 1.3) +
    NoLegend()
print(p1)

# Also show genotype split + cluster labels side by side
obj$kras_genotype <- factor(obj$kras_genotype, levels = c("Wild_type", "Mutant"))
cols6 <- c("#8936EF", "#F2CA19", "#FF00BD", "#E11845", "#0057E9", "#87E911")
p2a <- DimPlot(obj, reduction = "umap", group.by = "kras_genotype")
p2b <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE, cols = cols6)
print(p2a + p2b)
dev.off()

# ══════════════════════════════════════════════════════════════════════════
#  Fig 3B — DotPlot of marker genes
# ══════════════════════════════════════════════════════════════════════════
cat("  Fig 3B: DotPlot of marker genes...\n")
DefaultAssay(obj) <- "RNA"
Idents(obj) <- obj@meta.data$clusters
Idents(obj) <- factor(Idents(obj), levels = mylevel)
marker_genes <- cfg$marker_genes

pdf(file.path(plot_dir, "fig3B_dotplot.pdf"), width = 6, height = 7)
p2 <- DotPlot(obj, features = marker_genes, dot.scale = 10) +
    scale_colour_gradient2(low = "dodgerblue3", mid = "ghostwhite",
                           high = "firebrick", limits = c(-2.5, 2.5)) +
    xlab("Markers") + ylab("Clusters") + coord_flip() +
    theme(axis.text = element_text(size = 15), legend.position = "top")
print(p2)
dev.off()

# ══════════════════════════════════════════════════════════════════════════
#  Fig 3C — dittoBarPlot of cluster proportions by genotype
# ══════════════════════════════════════════════════════════════════════════
cat("  Fig 3C: dittoBarPlot...\n")
Idents(obj) <- obj@meta.data$clusters
Idents(obj) <- factor(Idents(obj), levels = mylevel)

pdf(file.path(plot_dir, "fig3C_dittobarplot.pdf"), width = 6, height = 6)
p3 <- dittoBarPlot(obj, Idents(obj), group.by = "kras_genotype",
                    x.reorder = c(2, 1),
                    var.labels.reorder = c(1, 6, 2, 5, 4, 3, 7),
                    color.panel = cluster_col)
print(p3)
dev.off()

# ══════════════════════════════════════════════════════════════════════════
#  Fig 3D — Wnt7b FeaturePlot split by genotype
# ══════════════════════════════════════════════════════════════════════════
cat("  Fig 3D: Wnt7b FeaturePlot (split by genotype)...\n")
pdf(file.path(plot_dir, "fig3D_Wnt7b_feature.pdf"), width = 12, height = 6)
p_wt  <- plot_featureplot(subset(obj, subset = kras_genotype == "Wild_type"), "Wnt7b")
p_mut <- plot_featureplot(subset(obj, subset = kras_genotype == "Mutant"),    "Wnt7b")
print(p_wt + p_mut)
dev.off()

# ══════════════════════════════════════════════════════════════════════════
#  Fig S3A — Individual marker gene FeaturePlots
# ══════════════════════════════════════════════════════════════════════════
cat("  Fig S3A: Marker gene FeaturePlots...\n")

# Exact min/max cutoffs from the original figure code
feature_params <- list(
    list("Muc5ac", 0, 3),   list("Gkn2", 0, 3),
    list("Tff1", 0.5, 6),   list("Gkn1", 0.5, 5),
    list("Smc2", 0.5, 3),   list("Top2a", 0.5, 4),
    list("Hmgb2", 1, 4),    list("Foxm1", 0, 2),
    list("Mki67", 0.3, 3),  list("Muc6", 0.3, 4.5),
    list("Cftr", 0.5, 3.5), list("Cd44", 2, 4.2),
    list("Glipr1", 1.5, 4.5), list("Lgr5", 0.5, 3.5)
)

pdf(file.path(plot_dir, "figS3A_featureplots.pdf"), width = 10, height = 10)
for (params in feature_params) {
    p <- plot_featureplot_v3(obj, params[[1]], params[[2]], params[[3]])
    print(p)
}
dev.off()

# ══════════════════════════════════════════════════════════════════════════
#  Fig S3B — Wnt family VlnPlot by genotype
# ══════════════════════════════════════════════════════════════════════════
cat("  Fig S3B: Wnt gene family VlnPlot...\n")
DefaultAssay(obj) <- "RNA"
wnt_genes <- cfg$wnt_genes

pdf(file.path(plot_dir, "figS3B_wnt_family_violin.pdf"), width = 16, height = 12)
p_wnt <- VlnPlot(object = obj, features = wnt_genes, group.by = "kras_genotype")
print(p_wnt)
dev.off()

cat("=== Step 5 complete ===\n")