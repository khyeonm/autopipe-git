# ==============================================================================
# 04_integrate_cluster.R — SCT integration, PCA, UMAP, clustering
# ==============================================================================
# Integrates wild-type and mutant objects on the SCT assay, then runs
# PCA → UMAP → FindNeighbors → FindClusters at two resolutions.
# Assigns cell-type labels and cl_genotype metadata.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(future)
  library(yaml)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
cfg  <- read_yaml(args[1])

plan("multicore", workers = cfg$threads)
options(future.globals.maxSize = cfg$max_memory_gb * 1024^3)
set.seed(1234)

outdir <- cfg$output_dir

# ---- Load normalised objects ------------------------------------------------
cat(">> Loading normalised objects\n")
wild   <- readRDS(file.path(outdir, "03_normalised", "wild_type_normalised.rds"))
mutant <- readRDS(file.path(outdir, "03_normalised", "mutant_normalised.rds"))

# ---- SCT-based integration -------------------------------------------------
cat(">> Integrating samples via SCT anchors\n")
DefaultAssay(wild)   <- "SCT"
DefaultAssay(mutant) <- "SCT"

obj.list <- list(wild, mutant)
features <- SelectIntegrationFeatures(obj.list)
obj.list <- PrepSCTIntegration(obj.list, anchor.features = features)
obj.anchors <- FindIntegrationAnchors(
  object.list          = obj.list,
  anchor.features      = features,
  normalization.method  = "SCT"
)
combined <- IntegrateData(
  anchorset            = obj.anchors,
  normalization.method  = "SCT"
)

# ---- PCA, UMAP, Clustering (initial resolution) ----------------------------
cat(">> Running PCA, UMAP, and clustering\n")
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = cfg$n_pcs, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:cfg$n_pcs)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:cfg$n_pcs)

# Initial clustering at resolution 0.3
combined <- FindClusters(combined, resolution = cfg$cluster_resolution_initial, algorithm = 1)

# Create cl_genotype metadata (cluster_genotype)
combined@meta.data$cl_genotype <- paste0(
  combined@meta.data$integrated_snn_res.0.3, "_",
  combined@meta.data$kras_genotype
)

# ---- Refined clustering (resolution 0.6) to split Pit / Pre-Pit ------------
cat(">> Refined clustering at resolution", cfg$cluster_resolution_refined, "\n")
combined <- FindClusters(combined, resolution = cfg$cluster_resolution_refined)

# Map cluster 9 (from res 0.6) → cluster 6 in cl_genotype (the Pit cluster)
combined@meta.data$cl_genotype[
  combined@meta.data$seurat_clusters == "9" &
  combined@meta.data$kras_genotype == "Wild_type"
] <- "6_Wild_type"

combined@meta.data$cl_genotype[
  combined@meta.data$seurat_clusters == "9" &
  combined@meta.data$kras_genotype == "Mutant"
] <- "6_Mutant"

# ---- Assign human-readable cell-type labels ---------------------------------
cat(">> Assigning cell-type labels\n")
cluster_labels <- cfg$cluster_labels  # mapping: "0" → "Neck", etc.

combined@meta.data$clusters <- NA
for (cl_num in names(cluster_labels)) {
  label <- cluster_labels[[cl_num]]
  for (geno in c("Wild_type", "Mutant")) {
    cl_geno <- paste0(cl_num, "_", geno)
    combined@meta.data$clusters[combined@meta.data$cl_genotype == cl_geno] <- label
  }
}

# Factor genotype
combined$kras_genotype <- factor(combined$kras_genotype,
                                  levels = c("Wild_type", "Mutant"))

# ---- Save -------------------------------------------------------------------
dir.create(file.path(outdir, "04_integrated"), recursive = TRUE, showWarnings = FALSE)
saveRDS(combined, file.path(outdir, "04_integrated", "RZ_RZK_integrated.rds"))

cat(">> Step 04 complete. Total cells:", ncol(combined), "\n")
cat("   Cluster distribution:\n")
print(table(combined@meta.data$clusters, combined@meta.data$kras_genotype))