#!/usr/bin/env Rscript
# =============================================================================
#  04_integrate_cluster.R — SCT-based integration + UMAP clustering
#
#  Integration strategy (from the original code):
#    1. SelectIntegrationFeatures on the two SCT-normalized objects
#    2. PrepSCTIntegration → FindIntegrationAnchors (SCT method)
#    3. IntegrateData → creates "integrated" assay
#    4. ScaleData → PCA (30 PCs) → UMAP → FindNeighbors → FindClusters
#       Initial clustering at resolution 0.3 → 6 clusters (0–5)
#    5. Refined clustering at resolution 0.6 to split Pit into Pre-Pit / Pit
#       (cluster 9 at res 0.6 maps to cluster 6 = "Pit")
#    6. Assign cell-type labels based on marker expression patterns
#    7. Create cl_genotype metadata column (cluster_genotype combinations)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(yaml)
  library(future)
})

config <- yaml::read_yaml("/pipeline/config.yaml")
plan("multicore", workers = config$workers)
options(future.globals.maxSize = config$max_ram_gb * 1024^3)
set.seed(config$seed)

cat("=== Step 04: Integration & Clustering ===\n")

# =============================================================================
#  1. Load normalized objects
# =============================================================================
wild   <- readRDS("/output/03_wild_normalized.rds")
mutant <- readRDS("/output/03_mutant_normalized.rds")

cat("Wild-type cells:", ncol(wild), "\n")
cat("Mutant cells:   ", ncol(mutant), "\n")

# =============================================================================
#  2. SCT-based integration
# =============================================================================
# Seurat's SCT integration workflow:
#   - SelectIntegrationFeatures picks genes variable across both datasets
#   - PrepSCTIntegration prepares SCT residuals for anchor finding
#   - FindIntegrationAnchors identifies mutual nearest neighbors
#   - IntegrateData corrects batch effects into the "integrated" assay
cat("  Selecting integration features...\n")
DefaultAssay(wild)   <- "SCT"
DefaultAssay(mutant) <- "SCT"
obj.list <- list(wild, mutant)

features <- SelectIntegrationFeatures(obj.list)
cat("  Integration features:", length(features), "\n")

cat("  Preparing SCT integration...\n")
obj.list <- PrepSCTIntegration(obj.list, anchor.features = features)

cat("  Finding integration anchors...\n")
obj.anchors <- FindIntegrationAnchors(
  object.list          = obj.list,
  anchor.features      = features,
  normalization.method = "SCT"
)

cat("  Integrating data...\n")
obj.combined <- IntegrateData(
  anchorset            = obj.anchors,
  normalization.method = "SCT"
)

rm(wild, mutant, obj.list, obj.anchors); gc()

# =============================================================================
#  3. Dimensionality reduction + initial clustering (resolution 0.3)
# =============================================================================
cat("  Scaling + PCA + UMAP + Clustering...\n")
DefaultAssay(obj.combined) <- "integrated"
obj.combined <- ScaleData(obj.combined, verbose = FALSE)
obj.combined <- RunPCA(obj.combined, npcs = config$clustering$npcs, verbose = FALSE)
obj.combined <- RunUMAP(
  obj.combined,
  reduction = "pca",
  dims = config$clustering$dims_start:config$clustering$dims_end
)
obj.combined <- FindNeighbors(
  obj.combined,
  reduction = "pca",
  dims = config$clustering$dims_start:config$clustering$dims_end
)
obj.combined <- FindClusters(
  obj.combined,
  resolution = config$clustering$resolution_initial,
  algorithm  = config$clustering$algorithm
)

cat("  Clusters at res", config$clustering$resolution_initial, ":",
    length(unique(Idents(obj.combined))), "\n")

# Create cl_genotype column: "cluster_genotype" (e.g., "0_Wild_type")
obj.combined@meta.data$cl_genotype <- paste0(
  obj.combined@meta.data$integrated_snn_res.0.3, "_",
  obj.combined@meta.data$kras_genotype
)

# =============================================================================
#  4. Refined clustering (resolution 0.6) to split Pit cells
# =============================================================================
# At resolution 0.3, cluster 5 = "Pit" includes both Pre-Pit and mature Pit.
# At resolution 0.6, cluster 9 separates out as a distinct Pit population.
# We re-map cluster 9 → cluster 6 and assign cell-type labels accordingly.
cat("  Refined clustering at res", config$clustering$resolution_refined, "...\n")
obj.combined <- FindClusters(
  obj.combined,
  resolution = config$clustering$resolution_refined
)

# Map cluster 9 (from res 0.6) to a new "6" identity for Pit cells
obj.combined@meta.data$cl_genotype[
  obj.combined@meta.data$seurat_clusters == "9" &
  obj.combined@meta.data$kras_genotype == "Wild_type"
] <- "6_Wild_type"

obj.combined@meta.data$cl_genotype[
  obj.combined@meta.data$seurat_clusters == "9" &
  obj.combined@meta.data$kras_genotype == "Mutant"
] <- "6_Mutant"

# =============================================================================
#  5. Assign cell-type labels
# =============================================================================
# Cell type identification based on marker gene expression:
#   Cluster 0 → Neck (Muc6, Cftr)
#   Cluster 1 → Lgr5+ stem cells (Lgr5)
#   Cluster 2 → SPEM (Glipr1)
#   Cluster 3 → Wnt7+ (Wnt7b)
#   Cluster 4 → Proliferating (Mki67, Stmn1, Top2a, Foxm1)
#   Cluster 5 → Pre-Pit (Gkn1, Tff1)
#   Cluster 6 → Pit (Muc5ac, Gkn2)
cell_types <- config$cell_types

for (cl_num in names(cell_types)) {
  label <- cell_types[[cl_num]]
  for (geno in c("Wild_type", "Mutant")) {
    cl_geno <- paste0(cl_num, "_", geno)
    obj.combined@meta.data$clusters[
      obj.combined@meta.data$cl_genotype == cl_geno
    ] <- label
  }
}

cat("  Cell types assigned:\n")
print(table(obj.combined@meta.data$clusters))

# =============================================================================
#  6. Save integrated + clustered object
# =============================================================================
saveRDS(obj.combined, "/output/04_integrated_clustered.rds")

cat("\n=== Step 04 complete ===\n")
cat("Total cells:", ncol(obj.combined), "\n")
cat("Clusters:   ", paste(unique(obj.combined@meta.data$clusters), collapse = ", "), "\n")