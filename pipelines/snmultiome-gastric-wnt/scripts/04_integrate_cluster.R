#!/usr/bin/env Rscript
# ===========================================================================
#  04_integrate_cluster.R — SCT-based integration + UMAP + clustering
#  Replicates: integrate_based_on_rna() + generate_umap_integrated()
#              + cluster refinement (res 0.6) + cell-type annotation
#
#  Integration (SCT anchors):
#    1. SelectIntegrationFeatures → PrepSCTIntegration
#    2. FindIntegrationAnchors (normalization.method = "SCT")
#    3. IntegrateData (normalization.method = "SCT")
#
#  UMAP + Clustering:
#    1. ScaleData on integrated assay
#    2. RunPCA (npcs = 30)
#    3. RunUMAP (dims 1:30)
#    4. FindNeighbors + FindClusters (resolution = 0.3 for initial 6 clusters)
#    5. Re-cluster at resolution = 0.6 to split Pit into Pre-Pit + Pit
#    6. Assign cell-type labels per cluster + genotype metadata
# ===========================================================================

suppressPackageStartupMessages({
    library(yaml)
    library(future)
    library(Seurat)
    library(Signac)
})

# ── Load configuration ────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
config_path <- args[which(args == "--config") + 1]
cfg <- read_yaml(config_path)

set.seed(cfg$seed)
plan("multicore", workers = cfg$workers)
options(future.globals.maxSize = cfg$future_ram_gb * 1024^3)
output_dir <- cfg$output_dir

cat("=== Step 4: Integrate & cluster ===\n")

# ── Load normalized objects ───────────────────────────────────────────────
wild   <- readRDS(file.path(output_dir, "rds", "wild_normalized.rds"))
mutant <- readRDS(file.path(output_dir, "rds", "mutant_normalized.rds"))

# ── SCT-based integration (integrate_based_on_rna) ───────────────────────
cat("Integrating wild-type and mutant via SCT anchors...\n")
DefaultAssay(wild)   <- "SCT"
DefaultAssay(mutant) <- "SCT"
obj.list <- list(wild, mutant)

features    <- SelectIntegrationFeatures(obj.list)
obj.list    <- PrepSCTIntegration(obj.list, anchor.features = features)
obj.anchors <- FindIntegrationAnchors(
    object.list          = obj.list,
    anchor.features      = features,
    normalization.method = "SCT"
)
obj.combined <- IntegrateData(
    anchorset            = obj.anchors,
    normalization.method = "SCT"
)
cat(sprintf("Integrated object: %d cells\n", ncol(obj.combined)))

# Free memory
rm(wild, mutant, obj.list, obj.anchors); gc()

# ── UMAP + initial clustering at resolution 0.3 ──────────────────────────
cat("Running PCA, UMAP, and clustering (res=0.3)...\n")
DefaultAssay(obj.combined) <- "integrated"
obj.combined <- ScaleData(obj.combined, verbose = FALSE)
obj.combined <- RunPCA(obj.combined, npcs = cfg$pca_dims, verbose = FALSE)
obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:cfg$pca_dims)
obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 1:cfg$pca_dims)
obj.combined <- FindClusters(obj.combined,
                             resolution = cfg$clustering_resolution_initial,
                             algorithm  = 1)

# Create cl_genotype metadata (cluster_genotype combination)
obj.combined@meta.data$cl_genotype <- paste0(
    obj.combined@meta.data$integrated_snn_res.0.3, "_",
    obj.combined@meta.data$kras_genotype
)
cat(sprintf("Initial clusters (res=0.3): %d\n",
            length(unique(obj.combined$seurat_clusters))))

# ── Refined clustering at resolution 0.6 (splits Pit cells) ──────────────
# As in the original figure code: FindClusters at 0.6, then map
# cluster 9 (from the higher resolution) into a new cluster "6"
cat("Refining clusters at resolution 0.6...\n")
obj.combined <- FindClusters(obj.combined,
                              resolution = cfg$clustering_resolution_refine)

# Map cluster 9 at res=0.6 to a new "6" genotype group (Pit cells)
obj.combined@meta.data$cl_genotype[
    obj.combined@meta.data$seurat_clusters == "9" &
    obj.combined@meta.data$kras_genotype == "Wild_type"
] <- "6_Wild_type"
obj.combined@meta.data$cl_genotype[
    obj.combined@meta.data$seurat_clusters == "9" &
    obj.combined@meta.data$kras_genotype == "Mutant"
] <- "6_Mutant"

# ── Assign cell-type labels ───────────────────────────────────────────────
# Mapping from config: 0=Neck, 1=Lgr5+, 2=SPEM, 3=Wnt7+, 4=Proliferating,
#                       5=Pre-Pit, 6=Pit
cat("Assigning cell-type labels...\n")
cluster_labels <- cfg$cluster_labels
for (cl_num in names(cluster_labels)) {
    label <- cluster_labels[[cl_num]]
    for (geno in c("Wild_type", "Mutant")) {
        cl_geno <- paste0(cl_num, "_", geno)
        obj.combined@meta.data$clusters[
            obj.combined@meta.data$cl_genotype == cl_geno
        ] <- label
    }
}

# Log cluster composition
cat("Cluster composition:\n")
print(table(obj.combined$clusters, obj.combined$kras_genotype))

# ── Save integrated + clustered object ────────────────────────────────────
saveRDS(obj.combined, file.path(output_dir, "rds", "integrated_clustered.rds"))
cat("=== Step 4 complete ===\n")