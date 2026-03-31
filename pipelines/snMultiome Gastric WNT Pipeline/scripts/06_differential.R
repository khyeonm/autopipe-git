# ==============================================================================
# 06_differential.R — Differential expression (RNA) & accessibility (ATAC)
# ==============================================================================
# 1) DEG: Mutant vs Wild_type (all cells, MAST test)
# 2) KRAS signaling overlap (Hallmark gene set via msigdbr + gprofiler2)
# 3) DA peaks: Wnt7+ vs Lgr5+ in mutant and wild-type (LR test)
# 4) Wnt7b expression statistics (Fisher's exact + Wilcoxon)
# 5) Wnt7b accessibility statistics (t-test on peak region counts)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(EnsDb.Mmusculus.v79)
  library(future)
  library(yaml)
  library(dplyr)
  library(msigdbr)
  library(gprofiler2)
})

args <- commandArgs(trailingOnly = TRUE)
cfg  <- read_yaml(args[1])

plan("multicore", workers = cfg$threads)
options(future.globals.maxSize = cfg$max_memory_gb * 1024^3)
set.seed(1234)

outdir <- cfg$output_dir
diffdir <- file.path(outdir, "06_differential")
dir.create(diffdir, recursive = TRUE, showWarnings = FALSE)

# ---- Load -------------------------------------------------------------------
cat(">> Loading integrated object\n")
obj <- readRDS(file.path(outdir, "04_integrated", "RZ_RZK_integrated.rds"))

# ---- Genome objects ---------------------------------------------------------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"
linkage_genome <- BSgenome.Mmusculus.UCSC.mm10

# ---- Join layers for RNA only (required by Seurat v5 before FindMarkers) ----
# Note: ChromatinAssay (ATAC/macs2) does NOT support JoinLayers — skip it.
cat(">> Joining layers for RNA assay\n")
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

# ==========================================================================
# 1) Differential gene expression: Mutant vs Wild_type (all cells, MAST)
# ==========================================================================
cat(">> DEG: Mutant vs Wild_type (MAST)\n")
DefaultAssay(obj) <- "RNA"
DEG <- FindMarkers(obj, group.by = "kras_genotype", ident.1 = "Mutant",
                   test.use = cfg$de_test)
write.csv(DEG, file.path(diffdir, "DEG_MAST_results.csv"))

# ==========================================================================
# 2) KRAS signaling overlap
# ==========================================================================
cat(">> KRAS signaling gene set overlap\n")
temp.deg <- DEG[DEG$p_val_adj < cfg$de_padj_cutoff &
                DEG$avg_log2FC > cfg$de_log2fc_cutoff, ]

hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
kras_signaling <- hallmark$gene_symbol[hallmark$gs_name == "HALLMARK_KRAS_SIGNALING_UP"]

# Convert mouse gene symbols to human via gprofiler2
if (nrow(temp.deg) > 0) {
  temp_sc_2 <- gconvert(query = rownames(temp.deg), organism = "hsapiens", target = "ENSG")
  overlap <- intersect(temp_sc_2$name, kras_signaling)
  kras_df <- data.frame(
    n_DEG_up = nrow(temp.deg),
    n_hallmark_kras_up = length(kras_signaling),
    n_overlap = length(overlap),
    overlap_genes = paste(overlap, collapse = ";")
  )
} else {
  kras_df <- data.frame(n_DEG_up = 0, n_hallmark_kras_up = length(kras_signaling),
                        n_overlap = 0, overlap_genes = "")
}
write.csv(kras_df, file.path(diffdir, "KRAS_signaling_overlap.csv"), row.names = FALSE)

# ==========================================================================
# 3) DA peaks: Wnt7+ (cluster 3) vs Lgr5+ (cluster 1)
# ==========================================================================
cat(">> DA peaks: Wnt7+ vs Lgr5+ (mutant)\n")
plan("multicore", workers = 1)   # LR test is not parallel-safe
DefaultAssay(obj) <- "ATAC"

da_peaks_mut <- FindMarkers(
  obj, ident.1 = "3_Mutant", ident.2 = "1_Mutant",
  group.by = "cl_genotype", min.pct = cfg$da_min_pct,
  logfc.threshold = cfg$da_logfc_threshold, test.use = cfg$da_test,
  latent.vars = cfg$da_latent_var
)

# Annotate closest genes
closest_mut <- ClosestFeature(obj, regions = rownames(da_peaks_mut))
rownames(closest_mut) <- closest_mut$query_region
da_peaks_mut$gene          <- closest_mut$gene_name
da_peaks_mut$gene_biotype  <- closest_mut$gene_biotype
da_peaks_mut$closest_region <- closest_mut$closest_region
da_peaks_mut$distance      <- closest_mut$distance
write.csv(da_peaks_mut, file.path(diffdir, "DA_peaks_Wnt7_vs_Lgr5_mutant.csv"))

cat(">> DA peaks: Wnt7+ vs Lgr5+ (wild-type)\n")
da_peaks_wt <- FindMarkers(
  obj, ident.1 = "3_Wild_type", ident.2 = "1_Wild_type",
  group.by = "cl_genotype", min.pct = cfg$da_min_pct,
  logfc.threshold = cfg$da_logfc_threshold, test.use = cfg$da_test,
  latent.vars = cfg$da_latent_var
)

closest_wt <- ClosestFeature(obj, regions = rownames(da_peaks_wt))
rownames(closest_wt) <- closest_wt$query_region
da_peaks_wt$gene          <- closest_wt$gene_name
da_peaks_wt$gene_biotype  <- closest_wt$gene_biotype
da_peaks_wt$closest_region <- closest_wt$closest_region
da_peaks_wt$distance      <- closest_wt$distance
write.csv(da_peaks_wt, file.path(diffdir, "DA_peaks_Wnt7_vs_Lgr5_wildtype.csv"))

# ==========================================================================
# 4) Wnt7b expression statistics (Fisher's exact per cluster)
# ==========================================================================
cat(">> Wnt7b expression statistics\n")
cl_pairs <- list(
  c("0_Wild_type", "0_Mutant"), c("1_Wild_type", "1_Mutant"),
  c("2_Wild_type", "2_Mutant"), c("3_Wild_type", "3_Mutant"),
  c("4_Wild_type", "4_Mutant"), c("5_Wild_type", "5_Mutant"),
  c("6_Wild_type", "6_Mutant")
)

# Use GetAssayData for Seurat v5 compatibility
rna_data <- GetAssayData(obj, assay = "RNA", layer = "data")

expr_stats <- data.frame()
for (pair in cl_pairs) {
  wt_cells  <- colnames(subset(obj, subset = cl_genotype == pair[1]))
  mut_cells <- colnames(subset(obj, subset = cl_genotype == pair[2]))

  wt_expr  <- rna_data["Wnt7b", wt_cells]
  mut_expr <- rna_data["Wnt7b", mut_cells]

  # Wilcoxon test
  wil <- wilcox.test(as.numeric(wt_expr), as.numeric(mut_expr))

  # Fisher's exact: proportion of Wnt7b+ cells
  wt_pos  <- sum(wt_expr != 0);  wt_neg  <- length(wt_expr) - wt_pos
  mut_pos <- sum(mut_expr != 0); mut_neg <- length(mut_expr) - mut_pos
  fisher_res <- fisher.test(matrix(c(wt_pos, mut_pos, wt_neg, mut_neg), nrow = 2),
                            alternative = "less")

  expr_stats <- rbind(expr_stats, data.frame(
    cluster_wt = pair[1], cluster_mut = pair[2],
    wt_n = length(wt_expr), mut_n = length(mut_expr),
    wt_pct_positive = wt_pos / length(wt_expr) * 100,
    mut_pct_positive = mut_pos / length(mut_expr) * 100,
    wilcox_pvalue = wil$p.value,
    fisher_pvalue = fisher_res$p.value,
    fisher_odds_ratio = fisher_res$estimate
  ))
}
write.csv(expr_stats, file.path(diffdir, "Wnt7b_expression_stats.csv"), row.names = FALSE)

# ==========================================================================
# 5) Wnt7b accessibility (peak region t-tests)
# ==========================================================================
cat(">> Wnt7b accessibility statistics\n")
acc_stats <- data.frame()
for (region_str in cfg$wnt7b_peak_regions) {
  region_counts <- CountsInRegion(
    object  = obj,
    assay   = "macs2",
    regions = StringToGRanges(region_str)
  )
  for (pair in cl_pairs) {
    wt_cells  <- colnames(subset(obj, subset = cl_genotype == pair[1]))
    mut_cells <- colnames(subset(obj, subset = cl_genotype == pair[2]))
    t_res <- t.test(region_counts[mut_cells], region_counts[wt_cells])

    acc_stats <- rbind(acc_stats, data.frame(
      region = region_str,
      cluster_wt = pair[1], cluster_mut = pair[2],
      mean_wt = t_res$estimate[2], mean_mut = t_res$estimate[1],
      t_pvalue = t_res$p.value
    ))
  }
}
write.csv(acc_stats, file.path(diffdir, "Wnt7b_accessibility_stats.csv"), row.names = FALSE)

cat(">> Step 06 complete.\n")