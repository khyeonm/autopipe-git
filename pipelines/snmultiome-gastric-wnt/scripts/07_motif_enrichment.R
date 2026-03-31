#!/usr/bin/env Rscript
# ===========================================================================
#  07_motif_enrichment.R — DA peaks + motif enrichment (Wnt7+ vs Lgr5+)
#  Replicates: FindMarkers on ATAC (LR test) → ClosestFeature annotation
#              → AddMotifs (JASPAR2020) → FindMotifs → enrichment bar plots
#
#  Performs the analysis separately for Mutant and Wild-type backgrounds,
#  exactly as in the original figure code (Fig 3I, S3C).
#
#  Also includes Wnt7b region ATAC enrichment testing (Fig 3H statistics)
#  and Wnt7b expression Wilcoxon tests per cluster.
# ===========================================================================

suppressPackageStartupMessages({
    library(yaml)
    library(future)
    library(Seurat)
    library(Signac)
    library(ggplot2)
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(JASPAR2020)
    library(TFBSTools)
    library(motifmatchr)
    library(chromVAR)
})

# ── Load configuration ────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
config_path <- args[which(args == "--config") + 1]
cfg <- read_yaml(config_path)

set.seed(cfg$seed)
# Use single worker for DA peak analysis (as in original code: plan("multicore", workers=1))
plan("multicore", workers = 1)
options(future.globals.maxSize = cfg$future_ram_gb * 1024^3)
output_dir <- cfg$output_dir
linkage_genome <- BSgenome.Mmusculus.UCSC.mm10

cat("=== Step 7: Motif enrichment analysis ===\n")

# ── Load integrated object ────────────────────────────────────────────────
obj <- readRDS(file.path(output_dir, "rds", "integrated_clustered.rds"))

# ══════════════════════════════════════════════════════════════════════════
#  Part A: Differential accessibility (Wnt7+ vs Lgr5+)
# ══════════════════════════════════════════════════════════════════════════
DefaultAssay(obj) <- "ATAC"

# ── DA peaks: Wnt7+ (cluster 3) vs Lgr5+ (cluster 1) in Mutant ──────────
cat("Finding DA peaks: Wnt7+ vs Lgr5+ in Mutant...\n")
da_peaks_mutant <- FindMarkers(
    obj,
    ident.1    = "3_Mutant",
    ident.2    = "1_Mutant",
    group.by   = "cl_genotype",
    min.pct    = cfg$da_min_pct,
    logfc.threshold = cfg$da_logfc_threshold,
    test.use   = cfg$da_test,
    latent.vars = "atac_peak_region_fragments"
)

# Annotate DA peaks with closest gene
closest_mut <- ClosestFeature(obj, regions = rownames(da_peaks_mutant))
rownames(closest_mut) <- closest_mut$query_region
da_peaks_mutant$gene         <- closest_mut$gene_name
da_peaks_mutant$gene_biotype <- closest_mut$gene_biotype
da_peaks_mutant$closest_region <- closest_mut$closest_region
da_peaks_mutant$distance     <- closest_mut$distance

write.csv(da_peaks_mutant,
          file.path(output_dir, "tables", "DA_peaks_Wnt_vs_Lgr5_mutant.csv"))
cat(sprintf("  Mutant DA peaks: %d total, %d significant (padj<0.05)\n",
            nrow(da_peaks_mutant),
            sum(da_peaks_mutant$p_val_adj < 0.05, na.rm = TRUE)))

# ── DA peaks: Wnt7+ vs Lgr5+ in Wild-type ────────────────────────────────
cat("Finding DA peaks: Wnt7+ vs Lgr5+ in Wild-type...\n")
da_peaks_wildtype <- FindMarkers(
    obj,
    ident.1    = "3_Wild_type",
    ident.2    = "1_Wild_type",
    group.by   = "cl_genotype",
    min.pct    = cfg$da_min_pct,
    logfc.threshold = cfg$da_logfc_threshold,
    test.use   = cfg$da_test,
    latent.vars = "atac_peak_region_fragments"
)

closest_wt <- ClosestFeature(obj, regions = rownames(da_peaks_wildtype))
rownames(closest_wt) <- closest_wt$query_region
da_peaks_wildtype$gene         <- closest_wt$gene_name
da_peaks_wildtype$gene_biotype <- closest_wt$gene_biotype
da_peaks_wildtype$closest_region <- closest_wt$closest_region
da_peaks_wildtype$distance     <- closest_wt$distance

write.csv(da_peaks_wildtype,
          file.path(output_dir, "tables", "DA_peaks_Wnt_vs_Lgr5_wildtype.csv"))
cat(sprintf("  Wild-type DA peaks: %d total, %d significant (padj<0.05)\n",
            nrow(da_peaks_wildtype),
            sum(da_peaks_wildtype$p_val_adj < 0.05, na.rm = TRUE)))

# ══════════════════════════════════════════════════════════════════════════
#  Part B: Motif enrichment in DA peaks (JASPAR2020)
# ══════════════════════════════════════════════════════════════════════════
cat("Adding JASPAR2020 motifs to object...\n")
pfm <- getMatrixSet(
    x    = JASPAR2020,
    opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)
obj <- AddMotifs(object = obj, genome = linkage_genome, pfm = pfm)

# ── Helper: run FindMotifs and generate barplot ───────────────────────────
run_motif_enrichment <- function(da_peaks, comparison_label, pdf_path, csv_path) {
    sig_peaks <- rownames(da_peaks[da_peaks$p_val_adj < cfg$da_pval_cutoff, ])
    cat(sprintf("  %s: %d peaks for motif enrichment (padj < %g)\n",
                comparison_label, length(sig_peaks), cfg$da_pval_cutoff))

    if (length(sig_peaks) == 0) {
        cat("  No significant peaks — skipping motif enrichment.\n")
        # Create empty outputs
        write.csv(data.frame(), csv_path, row.names = FALSE)
        pdf(pdf_path, width = 6, height = 6); plot.new(); dev.off()
        return(NULL)
    }

    enriched <- FindMotifs(object = obj, features = sig_peaks)
    enriched$logP <- -log10(enriched$p.adjust)

    # Save full enrichment table
    write.csv(enriched, csv_path, row.names = FALSE)

    # Significant motifs
    sig_motifs <- enriched[enriched$p.adjust < cfg$motif_padj_cutoff, ]
    cat(sprintf("  Significant motifs (padj < %g): %d\n",
                cfg$motif_padj_cutoff, nrow(sig_motifs)))

    # Plot top 15 motifs (as in original code)
    data_plot <- head(enriched, 15)
    y_orders  <- rev(data_plot$motif.name)

    p <- ggplot(data_plot, aes(x = logP, y = motif.name,
                                fill = percent.observed)) +
        geom_bar(stat = "identity", width = 0.2) +
        scale_fill_gradientn(
            colours = c("royalblue", "rosybrown2", "red"),
            limits  = c(0, 60)
        ) +
        scale_y_discrete(limits = y_orders) +
        xlim(0, 30) +
        ggtitle(comparison_label) +
        xlab("-log10(p.adjust)") +
        ylab("TF motif") +
        theme_linedraw() +
        theme(
            plot.title   = element_text(hjust = 0.5, color = "black",
                                        size = 15, face = "bold"),
            axis.title.x = element_text(color = "black", size = 12),
            axis.title.y = element_text(color = "black", size = 12)
        )

    pdf(pdf_path, width = 6, height = 6)
    print(p)
    dev.off()

    return(enriched)
}

# ── Mutant: Wnt7+ vs Lgr5+ motifs ────────────────────────────────────────
cat("Motif enrichment: Wnt7+ vs Lgr5+ (Mutant)...\n")
motifs_mut <- run_motif_enrichment(
    da_peaks_mutant,
    "Wnt7b+ vs Lgr5+\n(RZK)",
    file.path(output_dir, "plots", "fig3I_motifs_Wnt_vs_Lgr5_mutant.pdf"),
    file.path(output_dir, "tables", "enriched_motifs_mutant.csv")
)

# ── Wild-type: Wnt7+ vs Lgr5+ motifs ─────────────────────────────────────
cat("Motif enrichment: Wnt7+ vs Lgr5+ (Wild-type)...\n")
motifs_wt <- run_motif_enrichment(
    da_peaks_wildtype,
    "Wnt7b+ vs Lgr5+\n(RZ)",
    file.path(output_dir, "plots", "fig3I_motifs_Wnt_vs_Lgr5_wildtype.pdf"),
    file.path(output_dir, "tables", "enriched_motifs_wildtype.csv")
)

# ══════════════════════════════════════════════════════════════════════════
#  Part C: Wnt7b expression & ATAC region statistics (Fig 3H)
# ══════════════════════════════════════════════════════════════════════════
cat("\n--- Wnt7b expression Wilcoxon tests per cluster ---\n")
cluster_nums <- c("0", "1", "2", "3", "4", "5", "6")
DefaultAssay(obj) <- "RNA"

stats_results <- data.frame()
for (cl in cluster_nums) {
    wt_cells  <- colnames(subset(obj, subset = cl_genotype == paste0(cl, "_Wild_type")))
    mut_cells <- colnames(subset(obj, subset = cl_genotype == paste0(cl, "_Mutant")))

    if (length(wt_cells) > 1 && length(mut_cells) > 1) {
        wt_expr  <- obj[["RNA"]]@data["Wnt7b", wt_cells]
        mut_expr <- obj[["RNA"]]@data["Wnt7b", mut_cells]
        result   <- tryCatch(wilcox.test(mut_expr, wt_expr), error = function(e) NULL)
        if (!is.null(result)) {
            cat(sprintf("  Cluster %s: p = %g\n", cl, result$p.value))
            stats_results <- rbind(stats_results, data.frame(
                cluster = cl, p_value = result$p.value
            ))
        }
    }
}

# Wnt7b ATAC region enrichment tests
cat("\n--- Wnt7b ATAC region enrichment tests ---\n")
DefaultAssay(obj) <- "macs2"
for (region_str in cfg$wnt7b_regions) {
    cat(sprintf("\nRegion: %s\n", region_str))
    region_counts <- CountsInRegion(
        object  = obj,
        assay   = "macs2",
        regions = StringToGRanges(region_str)
    )
    for (cl in cluster_nums) {
        wt_id  <- paste0(cl, "_Wild_type")
        mut_id <- paste0(cl, "_Mutant")
        wt_cells  <- colnames(subset(obj, subset = cl_genotype == wt_id))
        mut_cells <- colnames(subset(obj, subset = cl_genotype == mut_id))
        if (length(wt_cells) > 1 && length(mut_cells) > 1) {
            result <- tryCatch(
                t.test(region_counts[mut_cells], region_counts[wt_cells]),
                error = function(e) NULL
            )
            if (!is.null(result)) {
                cat(sprintf("  Cluster %s (%s vs %s): p = %g\n",
                            cl, mut_id, wt_id, result$p.value))
            }
        }
    }
}

cat("\n=== Step 7 complete ===\n")