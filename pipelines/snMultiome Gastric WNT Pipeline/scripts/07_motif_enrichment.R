# ==============================================================================
# 07_motif_enrichment.R — TF motif enrichment on DA peaks (JASPAR2020)
# ==============================================================================
# Identifies enriched TF motifs in differentially accessible peaks
# between Wnt7+ and Lgr5+ populations, separately for mutant and WT.
# Generates horizontal bar plots (Fig 3I style).
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(EnsDb.Mmusculus.v79)
  library(motifmatchr)
  library(JASPAR2020)
  library(TFBSTools)
  library(chromVAR)
  library(ggplot2)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
cfg  <- read_yaml(args[1])

set.seed(1234)
outdir  <- cfg$output_dir
motifdir <- file.path(outdir, "07_motifs")
dir.create(motifdir, recursive = TRUE, showWarnings = FALSE)

linkage_genome <- BSgenome.Mmusculus.UCSC.mm10

# ---- Load -------------------------------------------------------------------
cat(">> Loading integrated object\n")
obj <- readRDS(file.path(outdir, "04_integrated", "RZ_RZK_integrated.rds"))

da_mut <- read.csv(file.path(outdir, "06_differential", "DA_peaks_Wnt7_vs_Lgr5_mutant.csv"),
                   row.names = 1)
da_wt  <- read.csv(file.path(outdir, "06_differential", "DA_peaks_Wnt7_vs_Lgr5_wildtype.csv"),
                   row.names = 1)

# ---- Add motifs to object ---------------------------------------------------
cat(">> Adding JASPAR motifs\n")
DefaultAssay(obj) <- "ATAC"
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = cfg$jaspar_collection,
              tax_group  = cfg$jaspar_tax_group,
              all_versions = FALSE)
)
obj <- AddMotifs(object = obj, genome = linkage_genome, pfm = pfm)

# ---- Helper: motif enrichment + bar plot ------------------------------------
run_motif_enrichment <- function(obj, da_df, padj_cutoff, title_str,
                                 csv_out, pdf_out) {
  sig_peaks <- rownames(da_df[da_df$p_val_adj < padj_cutoff, ])
  if (length(sig_peaks) == 0) {
    cat("   No significant DA peaks for", title_str, "— skipping.\n")
    write.csv(data.frame(), csv_out)
    pdf(pdf_out, width = 6, height = 6); plot.new(); dev.off()
    return(NULL)
  }

  enriched <- FindMotifs(object = obj, features = sig_peaks)
  write.csv(enriched, csv_out, row.names = FALSE)

  # Plot top 15
  enriched$logP <- -log10(enriched$p.adjust)
  sig_motifs <- enriched[enriched$p.adjust < cfg$motif_padj_cutoff, ]
  if (nrow(sig_motifs) == 0) {
    cat("   No significant motifs for", title_str, "\n")
    pdf(pdf_out, width = 6, height = 6); plot.new(); dev.off()
    return(enriched)
  }

  plot_data <- head(enriched, 15)
  y_order   <- rev(plot_data$motif.name)

  p <- ggplot(plot_data, aes(x = logP, y = motif.name, fill = percent.observed)) +
    geom_bar(stat = "identity", width = 0.2) +
    scale_fill_gradientn(colours = c("royalblue", "rosybrown2", "red"),
                         limits = c(0, 60)) +
    scale_y_discrete(limits = y_order) +
    xlim(0, 30) +
    ggtitle(title_str) + xlab("-log10(P)") + ylab("TF motif") +
    theme_linedraw() +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 15, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )

  ggsave(pdf_out, plot = p, width = 6, height = 6, dpi = 300)
  return(enriched)
}

# ---- Mutant: Wnt7+ vs Lgr5+ ------------------------------------------------
cat(">> Motif enrichment: Wnt7+ vs Lgr5+ (Mutant)\n")
run_motif_enrichment(
  obj, da_mut, cfg$da_padj_cutoff,
  "Wnt7b+ vs Lgr5+\n(RZK)",
  file.path(motifdir, "motif_enrichment_Wnt_vs_Lgr5_mutant.csv"),
  file.path(motifdir, "fig3I_motif_barplot_mutant.pdf")
)

# ---- Wild-type: Wnt7+ vs Lgr5+ ---------------------------------------------
cat(">> Motif enrichment: Wnt7+ vs Lgr5+ (Wild-type)\n")
run_motif_enrichment(
  obj, da_wt, cfg$da_padj_cutoff,
  "Wnt7b+ vs Lgr5+\n(RZ)",
  file.path(motifdir, "motif_enrichment_Wnt_vs_Lgr5_wildtype.csv"),
  file.path(motifdir, "fig3I_motif_barplot_wildtype.pdf")
)

cat(">> Step 07 complete.\n")