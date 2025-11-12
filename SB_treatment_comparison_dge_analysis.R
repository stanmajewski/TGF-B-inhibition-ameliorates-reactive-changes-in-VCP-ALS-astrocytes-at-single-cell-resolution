#!/usr/bin/env Rscript

# ============================================================================
# Treatment Comparison Differential Expression Analysis
# ============================================================================
# Description: Performs pseudobulk differential expression analysis comparing
#              treatment conditions in control using DESeq2, followed by GO enrichment
# 
# Author: Analysis Pipeline
# Last Updated: 2025-01-09
# ============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript treatment_comparison_dge_analysis.R <main_directory> <metadata_file> <seurat_object>
  
  Arguments:
    main_directory : Path to main analysis directory for outputs
    metadata_file  : Path to sample metadata CSV file (e.g., A_input_group_tab.csv)
    seurat_object  : Path to integrated Seurat object RDA file (e.g., B04/B04_seur_integr_labelled.rda)
  
  Example:
    Rscript treatment_comparison_dge_analysis.R \\
      /nemo/lab/patanir/home/users/majewss/seurat_int_all/ \\
      A_input_group_tab.csv \\
      B04/B04_seur_integr_labelled.rda
  ", call. = FALSE)
}

# Set environment and main directory
main_dir <- normalizePath(args[1], mustWork = TRUE)
metadata_file <- normalizePath(args[2], mustWork = TRUE)
seurat_file <- normalizePath(args[3], mustWork = TRUE)

# Record start time
start_time <- Sys.time()

setwd(main_dir)
options(bitmapType = "cairo")

# Log start time and parameters
cat("\n============================================================================\n")
cat("Treatment Comparison Differential Expression Analysis\n")
cat("============================================================================\n")
cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Main directory:", main_dir, "\n")
cat("Metadata file:", metadata_file, "\n")
cat("Seurat object:", seurat_file, "\n")
cat("R version:", R.version.string, "\n")
cat("============================================================================\n\n")

# Load packages
cat("Loading required packages...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(ggplot2)
  library(cowplot)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(dplyr)
  library(readr)
  library(Matrix)
  library(stringr)
  library(progeny)
  library(DESeq2)
  library(tidyr)
  library(tibble)
  library(ggrepel)
})
cat("All packages loaded successfully.\n\n")

# Create output directory structure
cat("Setting up output directories...\n")
out_dir <- paste0(main_dir, "/treatment_comparison_dge/")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Create subdirectories
subdirs <- c("deseq2_results", "volcano_plots", "go_enrichment")
for (subdir in subdirs) {
  dir_path <- paste0(out_dir, subdir, "/")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}
cat("Output directories created.\n\n")

# Load data
cat("Loading input data...\n")
gr_tab <- read_csv(metadata_file, show_col_types = FALSE)
cat("Loaded metadata for", nrow(gr_tab), "samples\n")

load(file = seurat_file)
cat("Loaded Seurat object with", ncol(seur), "cells\n\n")

# ============================================================================
# SUBSET AND PREPARE DATA
# ============================================================================

cat("============================================================================\n")
cat("Subsetting data for treatment comparison...\n")
cat("============================================================================\n")

# Subset the Seurat object based on the "group" column
seurat_subset <- subset(seur, subset = group %in% c("ctrl_u", "ctrl_SB"))
cat("Subsetted to", ncol(seurat_subset), "cells in ctrl_u and ctrl_SB groups\n\n")

# ============================================================================
# PSEUDOBULK AGGREGATION
# ============================================================================

cat("Performing pseudobulk aggregation...\n")
pseudobulk <- AggregateExpression(seurat_subset, group.by = "sample")
count_matrix <- pseudobulk$RNA
cat("Pseudobulk aggregation complete\n")
cat("Count matrix dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n\n")

# ============================================================================
# PREPARE SAMPLE METADATA
# ============================================================================

cat("Preparing sample metadata for DESeq2...\n")

# Automatically assign conditions based on sample names
gr_tab <- gr_tab %>%
  mutate(sample = gsub("_", "-", sample))  # Replace underscores with dashes

# Create sample info dataframe
sample_info <- gr_tab %>%
  select(sample, line, treatment, group) %>%
  filter(sample %in% colnames(count_matrix)) %>%
  column_to_rownames("sample")

# Ensure 'line' and 'group' are factors
sample_info$line <- as.factor(sample_info$line)
sample_info$group <- as.factor(sample_info$group)
sample_info$treatment <- as.factor(sample_info$treatment)

# Reorder to match count matrix columns
sample_info <- sample_info[colnames(count_matrix), , drop = FALSE]

cat("Sample info prepared for", nrow(sample_info), "samples\n")
cat("Treatment groups:", paste(unique(sample_info$treatment), collapse = ", "), "\n")
cat("Cell lines:", paste(unique(sample_info$line), collapse = ", "), "\n\n")

# ============================================================================
# DESEQ2 DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================

cat("============================================================================\n")
cat("Running DESeq2 analysis...\n")
cat("============================================================================\n")

# Create DESeq2 object
cat("Creating DESeq2 dataset...\n")
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix, 
  colData = sample_info, 
  design = ~ line + treatment
)

# Pre-filter low count genes
cat("Filtering low count genes...\n")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Retained", sum(keep), "genes after filtering (>= 10 total counts)\n\n")

# Set reference level for treatment
dds$treatment <- relevel(dds$treatment, ref = "untreated")
cat("Reference level set to 'untreated'\n\n")

# Run DESeq2
cat("Running DESeq2 normalization and statistical testing...\n")
dds <- DESeq(dds)
cat("DESeq2 analysis complete\n\n")

# Get results
cat("Extracting results for SB431542 vs untreated...\n")
results <- results(dds, contrast = c("treatment", "SB431542", "untreated"))

# Create results dataframe
results_df <- as.data.frame(results) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)

# Save full results
write.csv(results_df, paste0(out_dir, "deseq2_results/deseq2_full_results.csv"), row.names = FALSE)
cat("Full results saved\n\n")

# Filter for significant genes
significant_genes <- results_df %>%
  filter(!is.na(padj) & padj < 0.05) %>%
  arrange(padj)

cat("Significant genes (padj < 0.05):", nrow(significant_genes), "\n")
cat("  Upregulated:", sum(significant_genes$log2FoldChange > 0), "\n")
cat("  Downregulated:", sum(significant_genes$log2FoldChange < 0), "\n\n")

# Save significant genes
write.csv(significant_genes, paste0(out_dir, "deseq2_results/significant_genes.csv"), row.names = FALSE)

# Create summary table
summary_df <- data.frame(
  Comparison = "SB431542_vs_untreated",
  Total_Genes_Tested = nrow(results_df),
  Significant_Genes = nrow(significant_genes),
  Upregulated = sum(significant_genes$log2FoldChange > 0),
  Downregulated = sum(significant_genes$log2FoldChange < 0)
)
write.csv(summary_df, paste0(out_dir, "deseq2_results/summary.csv"), row.names = FALSE)

# ============================================================================
# VOLCANO PLOT VISUALIZATION
# ============================================================================

cat("============================================================================\n")
cat("Creating volcano plot...\n")
cat("============================================================================\n")

# Prepare data for volcano plot
volcano_df <- results_df %>%
  mutate(
    Significance = case_when(
      is.na(padj) ~ "Not Significant",
      padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

# Identify top genes to label
# Top 10 by fold change (up and down)
top_up_fc <- volcano_df %>%
  filter(Significance == "Upregulated") %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

top_down_fc <- volcano_df %>%
  filter(Significance == "Downregulated") %>%
  arrange(log2FoldChange) %>%
  head(10)

# Top 10 by p-value (up and down)
top_up_pval <- volcano_df %>%
  filter(Significance == "Upregulated") %>%
  arrange(padj) %>%
  head(10)

top_down_pval <- volcano_df %>%
  filter(Significance == "Downregulated") %>%
  arrange(padj) %>%
  head(10)

# Combine and create label column
labeled_genes <- unique(c(top_up_fc$gene_id, top_down_fc$gene_id, 
                          top_up_pval$gene_id, top_down_pval$gene_id))

volcano_df <- volcano_df %>%
  mutate(Label = ifelse(gene_id %in% labeled_genes, gene_id, NA))

# Define FDR threshold
log_fdr_threshold <- -log10(0.05)

# Create volcano plot
p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  # Scatter plot points
  geom_point(alpha = 0.7, size = 2) +
  
  # Label selected genes with black leader lines
  geom_text_repel(
    aes(label = Label),
    size = 4,
    color = "black",
    fontface = "plain",
    box.padding = 1.2,
    point.padding = 0.6,
    nudge_y = 8,
    segment.color = "black",
    segment.size = 0.5,
    min.segment.length = 0.1,
    max.overlaps = Inf,
    na.rm = TRUE
  ) +
  
  # Horizontal FDR significance threshold (0.05)
  geom_hline(yintercept = log_fdr_threshold, linetype = "dashed", color = "black", size = 1) +
  
  # Custom colors for categories
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  
  # Limit Y-axis to 70
  ylim(0, 70) +
  
  # Title and labels
  labs(
    title = "SB431542 vs Untreated",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  
  # Theme improvements + Remove Legend
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  
  # Add directional arrows and labels
  annotate("text", x = -5, y = 67, label = "Down in SB431542", 
           color = "black", size = 5, fontface = "bold") +
  annotate("text", x = 5, y = 67, label = "Up in SB431542", 
           color = "black", size = 5, fontface = "bold") +
  annotate("segment", x = -3, xend = -7, y = 65, yend = 65, 
           arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 1) +
  annotate("segment", x = 3, xend = 7, y = 65, yend = 65, 
           arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 1)

# Save volcano plot
ggsave(
  filename = paste0(out_dir, "volcano_plots/volcano_plot_labeled.pdf"),
  plot = p_volcano,
  width = 12,
  height = 10,
  dpi = 300
)

cat("Volcano plot saved\n\n")

# ============================================================================
# GO ENRICHMENT ANALYSIS
# ============================================================================

cat("============================================================================\n")
cat("Running GO enrichment analysis...\n")
cat("============================================================================\n")

# Ensure significant_genes is properly formatted
sig_genes <- significant_genes %>% filter(!is.na(log2FoldChange))

# Separate upregulated and downregulated genes
upregulated_genes <- sig_genes %>% filter(log2FoldChange > 0) %>% pull(gene_id)
downregulated_genes <- sig_genes %>% filter(log2FoldChange < 0) %>% pull(gene_id)

cat("Upregulated genes for GO:", length(upregulated_genes), "\n")
cat("Downregulated genes for GO:", length(downregulated_genes), "\n\n")

# Convert gene symbols to Entrez IDs
cat("Converting gene symbols to Entrez IDs...\n")
up_entrez_ids <- na.omit(mapIds(
  org.Hs.eg.db,
  keys = upregulated_genes,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
))

down_entrez_ids <- na.omit(mapIds(
  org.Hs.eg.db,
  keys = downregulated_genes,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
))

cat("Converted", length(up_entrez_ids), "upregulated genes to Entrez IDs\n")
cat("Converted", length(down_entrez_ids), "downregulated genes to Entrez IDs\n\n")

# Run GO enrichment analysis for upregulated genes
cat("Running GO enrichment for upregulated genes...\n")
go_up <- enrichGO(
  gene = up_entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Run GO enrichment analysis for downregulated genes
cat("Running GO enrichment for downregulated genes...\n")
go_down <- enrichGO(
  gene = down_entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Convert results to dataframes
go_up_results <- as.data.frame(go_up)
go_down_results <- as.data.frame(go_down)

cat("\nGO enrichment complete\n")
cat("Upregulated GO terms:", nrow(go_up_results), "\n")
cat("Downregulated GO terms:", nrow(go_down_results), "\n\n")

# Save results to CSV
write.csv(go_up_results, paste0(out_dir, "go_enrichment/GO_upregulated_results.csv"), row.names = FALSE)
write.csv(go_down_results, paste0(out_dir, "go_enrichment/GO_downregulated_results.csv"), row.names = FALSE)

# Create individual visualizations for upregulated genes
if (nrow(go_up_results) > 0) {
  cat("Creating GO visualizations for upregulated genes...\n")
  
  # Bar plot
  pdf(paste0(out_dir, "go_enrichment/GO_barplot_upregulated.pdf"), width = 12, height = 10)
  print(barplot(go_up, showCategory = 15, title = "Top 15 GO Biological Processes (Upregulated)"))
  dev.off()
  
  # Dot plot
  pdf(paste0(out_dir, "go_enrichment/GO_dotplot_upregulated.pdf"), width = 12, height = 10)
  print(dotplot(go_up, showCategory = 15, title = "GO Enrichment Dotplot (Upregulated)"))
  dev.off()
}

# Create individual visualizations for downregulated genes
if (nrow(go_down_results) > 0) {
  cat("Creating GO visualizations for downregulated genes...\n")
  
  # Bar plot
  pdf(paste0(out_dir, "go_enrichment/GO_barplot_downregulated.pdf"), width = 12, height = 10)
  print(barplot(go_down, showCategory = 15, title = "Top 15 GO Biological Processes (Downregulated)"))
  dev.off()
  
  # Dot plot
  pdf(paste0(out_dir, "go_enrichment/GO_dotplot_downregulated.pdf"), width = 12, height = 10)
  print(dotplot(go_down, showCategory = 15, title = "GO Enrichment Dotplot (Downregulated)"))
  dev.off()
}

cat("\n")

# ============================================================================
# COMBINED GO PLOT (MIRROR BAR PLOT)
# ============================================================================

cat("============================================================================\n")
cat("Creating combined GO mirror plot...\n")
cat("============================================================================\n")

# Add regulation direction to the GO results
go_up_results_labeled <- as.data.frame(go_up) %>% mutate(Regulation = "Upregulated")
go_down_results_labeled <- as.data.frame(go_down) %>% mutate(Regulation = "Downregulated")

# Merge results
go_combined <- bind_rows(go_up_results_labeled, go_down_results_labeled)

# Order by Adjusted P-Value
go_combined <- go_combined %>% arrange(p.adjust)

# Select top 20 most significant upregulated and downregulated GO terms
top_up <- go_combined %>%
  filter(Regulation == "Upregulated") %>%
  slice_min(p.adjust, n = 20)

top_down <- go_combined %>%
  filter(Regulation == "Downregulated") %>%
  slice_min(p.adjust, n = 20)

# Reverse order for downregulated terms so most significant is at the bottom
top_down <- top_down %>% arrange(desc(p.adjust))

# Combine top up and down terms
go_top_combined <- bind_rows(top_up, top_down)

# Reverse sign for Upregulated GO terms to make them face left
go_top_combined <- go_top_combined %>%
  mutate(AdjustedP = ifelse(Regulation == "Upregulated", -log10(p.adjust), log10(p.adjust)))

# Set factor levels for correct ordering
go_top_combined$Regulation <- factor(go_top_combined$Regulation, 
                                     levels = c("Downregulated", "Upregulated"))

# Define custom colors
color_palette <- c("Upregulated" = "#E60000", "Downregulated" = "#003366")

# Ensure correct ordering of GO terms
go_top_combined <- go_top_combined %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

# Create the mirror bar plot
p_go_mirror <- ggplot(go_top_combined, aes(y = Description, x = AdjustedP, fill = Regulation)) +
  geom_col(show.legend = TRUE) +
  scale_x_continuous(labels = function(x) abs(x)) +
  scale_fill_manual(values = color_palette) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  labs(
    title = "GO Biological Processes: SB431542 vs Untreated",
    x = "-Log10 Adjusted P-value",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save combined GO plot
ggsave(
  filename = paste0(out_dir, "go_enrichment/GO_combined_mirror_plot.pdf"),
  plot = p_go_mirror,
  width = 14,
  height = 12,
  dpi = 300
)

cat("Combined GO mirror plot saved\n\n")

# ============================================================================
# FINAL SUMMARY AND COMPLETION
# ============================================================================

cat("\n============================================================================\n")
cat("ANALYSIS SUMMARY\n")
cat("============================================================================\n")

cat("\nDESeq2 Results:\n")
cat("  Comparison: SB431542 vs untreated\n")
cat("  Total genes tested:", nrow(results_df), "\n")
cat("  Significant genes (padj < 0.05):", nrow(significant_genes), "\n")
cat("  Upregulated:", sum(significant_genes$log2FoldChange > 0), "\n")
cat("  Downregulated:", sum(significant_genes$log2FoldChange < 0), "\n")

cat("\nGO Enrichment Results:\n")
cat("  Upregulated GO terms:", nrow(go_up_results), "\n")
cat("  Downregulated GO terms:", nrow(go_down_results), "\n")

if (nrow(go_up_results) > 0) {
  cat("\n  Top 5 upregulated GO terms:\n")
  for (i in 1:min(5, nrow(go_up_results))) {
    cat(sprintf("    %d. %s (p.adjust = %.2e)\n",
                i,
                go_up_results$Description[i],
                go_up_results$p.adjust[i]))
  }
}

if (nrow(go_down_results) > 0) {
  cat("\n  Top 5 downregulated GO terms:\n")
  for (i in 1:min(5, nrow(go_down_results))) {
    cat(sprintf("    %d. %s (p.adjust = %.2e)\n",
                i,
                go_down_results$Description[i],
                go_down_results$p.adjust[i]))
  }
}

cat("\n============================================================================\n")
cat("OUTPUT FILES LOCATION\n")
cat("============================================================================\n")
cat("Main output directory:", out_dir, "\n")
cat("  - deseq2_results/    : DESeq2 results and summary tables\n")
cat("  - volcano_plots/     : Volcano plot with labeled genes\n")
cat("  - go_enrichment/     : GO enrichment results and visualizations\n")

cat("\n============================================================================\n")
cat("COMPLETION STATUS\n")
cat("============================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("Total elapsed time:", round(elapsed, 2), "minutes\n")

cat("\nAnalysis completed successfully!\n")
cat("============================================================================\n\n")

# Save session info
session_file <- paste0(out_dir, "session_info.txt")
writeLines(capture.output(sessionInfo()), session_file)
cat("Session information saved to:", session_file, "\n\n")
