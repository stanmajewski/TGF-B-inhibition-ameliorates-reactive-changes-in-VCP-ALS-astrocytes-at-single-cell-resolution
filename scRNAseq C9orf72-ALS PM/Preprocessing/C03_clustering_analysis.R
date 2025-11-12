#!/usr/bin/env Rscript

################################################################################
# Script: C03_clustering_analysis.R
# Description: Test clustering resolutions and characterize cell populations
#
# This script:
#   1. Tests multiple clustering resolutions (0.3 - 1.5)
#   2. Generates UMAP visualizations for each resolution
#   3. Creates marker gene dot plots for cluster characterization
#   4. Visualizes batch effects across samples
#   5. Creates cluster assignment template for manual annotation
#
# Inputs:
#   - results/02_integration/integrated_seurat.rda: Harmony-integrated dataset
#   - results/02_integration/filtered_metadata.csv: Sample metadata
#   - data/reference/cell_type_markers.csv: Cell type marker genes
#
# Outputs:
#   - results/03_clustering/clustered_seurat.rda: Dataset with cluster assignments
#   - results/03_clustering/umap_clustering_test.pdf: UMAP plots
#   - results/03_clustering/dotplot_cell_type_markers.pdf: Cell type marker expression
#   - results/03_clustering/dotplot_subtype_markers.pdf: Subtype marker expression
#   - results/03_clustering/umap_by_sample.pdf: Sample integration check
#   - results/03_clustering/umap_marker_expression.pdf: Individual marker UMAPs
#   - results/03_clustering/cluster_assignment_template.csv: Template for annotation
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(colorRamps)
})

# Initialize logging
cat("\n")
cat("========================================================================\n")
cat("C03: Clustering Analysis\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n\n")

# Set random seed for reproducibility
set.seed(1234)

# Define paths
project_dir <- getwd()
input_file <- "results/02_integration/integrated_seurat.rda"
metadata_file <- "results/02_integration/filtered_metadata.csv"
markers_file <- "data/reference/cell_type_markers.csv"
output_dir <- "results/03_clustering"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

################################################################################
# Define Functions
################################################################################

# Color palette function for categorical variables
generate_palette <- function(values) {
  n_values <- length(unique(values))
  
  if (n_values == 2) {
    colors <- c("grey", "blue")
  } else if (n_values == 3) {
    colors <- c("blue", "grey", "orange")
  } else if (n_values < 6) {
    colors <- matlab.like(6)[1:n_values]
  } else {
    colors <- matlab.like(n_values)
  }
  
  return(colors)
}

################################################################################
# Load Data
################################################################################

cat("\n>>> Loading integrated dataset...\n")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}
load(input_file)
cat("Loaded dataset with", ncol(filtered_seurat), "cells\n")

cat("\n>>> Loading sample metadata...\n")
sample_metadata <- read_csv(metadata_file, show_col_types = FALSE)

cat("\n>>> Loading marker gene lists...\n")
if (!file.exists(markers_file)) {
  warning("Marker file not found: ", markers_file)
  cat("Using default marker genes\n")
  
  # Default marker genes if file not available
  cell_type_markers <- c("RBFOX3", "SLC17A7", "GAD1", "GAD2", "GFAP", "AQP4", 
                          "AIF1", "MBP", "OLIG2", "PDGFRA", "PECAM1", "VWF")
  subtype_markers <- c("SATB2", "BCL11B", "FEZF2", "PCP4", "PVALB", "SST", 
                        "VIP", "SLC1A3", "GJA1")
} else {
  marker_data <- read_csv(markers_file, show_col_types = FALSE)
  cell_type_markers <- marker_data %>%
    filter(level %in% c("cell_types", "neuronal_lineage")) %>%
    pull(gene)
  subtype_markers <- marker_data %>%
    filter(level %in% c("NPC_patterning", "cortical_layers", "Interneuron_subtypes")) %>%
    pull(gene)
}

cat("Cell type markers:", length(cell_type_markers), "\n")
cat("Subtype markers:", length(subtype_markers), "\n")

################################################################################
# Test Multiple Clustering Resolutions
################################################################################

cat("\n>>> Testing clustering resolutions...\n")

# Define resolutions to test
test_resolutions <- c(0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5)
cat("Testing resolutions:", paste(test_resolutions, collapse = ", "), "\n")

# Store full dataset
full_seurat <- filtered_seurat

# Subsample for plotting if needed
if (ncol(filtered_seurat) > 100000) {
  cat("\nDataset is large (", ncol(filtered_seurat), "cells).\n")
  cat("Subsampling to 100,000 cells for visualization.\n")
  plot_cells <- sample(colnames(filtered_seurat), 100000)
  plot_seurat <- filtered_seurat[, plot_cells]
} else {
  plot_seurat <- filtered_seurat
}

# Initialize plot lists
umap_plots <- list()
dotplot_cell_types <- list()
dotplot_subtypes <- list()

# Generate UMAP plots by group and sample
cat("\n>>> Generating base UMAP visualizations...\n")

groups <- unique(sample_metadata$group)
samples <- unique(sample_metadata$sample)

umap_plots[["by_group"]] <- DimPlot(
  plot_seurat,
  group.by = "group",
  shuffle = TRUE,
  label = FALSE,
  reduction = "umap",
  pt.size = 0.01
) +
  scale_color_manual(limits = groups, values = generate_palette(groups)) +
  labs(title = "UMAP by Group")

umap_plots[["by_sample"]] <- DimPlot(
  plot_seurat,
  group.by = "sample",
  shuffle = TRUE,
  label = FALSE,
  reduction = "umap",
  pt.size = 0.01
) +
  scale_color_manual(limits = samples, values = generate_palette(samples)) +
  labs(title = "UMAP by Sample")

# Test each clustering resolution
cat("\n>>> Clustering at multiple resolutions...\n")

for (resolution in test_resolutions) {
  cat("  Testing resolution:", resolution, "\n")
  
  # Cluster full dataset
  full_seurat <- FindClusters(
    object = full_seurat,
    graph.name = "SCT_snn",
    algorithm = 1,
    resolution = resolution,
    verbose = FALSE
  )
  
  # Update plot dataset
  if (ncol(filtered_seurat) > 100000) {
    plot_seurat <- full_seurat[, plot_cells]
  } else {
    plot_seurat <- full_seurat
  }
  
  # UMAP plot for this resolution
  clusters <- unique(plot_seurat$seurat_clusters)
  res_name <- paste0("resolution_", resolution)
  
  umap_plots[[res_name]] <- DimPlot(
    plot_seurat,
    group.by = "seurat_clusters",
    shuffle = TRUE,
    label = TRUE,
    reduction = "umap",
    pt.size = 0.01
  ) +
    scale_color_manual(limits = clusters, values = generate_palette(clusters)) +
    NoLegend() +
    labs(title = paste("Clusters (resolution", resolution, ")"))
  
  # Dot plots for marker expression
  DefaultAssay(full_seurat) <- "SCT"
  
  # Cell type markers
  available_cell_markers <- intersect(cell_type_markers, rownames(full_seurat))
  if (length(available_cell_markers) > 0) {
    dotplot_cell_types[[res_name]] <- DotPlot(
      full_seurat,
      features = available_cell_markers,
      scale.by = "size"
    ) +
      RotatedAxis() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(title = paste("Cell Type Markers (resolution", resolution, ")"))
  }
  
  # Subtype markers
  available_subtype_markers <- intersect(subtype_markers, rownames(full_seurat))
  if (length(available_subtype_markers) > 0) {
    dotplot_subtypes[[res_name]] <- DotPlot(
      full_seurat,
      features = available_subtype_markers,
      scale.by = "size"
    ) +
      RotatedAxis() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(title = paste("Subtype Markers (resolution", resolution, ")"))
  }
}

# Restore full dataset
filtered_seurat <- full_seurat

# Save clustered dataset
cat("\n>>> Saving clustered dataset...\n")
clustered_file <- file.path(output_dir, "clustered_seurat.rda")
save(filtered_seurat, file = clustered_file)
cat("Saved to:", clustered_file, "\n")

################################################################################
# Save UMAP and Dot Plots
################################################################################

cat("\n>>> Saving visualization plots...\n")

# UMAP clustering plots
umap_file <- file.path(output_dir, "umap_clustering_test.pdf")
pdf(umap_file, width = 6, height = 5)
lapply(umap_plots, print)
dev.off()
cat("Saved UMAP plots to:", umap_file, "\n")

# Cell type marker dot plots
if (length(dotplot_cell_types) > 0) {
  dotplot_cell_file <- file.path(output_dir, "dotplot_cell_type_markers.pdf")
  pdf(dotplot_cell_file, width = 12, height = 10)
  lapply(dotplot_cell_types, print)
  dev.off()
  cat("Saved cell type marker plots to:", dotplot_cell_file, "\n")
}

# Subtype marker dot plots
if (length(dotplot_subtypes) > 0) {
  dotplot_subtype_file <- file.path(output_dir, "dotplot_subtype_markers.pdf")
  pdf(dotplot_subtype_file, width = 12, height = 10)
  lapply(dotplot_subtypes, print)
  dev.off()
  cat("Saved subtype marker plots to:", dotplot_subtype_file, "\n")
}

################################################################################
# Sample Integration Check
################################################################################

cat("\n>>> Creating sample integration visualization...\n")

# Subsample if needed
if (ncol(filtered_seurat) > 100000) {
  plot_seurat <- filtered_seurat[, sample(colnames(filtered_seurat), 100000)]
} else {
  plot_seurat <- filtered_seurat
}

# Extract UMAP coordinates
umap_data <- FetchData(plot_seurat, vars = c("umap_1", "umap_2", "group", "sample"))

# Create faceted UMAP by sample
sample_umap <- ggplot(umap_data, aes(x = umap_1, y = umap_2)) +
  geom_point(size = 0.05, alpha = 0.5) +
  facet_wrap(~ factor(sample, levels = samples), nrow = 2) +
  theme_bw() +
  labs(title = "UMAP by Sample (Integration Check)")

# Save sample integration plot
sample_file <- file.path(output_dir, "umap_by_sample.pdf")
pdf(sample_file, width = 10, height = 2.5)
print(sample_umap)
dev.off()
cat("Saved sample integration plot to:", sample_file, "\n")

################################################################################
# Individual Marker Expression on UMAP
################################################################################

cat("\n>>> Creating individual marker expression UMAPs...\n")

# Subsample for marker plots
if (ncol(filtered_seurat) > 10000) {
  plot_seurat <- filtered_seurat[, sample(colnames(filtered_seurat), 10000)]
} else {
  plot_seurat <- filtered_seurat
}

DefaultAssay(plot_seurat) <- "SCT"

# Combine all markers
all_markers <- unique(c(cell_type_markers, subtype_markers))
available_markers <- intersect(all_markers, rownames(plot_seurat))

if (length(available_markers) > 0) {
  marker_plots <- FeaturePlot(
    plot_seurat,
    features = available_markers,
    order = TRUE,
    reduction = "umap",
    pt.size = 0.01,
    ncol = 9
  )
  
  # Save marker expression plots
  marker_file <- file.path(output_dir, "umap_marker_expression.pdf")
  pdf(marker_file, width = 30, height = 15)
  print(marker_plots)
  dev.off()
  cat("Saved marker expression UMAPs to:", marker_file, "\n")
}

################################################################################
# Create Cluster Assignment Template
################################################################################

cat("\n>>> Creating cluster assignment template...\n")

# Use highest resolution for template
highest_res_col <- paste0("SCT_snn_res.", max(test_resolutions))
clusters <- unique(filtered_seurat@meta.data[[highest_res_col]])
clusters <- sort(as.numeric(as.character(clusters)))

# Create template
cluster_template <- tibble(
  cluster = clusters,
  cluster_name = paste0("cluster_", clusters),
  cell_type = "unassigned",
  cell_class = "unassigned"
)

# Save template
template_file <- file.path(output_dir, "cluster_assignment_template.csv")
write_csv(cluster_template, template_file)
cat("Saved cluster assignment template to:", template_file, "\n")
cat("\nPlease manually annotate clusters in this file based on marker expression.\n")
cat("Cluster labels will be used in subsequent analysis scripts.\n")

################################################################################
# Completion Summary
################################################################################

cat("\n")
cat("========================================================================\n")
cat("C03: Clustering Analysis Complete\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n")
cat("\nSummary:\n")
cat("  Total cells:", ncol(filtered_seurat), "\n")
cat("  Resolutions tested:", length(test_resolutions), "\n")
cat("  Clusters at max resolution:", length(clusters), "\n")
cat("  Cell type markers available:", length(intersect(cell_type_markers, rownames(filtered_seurat))), "\n")
cat("  Subtype markers available:", length(intersect(subtype_markers, rownames(filtered_seurat))), "\n")
cat("\nOutput files:\n")
cat("  -", clustered_file, "\n")
cat("  -", umap_file, "\n")
cat("  -", template_file, "\n")
cat("\nNext steps:\n")
cat("  1. Review clustering plots to select optimal resolution\n")
cat("  2. Annotate clusters in", basename(template_file), "\n")
cat("  3. Run C04_cluster_characterization.R with annotated clusters\n")
cat("========================================================================\n\n")
