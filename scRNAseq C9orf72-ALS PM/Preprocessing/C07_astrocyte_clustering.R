#!/usr/bin/env Rscript

################################################################################
# Script: C07_astrocyte_clustering.R
# Description: Multi-resolution clustering of astrocyte population
#
# This script:
#   1. Tests multiple clustering resolutions (0.1 - 0.6)
#   2. Generates UMAP visualizations for each resolution
#   3. Creates astrocyte subtype marker dot plots
#   4. Generates marker expression heatmaps by cluster
#   5. Uses clustree to visualize cluster relationships
#   6. Creates cluster assignment template for manual annotation
#
# Inputs:
#   - results/06_astrocyte_analysis/integrated_astrocytes.rda: Integrated astrocytes
#   - results/02_integration/filtered_metadata.csv: Sample metadata
#   - data/reference/cell_type_markers.csv: Marker gene lists (including astrocyte subtypes)
#
# Outputs:
#   - results/07_astrocyte_clustering/clustered_astrocytes.rda: Dataset with cluster assignments
#   - results/07_astrocyte_clustering/umap_clustering_test.pdf: UMAP plots
#   - results/07_astrocyte_clustering/dotplot_cell_type_markers.pdf: Cell type marker dot plots
#   - results/07_astrocyte_clustering/dotplot_astrocyte_subtype_markers.pdf: Subtype marker plots
#   - results/07_astrocyte_clustering/heatmap_subtype_markers.pdf: Marker heatmaps
#   - results/07_astrocyte_clustering/umap_by_sample.pdf: Sample integration check
#   - results/07_astrocyte_clustering/umap_marker_expression.pdf: Individual marker UMAPs
#   - results/07_astrocyte_clustering/astrocyte_cluster_assignment_template.csv: Annotation template
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(colorRamps)
  library(pheatmap)
  library(clustree)
})

# Initialize logging
cat("\n")
cat("========================================================================\n")
cat("C07: Astrocyte Clustering Analysis\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n\n")

# Set random seed for reproducibility
set.seed(1234)

# Define paths
project_dir <- getwd()
input_file <- "results/06_astrocyte_analysis/integrated_astrocytes.rda"
metadata_file <- "results/02_integration/filtered_metadata.csv"
markers_file <- "data/reference/cell_type_markers.csv"
output_dir <- "results/07_astrocyte_clustering"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

################################################################################
# Define Functions
################################################################################

# Color palette function
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

cat("\n>>> Loading integrated astrocyte dataset...\n")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}
load(input_file)
cat("Loaded dataset with", ncol(astrocytes), "astrocytes\n")

cat("\n>>> Loading sample metadata...\n")
sample_metadata <- read_csv(metadata_file, show_col_types = FALSE)

cat("\n>>> Loading marker gene lists...\n")
if (!file.exists(markers_file)) {
  warning("Marker file not found: ", markers_file)
  cat("Using default astrocyte marker genes\n")
  
  # Default astrocyte markers
  cell_type_markers <- c("GFAP", "AQP4", "ALDH1L1", "S100B", "GLUL", "SLC1A3")
  astrocyte_subtype_markers <- c("GJA1", "GFAP", "AQP4", "GLUL", "SLC1A2", 
                                  "SLC1A3", "ALDOC", "CLU", "MFGE8")
} else {
  marker_data <- read_csv(markers_file, show_col_types = FALSE)
  cell_type_markers <- marker_data %>%
    filter(level %in% c("cell_types", "neuronal_lineage")) %>%
    pull(gene)
  astrocyte_subtype_markers <- marker_data %>%
    filter(level == "Astrocyte_subtypes") %>%
    pull(gene)
}

cat("Cell type markers:", length(cell_type_markers), "\n")
cat("Astrocyte subtype markers:", length(astrocyte_subtype_markers), "\n")

################################################################################
# Test Multiple Clustering Resolutions
################################################################################

cat("\n>>> Testing clustering resolutions...\n")

# Define resolutions to test for astrocytes (lower range for subtypes)
test_resolutions <- c(0.1, 0.2, 0.3, 0.5, 0.6)
cat("Testing resolutions:", paste(test_resolutions, collapse = ", "), "\n")

# Store full dataset
full_astrocytes <- astrocytes

# Subsample for plotting if needed
if (ncol(astrocytes) > 100000) {
  cat("\nDataset is large (", ncol(astrocytes), "cells).\n")
  cat("Subsampling to 100,000 cells for visualization.\n")
  plot_cells <- sample(colnames(astrocytes), 100000)
  plot_astrocytes <- astrocytes[, plot_cells]
} else {
  plot_astrocytes <- astrocytes
}

# Initialize plot lists
umap_plots <- list()
dotplot_cell_types <- list()
dotplot_subtypes <- list()
marker_heatmap_list <- list()

# Extract expression matrix for heatmap plotting
DefaultAssay(full_astrocytes) <- "SCT"
expression_matrix <- full_astrocytes@assays$SCT$data
available_subtype_markers <- intersect(astrocyte_subtype_markers, rownames(expression_matrix))
marker_expression_matrix <- expression_matrix[available_subtype_markers, ]

# Generate base UMAP plots
cat("\n>>> Generating base UMAP visualizations...\n")

groups <- unique(sample_metadata$group)
samples <- unique(sample_metadata$sample)

umap_plots[["by_group"]] <- DimPlot(
  plot_astrocytes,
  group.by = "group",
  shuffle = TRUE,
  label = FALSE,
  reduction = "umap",
  pt.size = 0.01
) +
  scale_color_manual(limits = groups, values = generate_palette(groups)) +
  labs(title = "Astrocytes by Group")

umap_plots[["by_sample"]] <- DimPlot(
  plot_astrocytes,
  group.by = "sample",
  shuffle = TRUE,
  label = FALSE,
  reduction = "umap",
  pt.size = 0.01
) +
  scale_color_manual(limits = samples, values = generate_palette(samples)) +
  labs(title = "Astrocytes by Sample")

# Test each clustering resolution
cat("\n>>> Clustering at multiple resolutions...\n")

for (resolution in test_resolutions) {
  cat("  Testing resolution:", resolution, "\n")
  
  # Cluster full dataset
  full_astrocytes <- FindClusters(
    object = full_astrocytes,
    graph.name = "SCT_snn",
    algorithm = 1,
    resolution = resolution,
    verbose = FALSE
  )
  
  # Update plot dataset
  if (ncol(astrocytes) > 100000) {
    plot_astrocytes <- full_astrocytes[, plot_cells]
  } else {
    plot_astrocytes <- full_astrocytes
  }
  
  # UMAP plot for this resolution
  clusters <- unique(plot_astrocytes$seurat_clusters)
  res_name <- paste0("resolution_", resolution)
  
  umap_plots[[res_name]] <- DimPlot(
    plot_astrocytes,
    group.by = "seurat_clusters",
    shuffle = TRUE,
    label = TRUE,
    reduction = "umap",
    pt.size = 0.01
  ) +
    scale_color_manual(limits = clusters, values = generate_palette(clusters)) +
    NoLegend() +
    labs(title = paste("Astrocyte Clusters (resolution", resolution, ")"))
  
  # Dot plots for marker expression
  DefaultAssay(full_astrocytes) <- "SCT"
  
  # Cell type markers
  available_cell_markers <- intersect(cell_type_markers, rownames(full_astrocytes))
  if (length(available_cell_markers) > 0) {
    dotplot_cell_types[[res_name]] <- DotPlot(
      full_astrocytes,
      features = available_cell_markers,
      scale.by = "size"
    ) +
      RotatedAxis() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(title = paste("Cell Type Markers (resolution", resolution, ")"))
  }
  
  # Astrocyte subtype markers
  if (length(available_subtype_markers) > 0) {
    dotplot_subtypes[[res_name]] <- DotPlot(
      full_astrocytes,
      features = available_subtype_markers,
      scale.by = "size"
    ) +
      RotatedAxis() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(title = paste("Astrocyte Subtype Markers (resolution", resolution, ")"))
  }
  
  # Create marker by cluster matrix for heatmap
  cluster_matrix <- matrix(
    nrow = length(clusters),
    ncol = length(available_subtype_markers),
    dimnames = list(paste0("cluster_", clusters), available_subtype_markers)
  )
  
  for (cluster_id in clusters) {
    cluster_cells <- colnames(full_astrocytes)[full_astrocytes$seurat_clusters == cluster_id]
    cluster_expression <- marker_expression_matrix[, colnames(marker_expression_matrix) %in% cluster_cells]
    cluster_means <- apply(as.matrix(cluster_expression), 1, mean)
    cluster_matrix[paste0("cluster_", cluster_id), ] <- cluster_means
  }
  
  # Calculate Z-scores across clusters
  zscore_matrix <- apply(cluster_matrix, 2, scale)
  rownames(zscore_matrix) <- rownames(cluster_matrix)
  
  marker_heatmap_list[[res_name]] <- zscore_matrix
}

# Visualize cluster stability
cat("\n>>> Creating clustree visualization...\n")
clustree_plot <- clustree(full_astrocytes, prefix = "SCT_snn_res.")

# Restore full dataset
astrocytes <- full_astrocytes

# Save clustered dataset
cat("\n>>> Saving clustered dataset...\n")
clustered_file <- file.path(output_dir, "clustered_astrocytes.rda")
save(astrocytes, file = clustered_file)
cat("Saved to:", clustered_file, "\n")

################################################################################
# Save Visualizations
################################################################################

cat("\n>>> Saving visualization plots...\n")

# UMAP clustering plots
umap_file <- file.path(output_dir, "umap_clustering_test.pdf")
pdf(umap_file, width = 6, height = 5)
lapply(umap_plots, print)
print(clustree_plot)
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

# Astrocyte subtype marker dot plots
if (length(dotplot_subtypes) > 0) {
  dotplot_subtype_file <- file.path(output_dir, "dotplot_astrocyte_subtype_markers.pdf")
  pdf(dotplot_subtype_file, width = 36, height = 10)
  lapply(dotplot_subtypes, print)
  dev.off()
  cat("Saved astrocyte subtype marker plots to:", dotplot_subtype_file, "\n")
}

# Marker expression heatmaps
if (length(marker_heatmap_list) > 0) {
  heatmap_file <- file.path(output_dir, "heatmap_subtype_markers.pdf")
  pdf(heatmap_file, width = 36, height = 10)
  
  for (resol_name in names(marker_heatmap_list)) {
    zscore_matrix <- marker_heatmap_list[[resol_name]]
    plot_limits <- c(-max(abs(zscore_matrix)), max(abs(zscore_matrix)))
    
    pheatmap(
      zscore_matrix,
      show_rownames = TRUE,
      show_colnames = TRUE,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      clustering_distance_rows = "euclidean",
      clustering_method = "ward.D2",
      treeheight_row = 10,
      color = colorRampPalette(c("magenta", "black", "yellow"))(250),
      breaks = seq(plot_limits[1], plot_limits[2], length.out = 251),
      border_color = NA,
      fontsize = 10,
      cellwidth = 10,
      cellheight = 10,
      main = paste0(resol_name, " - Marker Z-score by Cluster")
    )
  }
  
  dev.off()
  cat("Saved marker heatmaps to:", heatmap_file, "\n")
}

################################################################################
# Sample Integration Check
################################################################################

cat("\n>>> Creating sample integration visualization...\n")

# Subsample if needed
if (ncol(astrocytes) > 100000) {
  plot_astrocytes <- astrocytes[, sample(colnames(astrocytes), 100000)]
} else {
  plot_astrocytes <- astrocytes
}

# Extract UMAP coordinates
umap_data <- FetchData(plot_astrocytes, vars = c("umap_1", "umap_2", "group", "sample"))

# Create faceted UMAP by sample
sample_umap <- ggplot(umap_data, aes(x = umap_1, y = umap_2)) +
  geom_point(size = 0.05, alpha = 0.5) +
  facet_wrap(~ factor(sample, levels = samples), nrow = 2) +
  theme_bw() +
  labs(title = "Astrocytes by Sample (Integration Check)")

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
if (ncol(astrocytes) > 10000) {
  plot_astrocytes <- astrocytes[, sample(colnames(astrocytes), 10000)]
} else {
  plot_astrocytes <- astrocytes
}

DefaultAssay(plot_astrocytes) <- "SCT"

# Combine all markers
all_markers <- unique(c(cell_type_markers, astrocyte_subtype_markers))
available_markers <- intersect(all_markers, rownames(plot_astrocytes))

if (length(available_markers) > 0) {
  marker_plots <- FeaturePlot(
    plot_astrocytes,
    features = available_markers,
    order = TRUE,
    reduction = "umap",
    pt.size = 0.01,
    ncol = 9
  )
  
  # Save marker expression plots
  marker_file <- file.path(output_dir, "umap_marker_expression.pdf")
  pdf(marker_file, width = 30, height = 45)
  print(marker_plots)
  dev.off()
  cat("Saved marker expression UMAPs to:", marker_file, "\n")
}

################################################################################
# Create Cluster Assignment Template
################################################################################

cat("\n>>> Creating astrocyte cluster assignment template...\n")

# Use highest resolution for template
highest_res_col <- paste0("SCT_snn_res.", max(test_resolutions))
clusters <- unique(astrocytes@meta.data[[highest_res_col]])
clusters <- sort(as.numeric(as.character(clusters)))

# Create template
cluster_template <- tibble(
  cluster = clusters,
  cluster_name = paste0("astrocyte_cluster_", clusters),
  astrocyte_subtype = "unassigned",
  astrocyte_state = "unassigned"
)

# Save template
template_file <- file.path(output_dir, "astrocyte_cluster_assignment_template.csv")
write_csv(cluster_template, template_file)
cat("Saved cluster assignment template to:", template_file, "\n")
cat("\nPlease manually annotate astrocyte clusters based on marker expression.\n")
cat("Consider astrocyte subtypes and states (e.g., homeostatic, reactive).\n")

################################################################################
# Completion Summary
################################################################################

cat("\n")
cat("========================================================================\n")
cat("C07: Astrocyte Clustering Analysis Complete\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n")
cat("\nSummary:\n")
cat("  Total astrocytes:", ncol(astrocytes), "\n")
cat("  Resolutions tested:", length(test_resolutions), "\n")
cat("  Clusters at max resolution:", length(clusters), "\n")
cat("  Astrocyte subtype markers available:", length(intersect(astrocyte_subtype_markers, rownames(astrocytes))), "\n")
cat("\nOutput files:\n")
cat("  -", clustered_file, "\n")
cat("  -", template_file, "\n")
cat("\nNext steps:\n")
cat("  1. Review clustering plots to select optimal resolution\n")
cat("  2. Annotate clusters in", basename(template_file), "\n")
cat("  3. Run C08_astrocyte_characterization.R with annotated clusters\n")
cat("========================================================================\n\n")
