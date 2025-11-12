#!/usr/bin/env Rscript

################################################################################
# Script: C04_cluster_characterization.R
# Description: Characterize clusters and perform differential abundance analysis
#
# This script:
#   1. Applies user-annotated cluster labels to dataset
#   2. Generates labeled UMAP visualizations
#   3. Creates marker expression dot plots with cluster names
#   4. Quantifies cell abundance across samples and clusters
#   5. Performs statistical tests for differential abundance
#   6. Runs sccomp analysis for compositional differences
#
# Inputs:
#   - results/03_clustering/clustered_seurat.rda: Clustered dataset
#   - results/02_integration/filtered_metadata.csv: Sample metadata
#   - results/03_clustering/cluster_assignment_template.csv: Annotated clusters
#   - data/reference/cell_type_markers.csv: Marker gene lists
#
# Outputs:
#   - results/04_characterization/labeled_seurat.rda: Dataset with cluster labels
#   - results/04_characterization/umap_labeled_clusters.pdf: Labeled UMAPs
#   - results/04_characterization/dotplot_labeled_clusters.pdf: Marker dot plots
#   - results/04_characterization/umap_by_sample_labeled.pdf: Sample comparison
#   - results/04_characterization/umap_marker_expression_labeled.pdf: Marker UMAPs
#   - results/04_characterization/abundance_by_sample_cluster.csv: Cell counts
#   - results/04_characterization/abundance_barplot.pdf: Abundance visualization
#   - results/04_characterization/abundance_heatmap.pdf: Abundance heatmap
#   - results/04_characterization/abundance_ttest.csv: T-test results
#   - results/04_characterization/sccomp_results.csv: sccomp analysis results
#   - results/04_characterization/sccomp_plots.pdf: sccomp visualizations
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(colorRamps)
  library(viridis)
  library(pheatmap)
  library(sccomp)
})

# Initialize logging
cat("\n")
cat("========================================================================\n")
cat("C04: Cluster Characterization\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n\n")

# Set random seed for reproducibility
set.seed(1234)

# Define paths
project_dir <- getwd()
input_file <- "results/03_clustering/clustered_seurat.rda"
metadata_file <- "results/02_integration/filtered_metadata.csv"
cluster_file <- "results/03_clustering/cluster_assignment_template.csv"
markers_file <- "data/reference/cell_type_markers.csv"
output_dir <- "results/04_characterization"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

################################################################################
# Define Functions
################################################################################

# Standard color palette
generate_palette <- function(values) {
  n_values <- length(unique(values))
  
  if (n_values == 2) {
    colors <- c("grey20", "dodgerblue")
  } else if (n_values == 3) {
    colors <- c("dodgerblue", "grey20", "orange")
  } else if (n_values < 6) {
    colors <- matlab.like(6)[1:n_values]
  } else {
    colors <- matlab.like(n_values)
  }
  
  return(colors)
}

# Distinct color palette (shuffled for better contrast)
generate_distinct_palette <- function(values) {
  n_values <- length(unique(values))
  
  if (n_values < 6) {
    colors <- matlab.like(6)[1:n_values]
  } else {
    colors <- matlab.like(n_values)
    set.seed(12)
    colors <- sample(colors)
  }
  
  return(colors)
}

# Heatmap function with viridis colors
create_abundance_heatmap <- function(matrix_data, title = "", 
                                     cluster_rows = FALSE, 
                                     cluster_cols = FALSE) {
  pheatmap(
    matrix_data,
    show_rownames = TRUE,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_colnames = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method = "ward.D2",
    treeheight_row = 50,
    color = viridis_pal(option = "magma")(250),
    breaks = seq(0, max(matrix_data), length.out = 251),
    border_color = NA,
    fontsize = 10,
    cellwidth = 10,
    cellheight = 10,
    main = title
  )
}

################################################################################
# Load Data
################################################################################

cat("\n>>> Loading clustered dataset...\n")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}
load(input_file)
cat("Loaded dataset with", ncol(filtered_seurat), "cells\n")

cat("\n>>> Loading sample metadata...\n")
sample_metadata <- read_csv(metadata_file, show_col_types = FALSE)

cat("\n>>> Loading cluster annotations...\n")
if (!file.exists(cluster_file)) {
  stop("Cluster annotation file not found: ", cluster_file,
       "\nPlease annotate clusters in cluster_assignment_template.csv first")
}
cluster_annotations <- read_csv(cluster_file, show_col_types = FALSE)
cat("Loaded annotations for", nrow(cluster_annotations), "clusters\n")

# Verify cluster annotations are complete
if (any(cluster_annotations$cell_type == "unassigned")) {
  warning("Some clusters are still marked as 'unassigned'")
}

cat("\n>>> Loading marker gene lists...\n")
if (!file.exists(markers_file)) {
  warning("Marker file not found: ", markers_file)
  cell_type_markers <- c("RBFOX3", "SLC17A7", "GAD1", "GFAP", "AIF1", "MBP", "PDGFRA")
  subtype_markers <- c("SATB2", "FEZF2", "PVALB", "SST", "VIP", "GJA1")
} else {
  marker_data <- read_csv(markers_file, show_col_types = FALSE)
  cell_type_markers <- marker_data %>%
    filter(level %in% c("cell_types", "neuronal_lineage")) %>%
    pull(gene)
  subtype_markers <- marker_data %>%
    filter(level %in% c("NPC_patterning", "cortical_layers", 
                         "Interneuron_subtypes", "Astrocyte_subtypes")) %>%
    pull(gene)
}

################################################################################
# Apply Cluster Labels
################################################################################

cat("\n>>> Applying cluster labels to dataset...\n")

# Set default clustering resolution (adjust based on your selection)
clustering_resolution <- 0.1
cluster_column <- paste0("SCT_snn_res.", clustering_resolution)

cat("Using clustering resolution:", clustering_resolution, "\n")

# Set active cluster identity
filtered_seurat$seurat_clusters <- filtered_seurat@meta.data[[cluster_column]]

# Add cluster annotations
filtered_seurat$cluster_name <- cluster_annotations$cluster_name[
  match(filtered_seurat$seurat_clusters, cluster_annotations$cluster)
]
filtered_seurat$cell_type <- cluster_annotations$cell_type[
  match(filtered_seurat$seurat_clusters, cluster_annotations$cluster)
]
filtered_seurat$cell_class <- cluster_annotations$cell_class[
  match(filtered_seurat$seurat_clusters, cluster_annotations$cluster)
]

cat("Applied labels:\n")
cat("  Cluster names:", length(unique(filtered_seurat$cluster_name)), "\n")
cat("  Cell types:", length(unique(filtered_seurat$cell_type)), "\n")
cat("  Cell classes:", length(unique(filtered_seurat$cell_class)), "\n")

# Save labeled dataset
labeled_file <- file.path(output_dir, "labeled_seurat.rda")
save(filtered_seurat, file = labeled_file)
cat("\nSaved labeled dataset to:", labeled_file, "\n")

################################################################################
# Generate Labeled UMAP Visualizations
################################################################################

cat("\n>>> Creating labeled UMAP visualizations...\n")

# Define grouping variables
groups <- unique(sample_metadata$group)
samples <- unique(sample_metadata$sample)
cluster_names <- unique(cluster_annotations$cluster_name)
cell_types <- unique(cluster_annotations$cell_type)
cell_classes <- unique(cluster_annotations$cell_class)

# Subsample for plotting if needed
full_seurat <- filtered_seurat
if (ncol(filtered_seurat) > 100000) {
  cat("Subsampling to 100,000 cells for visualization\n")
  plot_seurat <- filtered_seurat[, sample(colnames(filtered_seurat), 100000)]
} else {
  plot_seurat <- filtered_seurat
}

# Create UMAP plots
umap_plots <- list()

# By group
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

# By cluster number
umap_plots[["by_cluster"]] <- DimPlot(
  plot_seurat,
  group.by = "seurat_clusters",
  shuffle = TRUE,
  repel = TRUE,
  label = TRUE,
  reduction = "umap",
  pt.size = 0.01
) +
  scale_color_manual(values = generate_palette(plot_seurat$seurat_clusters)) +
  NoLegend() +
  labs(title = "Clusters")

# By cluster number (distinct colors)
umap_plots[["by_cluster_distinct"]] <- DimPlot(
  plot_seurat,
  group.by = "seurat_clusters",
  shuffle = TRUE,
  repel = TRUE,
  label = TRUE,
  reduction = "umap",
  pt.size = 0.01
) +
  scale_color_manual(values = generate_distinct_palette(plot_seurat$seurat_clusters)) +
  NoLegend() +
  labs(title = "Clusters (Distinct Colors)")

# By cluster name
umap_plots[["by_cluster_name"]] <- DimPlot(
  plot_seurat,
  group.by = "cluster_name",
  shuffle = TRUE,
  repel = TRUE,
  label = TRUE,
  reduction = "umap",
  pt.size = 0.01
) +
  scale_color_manual(limits = cluster_names, values = generate_palette(cluster_names)) +
  NoLegend() +
  labs(title = "Cluster Names")

# By cluster name (distinct colors)
umap_plots[["by_cluster_name_distinct"]] <- DimPlot(
  plot_seurat,
  group.by = "cluster_name",
  shuffle = TRUE,
  repel = TRUE,
  label = TRUE,
  reduction = "umap",
  pt.size = 0.01
) +
  scale_color_manual(limits = cluster_names, values = generate_distinct_palette(cluster_names)) +
  NoLegend() +
  labs(title = "Cluster Names (Distinct Colors)")

# By cell type
umap_plots[["by_cell_type"]] <- DimPlot(
  plot_seurat,
  group.by = "cell_type",
  shuffle = TRUE,
  repel = TRUE,
  label = TRUE,
  reduction = "umap",
  pt.size = 0.01
) +
  scale_color_manual(limits = cell_types, values = generate_palette(cell_types)) +
  NoLegend() +
  labs(title = "Cell Types")

# By cell type (distinct colors)
umap_plots[["by_cell_type_distinct"]] <- DimPlot(
  plot_seurat,
  group.by = "cell_type",
  shuffle = TRUE,
  repel = TRUE,
  label = TRUE,
  reduction = "umap",
  pt.size = 0.01
) +
  scale_color_manual(limits = cell_types, values = generate_distinct_palette(cell_types)) +
  NoLegend() +
  labs(title = "Cell Types (Distinct Colors)")

# By cell class
umap_plots[["by_cell_class"]] <- DimPlot(
  plot_seurat,
  group.by = "cell_class",
  shuffle = TRUE,
  repel = TRUE,
  label = TRUE,
  reduction = "umap",
  pt.size = 0.01
) +
  scale_color_manual(limits = cell_classes, values = generate_palette(cell_classes)) +
  NoLegend() +
  labs(title = "Cell Classes")

# By cell class (distinct colors)
umap_plots[["by_cell_class_distinct"]] <- DimPlot(
  plot_seurat,
  group.by = "cell_class",
  shuffle = TRUE,
  repel = TRUE,
  label = TRUE,
  reduction = "umap",
  pt.size = 0.01
) +
  scale_color_manual(limits = cell_classes, values = generate_distinct_palette(cell_classes)) +
  NoLegend() +
  labs(title = "Cell Classes (Distinct Colors)")

# Save UMAP plots
umap_file <- file.path(output_dir, "umap_labeled_clusters.pdf")
pdf(umap_file, width = 6, height = 5)
lapply(umap_plots, print)
dev.off()
cat("Saved UMAP plots to:", umap_file, "\n")

################################################################################
# Marker Expression Dot Plots
################################################################################

cat("\n>>> Creating marker expression dot plots...\n")

DefaultAssay(full_seurat) <- "SCT"

# Cell type markers
available_cell_markers <- intersect(cell_type_markers, rownames(full_seurat))
if (length(available_cell_markers) > 0) {
  dotplot_cell_types <- DotPlot(
    full_seurat,
    features = available_cell_markers,
    group.by = "cluster_name",
    scale.by = "size"
  ) +
    RotatedAxis() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_y_discrete(limits = cluster_names) +
    labs(title = "Cell Type Markers")
}

# Subtype markers
available_subtype_markers <- intersect(subtype_markers, rownames(full_seurat))
if (length(available_subtype_markers) > 0) {
  dotplot_subtypes <- DotPlot(
    full_seurat,
    features = available_subtype_markers,
    group.by = "cluster_name",
    scale.by = "size"
  ) +
    RotatedAxis() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_y_discrete(limits = cluster_names) +
    labs(title = "Subtype Markers")
}

# Save dot plots
dotplot_file <- file.path(output_dir, "dotplot_labeled_clusters.pdf")
pdf(dotplot_file, width = 10, height = 8)
if (exists("dotplot_cell_types")) print(dotplot_cell_types)
if (exists("dotplot_subtypes")) print(dotplot_subtypes)
dev.off()
cat("Saved dot plots to:", dotplot_file, "\n")

################################################################################
# Sample Integration Check
################################################################################

cat("\n>>> Creating sample comparison visualization...\n")

# Subsample if needed
if (ncol(full_seurat) > 100000) {
  plot_seurat <- full_seurat[, sample(colnames(full_seurat), 100000)]
} else {
  plot_seurat <- full_seurat
}

# Extract UMAP coordinates
umap_data <- FetchData(
  plot_seurat,
  vars = c("umap_1", "umap_2", "cluster_name", "group", "sample")
)

# Create faceted plot
sample_comparison <- ggplot(
  umap_data,
  aes(x = umap_1, y = umap_2, color = cluster_name)
) +
  geom_point(size = 0.05, alpha = 0.5) +
  scale_color_manual(limits = cluster_names, values = generate_distinct_palette(cluster_names)) +
  facet_wrap(~ factor(sample, levels = samples), nrow = 2) +
  theme_bw() +
  labs(title = "UMAP by Sample (Colored by Cluster)") +
  theme(legend.position = "right")

# Save sample comparison
sample_file <- file.path(output_dir, "umap_by_sample_labeled.pdf")
pdf(sample_file, width = 10, height = 5)
print(sample_comparison)
dev.off()
cat("Saved sample comparison to:", sample_file, "\n")

################################################################################
# Individual Marker Expression
################################################################################

cat("\n>>> Creating marker expression UMAPs...\n")

# Subsample for marker plots
if (ncol(full_seurat) > 10000) {
  plot_seurat <- full_seurat[, sample(colnames(full_seurat), 10000)]
} else {
  plot_seurat <- full_seurat
}

DefaultAssay(plot_seurat) <- "SCT"

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
  
  marker_file <- file.path(output_dir, "umap_marker_expression_labeled.pdf")
  pdf(marker_file, width = 30, height = 15)
  print(marker_plots)
  dev.off()
  cat("Saved marker expression UMAPs to:", marker_file, "\n")
}

################################################################################
# Cell Abundance Quantification
################################################################################

cat("\n>>> Quantifying cell abundance across samples and clusters...\n")

# Create abundance table
abundance_data <- tibble(
  cluster = rep(cluster_names, each = length(samples)),
  sample = rep(samples, length(cluster_names))
)

# Add group information
abundance_data$group <- sample_metadata$group[
  match(abundance_data$sample, sample_metadata$sample)
]

# Count cells per cluster per sample
cell_counts <- full_seurat@meta.data %>%
  group_by(cluster_name, sample) %>%
  summarise(n_cells = n(), .groups = "drop")

abundance_data$n_cells <- cell_counts$n_cells[
  match(
    paste0(abundance_data$cluster, abundance_data$sample),
    paste0(cell_counts$cluster_name, cell_counts$sample)
  )
]
abundance_data$n_cells[is.na(abundance_data$n_cells)] <- 0

# Calculate proportions
sample_totals <- abundance_data %>%
  group_by(sample) %>%
  summarise(total_cells = sum(n_cells))

abundance_data$total_sample <- sample_totals$total_cells[
  match(abundance_data$sample, sample_totals$sample)
]
abundance_data$fraction_of_sample <- abundance_data$n_cells / abundance_data$total_sample

cluster_totals <- abundance_data %>%
  group_by(cluster) %>%
  summarise(total_cells = sum(n_cells))

abundance_data$total_cluster <- cluster_totals$total_cells[
  match(abundance_data$cluster, cluster_totals$cluster)
]
abundance_data$fraction_of_cluster <- abundance_data$n_cells / abundance_data$total_cluster

# Save abundance data
abundance_file <- file.path(output_dir, "abundance_by_sample_cluster.csv")
write_csv(abundance_data, abundance_file)
cat("Saved abundance data to:", abundance_file, "\n")

################################################################################
# Abundance Visualizations
################################################################################

cat("\n>>> Creating abundance visualizations...\n")

# Bar plot with error bars
abundance_summary <- abundance_data %>%
  group_by(cluster, group) %>%
  summarise(
    mean_fraction = mean(fraction_of_sample),
    sd_fraction = sd(fraction_of_sample),
    .groups = "drop"
  )

abundance_barplot <- ggplot() +
  geom_col(
    data = abundance_summary,
    aes(x = cluster, y = mean_fraction, color = group),
    fill = "grey90",
    position = position_dodge(),
    width = 0.5,
    linewidth = 0.3
  ) +
  geom_errorbar(
    data = abundance_summary,
    aes(
      x = cluster,
      ymin = mean_fraction - sd_fraction,
      y = mean_fraction,
      ymax = mean_fraction + sd_fraction,
      color = group
    ),
    position = position_dodge(width = 0.5),
    width = 0.3,
    linewidth = 0.2
  ) +
  geom_point(
    data = abundance_data,
    aes(x = cluster, y = fraction_of_sample, color = group),
    position = position_dodge(width = 0.5),
    size = 0.5
  ) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits = cluster_names) +
  scale_color_manual(limits = groups, values = generate_palette(groups)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "Cell Abundance by Cluster and Group",
    x = "Cluster",
    y = "Fraction of Sample",
    color = "Group"
  )

# Save bar plot
barplot_file <- file.path(output_dir, "abundance_barplot.pdf")
pdf(barplot_file, width = 12, height = 7)
print(abundance_barplot)
dev.off()
cat("Saved abundance bar plot to:", barplot_file, "\n")

# Heatmap
abundance_matrix <- matrix(
  nrow = length(cluster_names),
  ncol = length(samples),
  dimnames = list(cluster_names, samples)
)

for (sample_id in samples) {
  sample_data <- abundance_data %>% filter(sample == sample_id)
  abundance_matrix[, sample_id] <- sample_data$fraction_of_sample[
    match(rownames(abundance_matrix), sample_data$cluster)
  ]
}
abundance_matrix[is.na(abundance_matrix)] <- 0

# Save heatmap
heatmap_file <- file.path(output_dir, "abundance_heatmap.pdf")
pdf(heatmap_file, width = 10, height = 8)
create_abundance_heatmap(abundance_matrix, title = "Cell Abundance by Sample and Cluster")
create_abundance_heatmap(
  abundance_matrix,
  title = "Cell Abundance (Clustered)",
  cluster_rows = TRUE,
  cluster_cols = TRUE
)
dev.off()
cat("Saved abundance heatmap to:", heatmap_file, "\n")

################################################################################
# Statistical Testing
################################################################################

cat("\n>>> Performing differential abundance t-tests...\n")

# T-test for each cluster
ttest_results <- tibble(
  cluster = cluster_names,
  p_value = NA_real_,
  p_adjusted = NA_real_
)

for (cluster_id in cluster_names) {
  cluster_data <- abundance_data %>% filter(cluster == cluster_id)
  
  if (length(unique(cluster_data$group)) >= 2) {
    test_result <- pairwise.t.test(
      cluster_data$fraction_of_sample,
      cluster_data$group,
      p.adjust.method = "none"
    )
    ttest_results$p_value[ttest_results$cluster == cluster_id] <- test_result$p.value[1, 1]
  }
}

# Adjust p-values
ttest_results$p_adjusted <- p.adjust(ttest_results$p_value, method = "BH")

# Save t-test results
ttest_file <- file.path(output_dir, "abundance_ttest.csv")
write_csv(ttest_results, ttest_file)
cat("Saved t-test results to:", ttest_file, "\n")

################################################################################
# sccomp Differential Abundance Analysis
################################################################################

cat("\n>>> Running sccomp differential abundance analysis...\n")

sccomp_results <- full_seurat %>%
  sccomp_estimate(
    formula_composition = ~ group,
    .sample = sample,
    .cell_group = cluster_name,
    bimodal_mean_variability_association = TRUE,
    cores = 1,
    verbose = FALSE
  ) %>%
  sccomp_test()

# Save sccomp results
sccomp_file <- file.path(output_dir, "sccomp_results.csv")
write_csv(sccomp_results, sccomp_file)
cat("Saved sccomp results to:", sccomp_file, "\n")

# Create sccomp visualizations
sccomp_plot_file <- file.path(output_dir, "sccomp_plots.pdf")
pdf(sccomp_plot_file, width = 8, height = 8)

# Interval plot
print(sccomp_results %>% plot_1D_intervals())

# Boxplot
print(sccomp_results %>% sccomp_boxplot(factor = "group"))

dev.off()
cat("Saved sccomp plots to:", sccomp_plot_file, "\n")

################################################################################
# Session Info
################################################################################

cat("\n>>> Recording session information...\n")
session_file <- file.path(output_dir, "session_info.txt")
writeLines(capture.output(sessionInfo()), session_file)
cat("Saved session info to:", session_file, "\n")

################################################################################
# Completion Summary
################################################################################

cat("\n")
cat("========================================================================\n")
cat("C04: Cluster Characterization Complete\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n")
cat("\nSummary:\n")
cat("  Total cells:", ncol(full_seurat), "\n")
cat("  Clusters:", length(cluster_names), "\n")
cat("  Cell types:", length(unique(full_seurat$cell_type)), "\n")
cat("  Samples:", length(samples), "\n")
cat("  Groups:", length(groups), "\n")
cat("\nDifferential abundance:\n")
cat("  Significant clusters (p < 0.05):", 
    sum(ttest_results$p_adjusted < 0.05, na.rm = TRUE), "\n")
cat("\nOutput files:\n")
cat("  -", labeled_file, "\n")
cat("  -", abundance_file, "\n")
cat("  -", ttest_file, "\n")
cat("  -", sccomp_file, "\n")
cat("========================================================================\n\n")
