#!/usr/bin/env Rscript

################################################################################
# Script: C05_sctype_annotation.R
# Description: Automated cell type annotation using scType
#
# This script:
#   1. Loads the scType annotation framework from GitHub
#   2. Prepares brain tissue-specific marker gene sets
#   3. Performs automated cell type classification
#   4. Generates comprehensive visualization suite
#   5. Extracts astrocyte population for downstream analysis
#
# Inputs:
#   - results/04_characterization/labeled_seurat.rda: Labeled dataset
#   - results/02_integration/filtered_metadata.csv: Sample metadata
#   - scType reference database (loaded from GitHub)
#
# Outputs:
#   - results/05_sctype/sctype_annotated_seurat.rda: Dataset with scType labels
#   - results/05_sctype/astrocytes_subset.rda: Astrocyte-only dataset
#   - results/05_sctype/sctype_umap.pdf: UMAP with cell type annotations
#   - results/05_sctype/cell_type_distribution.pdf: Bar plot of cell types
#   - results/05_sctype/marker_expression_heatmap.pdf: Marker expression dot plot
#   - results/05_sctype/violin_plots.pdf: Key marker violin plots
#   - results/05_sctype/feature_plots.pdf: Individual marker UMAPs
#   - results/05_sctype/cluster_composition.pdf: Cell type by cluster
#   - results/05_sctype/sctype_complete_report.pdf: Comprehensive PDF report
#   - results/05_sctype/cell_type_summary.csv: Cell type statistics
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(colorRamps)
  library(viridis)
  library(patchwork)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(HGNChelper)
  library(openxlsx)
})

# Initialize logging
cat("\n")
cat("========================================================================\n")
cat("C05: scType Automated Annotation\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n\n")

# Set random seed for reproducibility
set.seed(1234)

# Define paths
project_dir <- getwd()
input_file <- "results/04_characterization/labeled_seurat.rda"
metadata_file <- "results/02_integration/filtered_metadata.csv"
output_dir <- "results/05_sctype"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

################################################################################
# Load scType Functions
################################################################################

cat("\n>>> Loading scType framework from GitHub...\n")

# Load scType gene sets preparation function
sctype_url_1 <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"
source(sctype_url_1)

# Load scType scoring function
sctype_url_2 <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R"
source(sctype_url_2)

cat("scType functions loaded successfully\n")

################################################################################
# Define Color Palette Functions
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

# Distinct color palette
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

################################################################################
# Load Data
################################################################################

cat("\n>>> Loading labeled dataset...\n")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}
load(input_file)
cat("Loaded dataset with", ncol(filtered_seurat), "cells\n")

cat("\n>>> Loading sample metadata...\n")
sample_metadata <- read_csv(metadata_file, show_col_types = FALSE)

################################################################################
# Prepare Brain-Specific Marker Gene Sets
################################################################################

cat("\n>>> Preparing brain tissue marker gene sets...\n")

# Load scType database for brain tissue
db_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
cat("Loading scType database from:", db_url, "\n")

sctype_db <- gene_sets_prepare(db_url, "Brain")
cat("Loaded scType brain markers\n")

# Define custom brain marker sets (more comprehensive)
brain_markers <- list(
  gs_positive = list(
    "Excitatory neurons" = c("SLC17A7", "CAMK2A", "RBFOX3", "NRGN", "NEFL", 
                              "SNAP25", "SYT1", "STMN2", "THY1"),
    "Inhibitory neurons" = c("GAD1", "GAD2", "SLC32A1", "PVALB", "SST", "VIP"),
    "Astrocytes" = c("GFAP", "AQP4", "ALDH1L1", "S100B", "GLUL", "SLC1A3"),
    "Microglia" = c("AIF1", "CD68", "CX3CR1", "TMEM119", "P2RY12", "CSF1R"),
    "Oligodendrocytes" = c("MBP", "MOG", "OLIG2", "PLP1", "MAG", "MOBP"),
    "OPC" = c("PDGFRA", "CSPG4", "SOX10", "OLIG2"),
    "Endothelial" = c("PECAM1", "VWF", "CLDN5", "FLT1"),
    "Pericytes" = c("PDGFRB", "RGS5", "NOTCH3", "ACTA2")
  ),
  gs_negative = list(
    "Excitatory neurons" = c("GFAP", "AIF1", "MBP", "PECAM1"),
    "Inhibitory neurons" = c("GFAP", "AIF1", "MBP", "PECAM1"),
    "Astrocytes" = c("SLC17A7", "AIF1", "MBP", "PECAM1"),
    "Microglia" = c("GFAP", "SLC17A7", "MBP", "PECAM1"),
    "Oligodendrocytes" = c("GFAP", "AIF1", "SLC17A7", "PECAM1"),
    "OPC" = c("GFAP", "AIF1", "MBP", "PECAM1"),
    "Endothelial" = c("GFAP", "AIF1", "SLC17A7", "MBP"),
    "Pericytes" = c("GFAP", "AIF1", "SLC17A7", "MBP")
  )
)

cat("Custom brain marker sets prepared\n")

################################################################################
# Run scType Classification
################################################################################

cat("\n>>> Running scType classification...\n")

# Extract scaled data
DefaultAssay(filtered_seurat) <- "SCT"
scaled_data <- as.matrix(filtered_seurat[["SCT"]]@scale.data)
cat("Extracted scaled expression matrix:", 
    nrow(scaled_data), "genes x", ncol(scaled_data), "cells\n")

# Calculate scType scores
cat("Calculating cell type scores...\n")
sctype_scores <- sctype_score(
  scRNAseqData = scaled_data,
  scaled = TRUE,
  gs = brain_markers$gs_positive,
  gs2 = brain_markers$gs_negative
)

cat("Scores calculated for", nrow(sctype_scores), "cell types\n")

# Assign cell types per cluster
cat("Assigning cell types to clusters...\n")

# Use resolution 0.1 clustering (adjust if different)
cluster_column <- "SCT_snn_res.0.1"
clusters <- unique(filtered_seurat@meta.data[[cluster_column]])

classification_results <- do.call("rbind", lapply(clusters, function(cluster_id) {
  # Get cells in this cluster
  cluster_cells <- rownames(filtered_seurat@meta.data[
    filtered_seurat@meta.data[[cluster_column]] == cluster_id,
  ])
  
  # Sum scores for cells in cluster
  cluster_scores <- sort(
    rowSums(sctype_scores[, cluster_cells, drop = FALSE]),
    decreasing = TRUE
  )
  
  # Return top cell types for this cluster
  head(data.frame(
    cluster = cluster_id,
    cell_type = names(cluster_scores),
    score = cluster_scores,
    n_cells = length(cluster_cells)
  ), 10)
}))

# Select top-scoring cell type per cluster
cluster_assignments <- classification_results %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = score)

cat("Cell type assignments:\n")
print(cluster_assignments)

# Add scType classification to dataset
filtered_seurat$sctype_classification <- ""
for (cluster_id in unique(cluster_assignments$cluster)) {
  assigned_type <- cluster_assignments %>%
    filter(cluster == cluster_id) %>%
    pull(cell_type) %>%
    first()
  
  filtered_seurat$sctype_classification[
    filtered_seurat@meta.data[[cluster_column]] == cluster_id
  ] <- as.character(assigned_type)
}

cat("\nscType classification complete\n")
cat("Cell types identified:", 
    length(unique(filtered_seurat$sctype_classification)), "\n")

# Save annotated dataset
annotated_file <- file.path(output_dir, "sctype_annotated_seurat.rda")
save(filtered_seurat, file = annotated_file)
cat("Saved annotated dataset to:", annotated_file, "\n")

################################################################################
# Extract Astrocyte Subset
################################################################################

cat("\n>>> Extracting astrocyte population...\n")

astrocytes <- filtered_seurat[, filtered_seurat$sctype_classification == "Astrocytes"]
cat("Astrocyte cells:", ncol(astrocytes), "\n")

if (ncol(astrocytes) > 0) {
  astrocyte_file <- file.path(output_dir, "astrocytes_subset.rda")
  save(astrocytes, file = astrocyte_file)
  cat("Saved astrocyte subset to:", astrocyte_file, "\n")
}

################################################################################
# Generate Visualizations
################################################################################

cat("\n>>> Generating visualization suite...\n")

# Generate color palettes
cell_type_colors <- generate_distinct_palette(filtered_seurat$sctype_classification)
cluster_colors <- generate_palette(filtered_seurat@meta.data[[cluster_column]])

# 1. UMAP with scType annotations
cat("Creating UMAP plots...\n")
umap_sctype <- DimPlot(
  filtered_seurat,
  reduction = "umap",
  group.by = "sctype_classification",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5,
  cols = cell_type_colors
) +
  ggtitle("scType Cell Type Annotations") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

# 2. UMAP with clusters
umap_clusters <- DimPlot(
  filtered_seurat,
  reduction = "umap",
  group.by = cluster_column,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5,
  cols = cluster_colors
) +
  ggtitle("Clusters") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

# 3. Combined UMAP plot
umap_combined <- umap_sctype | umap_clusters

# Save UMAP plots
umap_file <- file.path(output_dir, "sctype_umap.pdf")
ggsave(umap_file, umap_combined, width = 12, height = 6, dpi = 300)
cat("Saved UMAP plots to:", umap_file, "\n")

# 4. Cell type distribution bar plot
cat("Creating cell type distribution plot...\n")
cell_type_counts <- as.data.frame(table(filtered_seurat$sctype_classification))
colnames(cell_type_counts) <- c("Cell_Type", "Count")

bar_colors <- generate_distinct_palette(cell_type_counts$Cell_Type)

distribution_plot <- ggplot(
  cell_type_counts,
  aes(x = reorder(Cell_Type, Count), y = Count, fill = Cell_Type)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = bar_colors) +
  coord_flip() +
  labs(title = "Cell Type Distribution", x = "Cell Type", y = "Number of Cells") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 10)
  )

# Save distribution plot
distribution_file <- file.path(output_dir, "cell_type_distribution.pdf")
ggsave(distribution_file, distribution_plot, width = 8, height = 6, dpi = 300)
cat("Saved distribution plot to:", distribution_file, "\n")

# 5. Marker expression heatmap
cat("Creating marker expression heatmap...\n")
key_markers <- c("NRGN", "NEFL", "SNAP25", "GFAP", "AQP4", "AIF1", 
                  "CX3CR1", "MBP", "OLIG2")
available_markers <- intersect(key_markers, rownames(filtered_seurat))

if (length(available_markers) > 0) {
  marker_heatmap <- DotPlot(
    filtered_seurat,
    features = available_markers,
    group.by = "sctype_classification",
    cols = c("lightgrey", "red")
  ) +
    RotatedAxis() +
    ggtitle("Key Marker Expression by Cell Type") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  heatmap_file <- file.path(output_dir, "marker_expression_heatmap.pdf")
  ggsave(heatmap_file, marker_heatmap, width = 10, height = 6, dpi = 300)
  cat("Saved marker heatmap to:", heatmap_file, "\n")
}

# 6. Violin plots for astrocyte markers
cat("Creating violin plots for key markers...\n")
astrocyte_markers <- c("GJA1", "GFAP", "AQP4", "S100B")
available_astro_markers <- intersect(astrocyte_markers, rownames(filtered_seurat))

if (length(available_astro_markers) > 0) {
  violin_plots <- VlnPlot(
    filtered_seurat,
    features = available_astro_markers,
    group.by = "sctype_classification",
    cols = cell_type_colors,
    ncol = 2
  ) +
    theme(legend.position = "none")
  
  violin_file <- file.path(output_dir, "violin_plots.pdf")
  ggsave(violin_file, violin_plots, width = 12, height = 8, dpi = 300)
  cat("Saved violin plots to:", violin_file, "\n")
}

# 7. Feature plots for key markers
cat("Creating feature plots...\n")
if (length(available_astro_markers) > 0) {
  feature_plots <- FeaturePlot(
    filtered_seurat,
    features = available_astro_markers,
    ncol = 2,
    pt.size = 0.3
  ) &
    theme_minimal() &
    theme(plot.title = element_text(size = 12, face = "bold"))
  
  feature_file <- file.path(output_dir, "feature_plots.pdf")
  ggsave(feature_file, feature_plots, width = 10, height = 8, dpi = 300)
  cat("Saved feature plots to:", feature_file, "\n")
}

# 8. Cluster composition by cell type
cat("Creating cluster composition plot...\n")
cluster_composition <- table(
  filtered_seurat@meta.data[[cluster_column]],
  filtered_seurat$sctype_classification
)
composition_df <- as.data.frame.matrix(cluster_composition)
composition_df$Cluster <- rownames(composition_df)

composition_long <- composition_df %>%
  pivot_longer(cols = -Cluster, names_to = "Cell_Type", values_to = "Count") %>%
  group_by(Cluster) %>%
  mutate(Percentage = Count / sum(Count) * 100)

composition_plot <- ggplot(
  composition_long,
  aes(x = Cluster, y = Percentage, fill = Cell_Type)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_colors) +
  labs(
    title = "Cell Type Composition by Cluster",
    x = "Cluster",
    y = "Percentage (%)",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

composition_file <- file.path(output_dir, "cluster_composition.pdf")
ggsave(composition_file, composition_plot, width = 10, height = 6, dpi = 300)
cat("Saved composition plot to:", composition_file, "\n")

################################################################################
# Generate Summary Statistics
################################################################################

cat("\n>>> Generating summary statistics...\n")

cell_type_summary <- as.data.frame(table(filtered_seurat$sctype_classification))
colnames(cell_type_summary) <- c("Cell_Type", "Number_of_Cells")
cell_type_summary$Percentage <- round(
  (cell_type_summary$Number_of_Cells / sum(cell_type_summary$Number_of_Cells)) * 100,
  1
)
cell_type_summary <- cell_type_summary[
  order(cell_type_summary$Number_of_Cells, decreasing = TRUE),
]

# Save summary
summary_file <- file.path(output_dir, "cell_type_summary.csv")
write_csv(cell_type_summary, summary_file)
cat("Saved summary statistics to:", summary_file, "\n")

print(cell_type_summary)

################################################################################
# Create Comprehensive PDF Report
################################################################################

cat("\n>>> Creating comprehensive PDF report...\n")

report_file <- file.path(output_dir, "sctype_complete_report.pdf")
pdf(report_file, width = 12, height = 8)

# Page 1: Combined UMAP
print(umap_combined)

# Page 2: Cell type distribution
print(distribution_plot)

# Page 3: Marker heatmap
if (exists("marker_heatmap")) print(marker_heatmap)

# Page 4: Violin plots
if (exists("violin_plots")) print(violin_plots)

# Page 5: Feature plots
if (exists("feature_plots")) print(feature_plots)

# Page 6: Cluster composition
print(composition_plot)

# Page 7: Summary table
summary_grob <- tableGrob(cell_type_summary)
grid.arrange(summary_grob, top = "Cell Type Summary Statistics")

dev.off()
cat("Saved comprehensive report to:", report_file, "\n")

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
cat("C05: scType Annotation Complete\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n")
cat("\nSummary:\n")
cat("  Total cells:", ncol(filtered_seurat), "\n")
cat("  Cell types identified:", length(unique(filtered_seurat$sctype_classification)), "\n")
cat("  Astrocytes:", ncol(astrocytes), "\n")
cat("\nCell type distribution:\n")
print(cell_type_summary)
cat("\nOutput files:\n")
cat("  -", annotated_file, "\n")
if (ncol(astrocytes) > 0) cat("  -", astrocyte_file, "\n")
cat("  -", report_file, "\n")
cat("========================================================================\n\n")
