#!/usr/bin/env Rscript

################################################################################
# Script: C02_harmony_integration.R
# Description: Filter low-quality cells and integrate samples with Harmony
#
# This script:
#   1. Filters cells based on QC metrics (UMI counts, mitochondrial %)
#   2. Removes samples with excessive cell loss or insufficient cells
#   3. Normalizes data using SCTransform
#   4. Integrates samples using Harmony batch correction
#   5. Performs dimensionality reduction (PCA, UMAP) and clustering
#
# Inputs:
#   - results/01_data_loading/merged_seurat.rda: Merged unfiltered dataset
#   - data/metadata/sample_metadata.csv: Sample information
#
# Outputs:
#   - results/02_integration/filter_statistics.csv: Cell filtering summary
#   - results/02_integration/filtered_metadata.csv: Metadata for retained samples
#   - results/02_integration/integrated_seurat.rda: Harmony-integrated dataset
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(harmony)
  library(future)
})

# Initialize logging
cat("\n")
cat("========================================================================\n")
cat("C02: Harmony Integration\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n\n")

# Set random seed for reproducibility
set.seed(1234)

# Define paths
project_dir <- getwd()
input_file <- "results/01_data_loading/merged_seurat.rda"
metadata_file <- "data/metadata/sample_metadata.csv"
output_dir <- "results/02_integration"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

################################################################################
# Load Data
################################################################################

cat("\n>>> Loading merged dataset...\n")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

load(input_file)
cat("Loaded dataset with", ncol(merged_seurat), "cells\n")

cat("\n>>> Loading sample metadata...\n")
sample_metadata <- read_csv(metadata_file, show_col_types = FALSE)

################################################################################
# Filter Low-Quality Cells
################################################################################

cat("\n>>> Filtering low-quality cells...\n")

# Set default assay
DefaultAssay(merged_seurat) <- "RNA"

# Calculate pre-filter statistics
pre_filter_stats <- merged_seurat@meta.data %>%
  group_by(sample) %>%
  summarise(cells_unfiltered = n())

cat("Initial cell counts per sample:\n")
print(pre_filter_stats)

# Apply quality filters
cat("\nApplying quality filters:\n")
cat("  - Total UMI counts: 500 - 30,000\n")
cat("  - Mitochondrial percentage: < 10%\n\n")

filtered_seurat <- subset(
  x = merged_seurat,
  subset = nCount_RNA > 500 &
           nCount_RNA < 30000 &
           percent_mt < 10
)

# Calculate post-filter statistics
post_filter_stats <- filtered_seurat@meta.data %>%
  group_by(sample) %>%
  summarise(cells_filtered = n())

# Combine filter statistics
filter_stats <- pre_filter_stats %>%
  left_join(post_filter_stats, by = "sample") %>%
  mutate(
    cells_filtered = replace_na(cells_filtered, 0),
    cells_removed = cells_unfiltered - cells_filtered,
    fraction_removed = cells_removed / cells_unfiltered
  )

cat("Post-filter cell counts:\n")
print(filter_stats)

################################################################################
# Remove Low-Quality Samples
################################################################################

cat("\n>>> Identifying samples to retain...\n")
cat("Criteria:\n")
cat("  - Fraction removed < 50%\n")
cat("  - Remaining cells >= 1,000\n\n")

samples_to_retain <- filter_stats %>%
  filter(
    fraction_removed <= 0.5,
    cells_filtered >= 1000,
    sample %in% sample_metadata$sample
  ) %>%
  pull(sample)

cat("Retaining", length(samples_to_retain), "of", nrow(filter_stats), "samples\n")

# Update filter statistics
filter_stats <- filter_stats %>%
  mutate(cells_retained = if_else(sample %in% samples_to_retain, cells_filtered, 0))

# Save filter statistics
filter_stats_file <- file.path(output_dir, "filter_statistics.csv")
write_csv(filter_stats, filter_stats_file)
cat("Saved filter statistics to:", filter_stats_file, "\n")

# Filter dataset to retained samples
filtered_seurat <- subset(
  x = filtered_seurat,
  subset = sample %in% samples_to_retain
)

cat("Final dataset:", ncol(filtered_seurat), "cells from", 
    length(unique(filtered_seurat$sample)), "samples\n")

# Update and save filtered metadata
filtered_metadata <- sample_metadata %>%
  filter(sample %in% samples_to_retain)

filtered_metadata_file <- file.path(output_dir, "filtered_metadata.csv")
write_csv(filtered_metadata, filtered_metadata_file)
cat("Saved filtered metadata to:", filtered_metadata_file, "\n")

# Garbage collection
gc(verbose = FALSE, full = TRUE)

################################################################################
# SCTransform Normalization with Sample Splitting
################################################################################

cat("\n>>> Performing SCTransform normalization (layer-wise)...\n")

# Increase memory limit for SCTransform
options(future.globals.maxSize = 1000 * 1024^2)

# Split RNA assay by sample for layer-wise processing
cat("Splitting RNA layers by sample...\n")
filtered_seurat[["RNA"]] <- split(filtered_seurat[["RNA"]], f = filtered_seurat$sample)

# Perform SCTransform
cat("Running SCTransform normalization...\n")
filtered_seurat <- SCTransform(
  filtered_seurat,
  assay = "RNA",
  new.assay.name = "SCT",
  ncells = 3000,
  variable.features.n = 2000,
  conserve.memory = TRUE,
  verbose = TRUE
)

cat("SCTransform normalization complete\n")

################################################################################
# PCA and Harmony Integration
################################################################################

cat("\n>>> Running PCA...\n")
filtered_seurat <- RunPCA(
  filtered_seurat,
  assay = "SCT",
  reduction.name = "pca",
  verbose = FALSE
)

cat("\n>>> Running Harmony integration...\n")
filtered_seurat <- IntegrateLayers(
  object = filtered_seurat,
  method = HarmonyIntegration,
  assay = "SCT",
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = TRUE
)

cat("Harmony integration complete\n")

# Save intermediate result
intermediate_file <- file.path(output_dir, "integrated_seurat_intermediate.rda")
save(filtered_seurat, file = intermediate_file)
cat("Saved intermediate result to:", intermediate_file, "\n")

################################################################################
# Join Layers and Final Normalization
################################################################################

cat("\n>>> Joining RNA layers and performing final normalization...\n")

# Increase memory limit
options(future.globals.maxSize = 30 * 1024^3)

# Join RNA layers
cat("Joining RNA layers...\n")
filtered_seurat <- JoinLayers(filtered_seurat, assay = "RNA")

# Re-run SCTransform on joined data
cat("Running final SCTransform normalization...\n")
filtered_seurat <- SCTransform(
  filtered_seurat,
  assay = "RNA",
  new.assay.name = "SCT",
  ncells = 3000,
  variable.features.n = 2000,
  conserve.memory = TRUE,
  verbose = TRUE
)

# Re-run PCA
cat("Running final PCA...\n")
filtered_seurat <- RunPCA(
  filtered_seurat,
  assay = "SCT",
  reduction.name = "pca",
  verbose = FALSE
)

################################################################################
# UMAP and Clustering Preparation
################################################################################

cat("\n>>> Generating UMAP embedding...\n")
filtered_seurat <- RunUMAP(
  filtered_seurat,
  reduction = "harmony",
  dims = 1:30,
  return.model = TRUE,
  verbose = FALSE
)

cat("\n>>> Computing nearest neighbor graph...\n")
filtered_seurat <- FindNeighbors(
  filtered_seurat,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

cat("Nearest neighbor graph computed\n")

################################################################################
# Save Final Integrated Dataset
################################################################################

cat("\n>>> Saving final integrated dataset...\n")

output_file <- file.path(output_dir, "integrated_seurat.rda")
save(filtered_seurat, file = output_file)
cat("Saved to:", output_file, "\n")

# Garbage collection
gc(verbose = FALSE, full = TRUE)

################################################################################
# Completion Summary
################################################################################

cat("\n")
cat("========================================================================\n")
cat("C02: Harmony Integration Complete\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n")
cat("\nSummary:\n")
cat("  Input cells:", nrow(pre_filter_stats %>% summarise(total = sum(cells_unfiltered))), "\n")
cat("  Cells after QC filtering:", ncol(filtered_seurat), "\n")
cat("  Samples retained:", length(unique(filtered_seurat$sample)), "\n")
cat("  Genes:", nrow(filtered_seurat), "\n")
cat("  PCA dimensions computed: 50\n")
cat("  Harmony dimensions: 30\n")
cat("\nOutput files:\n")
cat("  -", filter_stats_file, "\n")
cat("  -", filtered_metadata_file, "\n")
cat("  -", output_file, "\n")
cat("========================================================================\n\n")
