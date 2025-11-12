#!/usr/bin/env Rscript

################################################################################
# Script: C06_astrocyte_integration.R
# Description: Re-integrate and normalize astrocyte-only population
#
# This script:
#   1. Loads astrocyte subset from scType annotation (C05)
#   2. Re-normalizes astrocyte-specific data with SCTransform
#   3. Re-integrates astrocytes across samples using Harmony
#   4. Generates new UMAP embedding for astrocyte population
#   5. Prepares for astrocyte-specific clustering
#
# Inputs:
#   - results/05_sctype/astrocytes_subset.rda: Astrocyte-only dataset
#   - results/02_integration/filtered_metadata.csv: Sample metadata
#
# Outputs:
#   - results/06_astrocyte_analysis/integrated_astrocytes.rda: Integrated astrocytes
#   - Session info and progress logs
#
# Note: This script performs a fresh integration specifically for astrocytes
#       to better resolve astrocyte subtypes and states
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
cat("C06: Astrocyte Re-Integration\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n\n")

# Set random seed for reproducibility
set.seed(1234)

# Define paths
project_dir <- getwd()
input_file <- "results/05_sctype/astrocytes_subset.rda"
metadata_file <- "results/02_integration/filtered_metadata.csv"
output_dir <- "results/06_astrocyte_analysis"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

################################################################################
# Load Data
################################################################################

cat("\n>>> Loading astrocyte subset...\n")

if (!file.exists(input_file)) {
  stop("Astrocyte subset not found: ", input_file,
       "\nPlease run C05_sctype_annotation.R first")
}

load(input_file)
cat("Loaded astrocyte dataset with", ncol(astrocytes), "cells\n")

cat("\n>>> Loading sample metadata...\n")
sample_metadata <- read_csv(metadata_file, show_col_types = FALSE)

################################################################################
# Re-normalize Astrocyte Data
################################################################################

cat("\n>>> Performing SCTransform normalization on astrocytes...\n")

# Increase memory limit for SCTransform
options(future.globals.maxSize = 1000 * 1024^2)

# Split RNA assay by sample for layer-wise processing
cat("Splitting RNA layers by sample...\n")
astrocytes[["RNA"]] <- split(astrocytes[["RNA"]], f = astrocytes$sample)

# Perform SCTransform normalization
cat("Running SCTransform normalization...\n")
astrocytes <- SCTransform(
  astrocytes,
  assay = "RNA",
  new.assay.name = "SCT",
  ncells = 3000,
  variable.features.n = 2000,
  conserve.memory = TRUE,
  verbose = TRUE
)

cat("SCTransform normalization complete\n")

# Save intermediate result
intermediate_file <- file.path(output_dir, "astrocytes_normalized_intermediate.rda")
save(astrocytes, file = intermediate_file)
cat("Saved intermediate normalization to:", intermediate_file, "\n")

################################################################################
# PCA and Harmony Integration
################################################################################

cat("\n>>> Running PCA on astrocyte data...\n")
astrocytes <- RunPCA(
  astrocytes,
  assay = "SCT",
  reduction.name = "pca",
  verbose = FALSE
)

cat("\n>>> Running Harmony integration on astrocytes...\n")
astrocytes <- IntegrateLayers(
  object = astrocytes,
  method = HarmonyIntegration,
  assay = "SCT",
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = TRUE
)

cat("Harmony integration complete\n")

# Save post-integration
integration_file <- file.path(output_dir, "astrocytes_harmony_integrated.rda")
save(astrocytes, file = integration_file)
cat("Saved Harmony integration to:", integration_file, "\n")

################################################################################
# Join Layers and Final Normalization
################################################################################

cat("\n>>> Joining RNA layers and performing final normalization...\n")

# Increase memory limit
options(future.globals.maxSize = 30 * 1024^3)

# Join RNA layers
cat("Joining RNA layers...\n")
astrocytes <- JoinLayers(astrocytes, assay = "RNA")

# Re-run SCTransform on joined data
cat("Running final SCTransform normalization...\n")
astrocytes <- SCTransform(
  astrocytes,
  assay = "RNA",
  new.assay.name = "SCT",
  ncells = 3000,
  variable.features.n = 2000,
  conserve.memory = TRUE,
  verbose = TRUE
)

# Re-run PCA
cat("Running final PCA...\n")
astrocytes <- RunPCA(
  astrocytes,
  assay = "SCT",
  reduction.name = "pca",
  verbose = FALSE
)

################################################################################
# UMAP and Clustering Preparation
################################################################################

cat("\n>>> Generating UMAP embedding for astrocytes...\n")
astrocytes <- RunUMAP(
  astrocytes,
  reduction = "harmony",
  dims = 1:30,
  return.model = TRUE,
  verbose = FALSE
)

cat("\n>>> Computing nearest neighbor graph...\n")
astrocytes <- FindNeighbors(
  astrocytes,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

cat("Nearest neighbor graph computed\n")

################################################################################
# Save Final Integrated Astrocyte Dataset
################################################################################

cat("\n>>> Saving final integrated astrocyte dataset...\n")

output_file <- file.path(output_dir, "integrated_astrocytes.rda")
save(astrocytes, file = output_file)
cat("Saved to:", output_file, "\n")

# Garbage collection
gc(verbose = FALSE, full = TRUE)

################################################################################
# Completion Summary
################################################################################

cat("\n")
cat("========================================================================\n")
cat("C06: Astrocyte Re-Integration Complete\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n")
cat("\nSummary:\n")
cat("  Input astrocytes:", ncol(astrocytes), "\n")
cat("  Samples:", length(unique(astrocytes$sample)), "\n")
cat("  Genes:", nrow(astrocytes), "\n")
cat("  Integration method: Harmony\n")
cat("  Harmony dimensions: 30\n")
cat("  Ready for clustering: Yes\n")
cat("\nOutput files:\n")
cat("  -", output_file, "\n")
cat("\nNext steps:\n")
cat("  1. Run C07_astrocyte_clustering.R to identify astrocyte subtypes\n")
cat("  2. Review clustering at multiple resolutions\n")
cat("  3. Annotate astrocyte clusters based on marker expression\n")
cat("========================================================================\n\n")
