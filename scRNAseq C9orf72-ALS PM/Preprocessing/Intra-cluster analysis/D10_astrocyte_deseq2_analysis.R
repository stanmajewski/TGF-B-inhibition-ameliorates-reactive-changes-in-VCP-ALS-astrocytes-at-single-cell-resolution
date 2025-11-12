#!/usr/bin/env Rscript

################################################################################
# Script: D10_astrocyte_deseq2_analysis.R
# Description: DESeq2 differential expression analysis on astrocyte pseudobulks
#
# This script:
#   1. Loads pseudobulk data from D09
#   2. Filters clusters with sufficient replication (≥2 pseudobulks per group)
#   3. Removes lowly expressed genes (<0.1 counts/cell)
#   4. Runs DESeq2 Wald test with cluster_group design
#   5. Extracts differential expression results for each cluster
#   6. Identifies significant DEGs (padj < 0.05, |log2FC| > log2(1.5))
#   7. Separates upregulated and downregulated genes
#   8. Annotates transcription factors
#
# Inputs:
#   - results/09_astrocyte_de/pseudobulk_data.rda: Pseudobulk counts
#   - data/reference/transcription_factors.csv: TF annotations (optional)
#
# Outputs:
#   - results/09_astrocyte_de/deseq2_dataset.rda: DESeq2 results object
#   - results/09_astrocyte_de/deg_lists.csv: DEG lists by cluster
#   - results/09_astrocyte_de/deg_counts.csv: Number of DEGs per cluster
#
# Note: Uses log2(1.5) = 0.58 fold-change threshold, common in RNA-seq
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
})

# Initialize logging
cat("\n")
cat("========================================================================\n")
cat("D10: Astrocyte DESeq2 Differential Expression Analysis\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n\n")

# Set random seed for reproducibility
set.seed(12345)

# Define paths
project_dir <- getwd()
input_file <- "results/09_astrocyte_de/pseudobulk_data.rda"
tf_file <- "data/reference/transcription_factors.csv"
output_dir <- "results/09_astrocyte_de"

################################################################################
# Load Data
################################################################################

cat("\n>>> Loading pseudobulk data...\n")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file,
       "\nPlease run D09_astrocyte_pseudobulk_generation.R first")
}
load(input_file)

cat("Loaded pseudobulk data:\n")
cat("  Pseudobulks:", ncol(pseudobulk_data$counts), "\n")
cat("  Genes:", nrow(pseudobulk_data$counts), "\n")

# Load transcription factor annotations (optional)
if (file.exists(tf_file)) {
  tf_data <- read_csv(tf_file, show_col_types = FALSE)
  transcription_factors <- unique(tf_data$Symbol)
  cat("\nLoaded", length(transcription_factors), "transcription factors\n")
} else {
  cat("\nTranscription factor file not found, will skip TF annotation\n")
  transcription_factors <- character(0)
}

################################################################################
# Prepare Data for DESeq2
################################################################################

cat("\n>>> Preparing data for DESeq2 analysis...\n")

# Define comparison
comparison_name <- "C9_ALS_vs_Control"
comparison_groups <- c("C9-ALS", "Control")

cat("Comparison:", comparison_groups[1], "vs", comparison_groups[2], "\n")

# Fix sample names in count matrix if needed (underscores, not dots)
count_colnames <- colnames(pseudobulk_data$counts)
count_colnames <- str_replace_all(count_colnames, "\\.", "-")
colnames(pseudobulk_data$counts) <- count_colnames

# Check which clusters have sufficient replication
metadata <- pseudobulk_data$metadata
clusters_with_replicates <- metadata %>%
  group_by(cluster_name, group) %>%
  summarise(n_pseudobulks = n(), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = n_pseudobulks, values_fill = 0)

cat("\nPseudobulks per cluster per group:\n")
print(as.data.frame(clusters_with_replicates))

# Keep only clusters with ≥2 pseudobulks in each group
valid_clusters <- clusters_with_replicates %>%
  filter(if_all(where(is.numeric), ~ . >= 2)) %>%
  pull(cluster_name)

if (length(valid_clusters) == 0) {
  stop("No clusters have sufficient replication (≥2 pseudobulks per group)")
}

cat("\nClusters with sufficient replication:", length(valid_clusters), "\n")
cat("  -", paste(valid_clusters, collapse = "\n  - "), "\n")

# Filter metadata and counts
filtered_metadata <- metadata %>%
  filter(cluster_name %in% valid_clusters)

cat("\nFiltered to", nrow(filtered_metadata), "pseudobulks\n")

# Create cluster_group variable for DESeq2 design
filtered_metadata$cluster_group <- paste0(
  filtered_metadata$cluster_name, "_", filtered_metadata$group
)

# Set rownames for DESeq2
rownames(filtered_metadata) <- filtered_metadata$pseudobulk_id

################################################################################
# Filter Low Expression Genes
################################################################################

cat("\n>>> Filtering lowly expressed genes...\n")

# Calculate counts per cell for each pseudobulk
pseudobulks <- filtered_metadata$pseudobulk_id
count_matrix <- pseudobulk_data$counts[, pseudobulks]
cells_per_bulk <- lengths(pseudobulk_data$cell_list)[pseudobulks]

# Normalize by cell count
normalized_counts <- count_matrix
for (bulk_id in pseudobulks) {
  normalized_counts[, bulk_id] <- count_matrix[, bulk_id] / cells_per_bulk[bulk_id]
}

# Keep genes with >0.1 counts/cell in at least one pseudobulk
min_counts_per_cell <- 0.1
keep_genes <- rownames(normalized_counts)[apply(normalized_counts, 1, max) > min_counts_per_cell]

cat("Genes before filtering:", nrow(count_matrix), "\n")
cat("Genes after filtering (>", min_counts_per_cell, "counts/cell):", length(keep_genes), "\n")

# Filter count matrix
filtered_counts <- count_matrix[keep_genes, ]

################################################################################
# Run DESeq2 Analysis
################################################################################

cat("\n>>> Running DESeq2 analysis...\n")
cat("This may take several minutes...\n")

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = filtered_metadata,
  design = ~ cluster_group
)

# Run DESeq2 (Wald test)
dds <- DESeq(dds, quiet = FALSE)

cat("DESeq2 analysis complete\n")

# Generate variance-stabilized transformation for visualization
cat("\n>>> Generating variance-stabilized transformation...\n")
vst_matrix <- assay(vst(dds, blind = FALSE))

################################################################################
# Extract DESeq2 Results
################################################################################

cat("\n>>> Extracting differential expression results...\n")

# Set fold-change threshold
fc_threshold <- log2(1.5)  # log2(1.5) ≈ 0.58
padj_threshold <- 0.05

cat("Thresholds:\n")
cat("  Adjusted p-value:", padj_threshold, "\n")
cat("  |log2 fold-change|:", round(fc_threshold, 2), "\n")

# Extract results for each cluster
deseq_results_list <- list()
deg_lists <- list()

for (cluster_id in valid_clusters) {
  cat("\n  Extracting results for", cluster_id, "...\n")
  
  # Define contrast
  group1 <- paste0(cluster_id, "_", comparison_groups[1])
  group2 <- paste0(cluster_id, "_", comparison_groups[2])
  
  # Extract results
  res <- results(dds, contrast = c("cluster_group", group1, group2))
  res_df <- as.data.frame(res)
  
  # Store results
  comparison_id <- paste0(cluster_id, "_", comparison_name)
  deseq_results_list[[comparison_id]] <- res_df
  
  # Identify upregulated genes
  upregulated <- rownames(res_df)[
    !is.na(res_df$padj) &
    res_df$padj <= padj_threshold &
    res_df$log2FoldChange >= fc_threshold
  ]
  deg_lists[[paste0(comparison_id, "_up")]] <- upregulated
  
  # Identify downregulated genes
  downregulated <- rownames(res_df)[
    !is.na(res_df$padj) &
    res_df$padj <= padj_threshold &
    res_df$log2FoldChange <= -fc_threshold
  ]
  deg_lists[[paste0(comparison_id, "_down")]] <- downregulated
  
  cat("    Upregulated:", length(upregulated), "\n")
  cat("    Downregulated:", length(downregulated), "\n")
}

################################################################################
# Annotate Transcription Factors
################################################################################

if (length(transcription_factors) > 0) {
  cat("\n>>> Annotating transcription factors...\n")
  
  tf_deg_lists <- lapply(deg_lists, function(genes) {
    genes[genes %in% transcription_factors]
  })
  names(tf_deg_lists) <- paste0(names(deg_lists), "_TF")
  
  # Add TF lists to deg_lists
  deg_lists <- c(deg_lists, tf_deg_lists)
}

################################################################################
# Save Results
################################################################################

cat("\n>>> Saving results...\n")

# Create results object
deseq_data <- list(
  metadata = filtered_metadata,
  cell_list = pseudobulk_data$cell_list[pseudobulks],
  counts = filtered_counts,
  deseq_dataset = dds,
  vst_matrix = vst_matrix,
  deseq_results = deseq_results_list,
  deg_lists = deg_lists,
  thresholds = list(
    padj = padj_threshold,
    log2fc = fc_threshold
  )
)

# Save DESeq2 results
output_file <- file.path(output_dir, "deseq2_dataset.rda")
save(deseq_data, file = output_file)
cat("Saved DESeq2 results to:", output_file, "\n")

# Save DEG lists as CSV
deg_matrix <- matrix(
  nrow = max(lengths(deg_lists)),
  ncol = length(deg_lists)
)
colnames(deg_matrix) <- names(deg_lists)

for (i in seq_along(deg_lists)) {
  genes <- deg_lists[[i]]
  if (length(genes) > 0) {
    deg_matrix[1:length(genes), i] <- genes
  }
}
deg_matrix[is.na(deg_matrix)] <- ""

deg_file <- file.path(output_dir, "deg_lists.csv")
write_csv(as_tibble(deg_matrix), deg_file)
cat("Saved DEG lists to:", deg_file, "\n")

# Save DEG counts
deg_counts <- tibble(
  gene_set = names(deg_lists),
  n_genes = lengths(deg_lists)
)

counts_file <- file.path(output_dir, "deg_counts.csv")
write_csv(deg_counts, counts_file)
cat("Saved DEG counts to:", counts_file, "\n")

################################################################################
# Summary Statistics
################################################################################

cat("\n>>> Computing summary statistics...\n")

# Total DEGs per cluster
total_degs_per_cluster <- tibble(
  cluster = valid_clusters,
  upregulated = NA_integer_,
  downregulated = NA_integer_,
  total = NA_integer_
)

for (i in seq_along(valid_clusters)) {
  cluster_id <- valid_clusters[i]
  comparison_id <- paste0(cluster_id, "_", comparison_name)
  
  total_degs_per_cluster$upregulated[i] <- length(deg_lists[[paste0(comparison_id, "_up")]])
  total_degs_per_cluster$downregulated[i] <- length(deg_lists[[paste0(comparison_id, "_down")]])
  total_degs_per_cluster$total[i] <- total_degs_per_cluster$upregulated[i] + 
                                      total_degs_per_cluster$downregulated[i]
}

cat("\nDEGs per astrocyte cluster:\n")
print(as.data.frame(total_degs_per_cluster))

################################################################################
# Completion Summary
################################################################################

cat("\n")
cat("========================================================================\n")
cat("D10: Astrocyte DESeq2 Analysis Complete\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n")
cat("\nSummary:\n")
cat("  Clusters analyzed:", length(valid_clusters), "\n")
cat("  Genes tested:", nrow(filtered_counts), "\n")
cat("  Total DEGs:", sum(total_degs_per_cluster$total), "\n")
cat("  Total upregulated:", sum(total_degs_per_cluster$upregulated), "\n")
cat("  Total downregulated:", sum(total_degs_per_cluster$downregulated), "\n")
cat("\nThresholds:\n")
cat("  Adjusted p-value: <", padj_threshold, "\n")
cat("  |log2 fold-change|: >", round(fc_threshold, 2), "\n")
cat("\nOutput files:\n")
cat("  -", output_file, "\n")
cat("  -", deg_file, "\n")
cat("  -", counts_file, "\n")
cat("\nNext steps:\n")
cat("  1. Run D11_astrocyte_pathway_enrichment.R for functional analysis\n")
cat("  2. Perform GO, KEGG, and MSigDB enrichment on DEGs\n")
cat("  3. Create visualizations (volcano plots, heatmaps)\n")
cat("========================================================================\n\n")
