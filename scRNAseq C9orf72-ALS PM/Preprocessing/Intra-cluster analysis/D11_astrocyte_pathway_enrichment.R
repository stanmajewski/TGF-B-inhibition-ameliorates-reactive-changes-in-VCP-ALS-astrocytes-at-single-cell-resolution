#!/usr/bin/env Rscript

################################################################################
# Script: D11_astrocyte_pathway_enrichment.R
# Description: Pathway enrichment analysis and visualization of astrocyte DEGs
#
# This script:
#   1. Loads DESeq2 results from D10
#   2. Performs GO biological process enrichment on DEGs
#   3. Analyzes MSigDB pathway enrichment (Hallmark, C2, C5)
#   4. Creates enrichment heatmaps showing DEG overlap by cluster
#   5. Generates DEG expression heatmaps (Z-scores)
#   6. Plots DEG counts per cluster
#   7. Creates functional heatmaps by GO terms
#
# Inputs:
#   - results/09_astrocyte_de/deseq2_dataset.rda: DESeq2 results
#
# Outputs:
#   - results/09_astrocyte_de/go_enrichment_results.csv: GO enrichment
#   - results/09_astrocyte_de/msigdb_enrichment_results.csv: MSigDB enrichment
#   - results/09_astrocyte_de/deg_counts_barplot.pdf: DEG count visualization
#   - results/09_astrocyte_de/go_enrichment_heatmaps.pdf: GO term heatmaps
#   - results/09_astrocyte_de/msigdb_enrichment_heatmaps.pdf: Pathway heatmaps
#   - results/09_astrocyte_de/deg_expression_heatmaps.pdf: Expression heatmaps
#   - results/09_astrocyte_de/go_functional_heatmaps.pdf: GO-based expression
#
# Note: Requires clusterProfiler and msigdbr packages
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(colorRamps)
  library(viridis)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
})

# Initialize logging
cat("\n")
cat("========================================================================\n")
cat("D11: Astrocyte Pathway Enrichment and Visualization\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n\n")

# Set random seed
set.seed(123)

# Define paths
project_dir <- getwd()
input_file <- "results/09_astrocyte_de/deseq2_dataset.rda"
output_dir <- "results/09_astrocyte_de"

################################################################################
# Define Functions
################################################################################

# Color palette function
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

# Heatmap function for expression data
create_expression_heatmap <- function(matrix_data, metadata, 
                                      annotation_cols = NULL,
                                      cluster_rows = TRUE, 
                                      cluster_cols = FALSE,
                                      color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                                      limits = NULL,
                                      cellwidth = 15, 
                                      cellheight = 10,
                                      fontsize = 10, 
                                      title = "Gene Expression Z-score") {
  
  # Define annotation bars
  if (!is.null(annotation_cols)) {
    annot_df <- data.frame(row.names = colnames(matrix_data))
    
    for (col in annotation_cols) {
      values <- metadata[[col]][match(colnames(matrix_data), metadata$cluster_name)]
      values <- factor(values, levels = unique(metadata[[col]]))
      annot_df[[col]] <- values
    }
    
    # Create colors for annotations
    annot_colors <- lapply(annotation_cols, function(col) {
      vals <- generate_palette(unique(annot_df[[col]]))
      names(vals) <- levels(annot_df[[col]])
      return(vals)
    })
    names(annot_colors) <- annotation_cols
  } else {
    annot_df <- NULL
    annot_colors <- NULL
  }
  
  # Define plot limits
  if (is.null(limits)) {
    max_val <- 0.7 * max(abs(matrix_data), na.rm = TRUE)
    limits <- c(-max_val, max_val)
  }
  
  # Create plot
  pheatmap(
    matrix_data,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    color = color,
    breaks = seq(limits[1], limits[2], length.out = length(color) + 1),
    annotation_col = annot_df,
    annotation_colors = annot_colors,
    border_color = NA,
    cellwidth = cellwidth,
    cellheight = cellheight,
    fontsize = fontsize,
    main = title
  )
}

################################################################################
# Load Data
################################################################################

cat("\n>>> Loading DESeq2 results...\n")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}
load(input_file)

cat("Loaded DESeq2 data:\n")
cat("  Clusters:", length(unique(deseq_data$metadata$cluster_name)), "\n")
cat("  DEG lists:", length(deseq_data$deg_lists), "\n")

metadata <- deseq_data$metadata
cluster_names <- unique(metadata$cluster_name)
comparison_name <- "C9_ALS_vs_Control"

################################################################################
# Plot DEG Counts
################################################################################

cat("\n>>> Creating DEG count visualizations...\n")

# Prepare data
deg_summary <- tibble(
  comparison = names(deseq_data$deg_lists),
  cluster_name = NA_character_,
  direction = NA_character_,
  n_genes = lengths(deseq_data$deg_lists)
) %>%
  filter(!grepl("_TF$", comparison))  # Exclude TF-only lists

# Extract cluster names and directions
for (cluster_id in cluster_names) {
  deg_summary$cluster_name[grepl(cluster_id, deg_summary$comparison)] <- cluster_id
}
deg_summary$direction[grepl("_up$", deg_summary$comparison)] <- "Upregulated"
deg_summary$direction[grepl("_down$", deg_summary$comparison)] <- "Downregulated"

# Create bar plot
deg_barplot <- ggplot(deg_summary, aes(x = cluster_name, y = n_genes, fill = direction)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.7)) +
  scale_x_discrete(limits = cluster_names) +
  scale_fill_manual(
    values = c("Downregulated" = "blue", "Upregulated" = "red")
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "Differentially Expressed Genes in C9-ALS vs Control",
    subtitle = paste("Astrocyte clusters (padj < 0.05, |log2FC| >", 
                     round(deseq_data$thresholds$log2fc, 2), ")"),
    x = "Astrocyte Cluster",
    y = "Number of DEGs",
    fill = "Direction"
  )

# Save plot
barplot_file <- file.path(output_dir, "deg_counts_barplot.pdf")
pdf(barplot_file, width = 5, height = 4)
print(deg_barplot)
dev.off()
cat("Saved DEG count plot to:", barplot_file, "\n")

################################################################################
# GO Enrichment Analysis
################################################################################

cat("\n>>> Performing GO enrichment analysis...\n")

go_results_list <- list()

for (cluster_id in cluster_names) {
  cat("  Analyzing", cluster_id, "...\n")
  
  # Get DEGs for this cluster (up + down)
  comparison_id <- paste0(cluster_id, "_", comparison_name)
  deg_genes <- unique(c(
    deseq_data$deg_lists[[paste0(comparison_id, "_up")]],
    deseq_data$deg_lists[[paste0(comparison_id, "_down")]]
  ))
  
  if (length(deg_genes) < 5) {
    cat("    Skipping (too few DEGs)\n")
    next
  }
  
  # Perform GO enrichment
  tryCatch({
    ego <- enrichGO(
      gene = deg_genes,
      OrgDb = org.Hs.eg.db,
      keyType = 'SYMBOL',
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.01,
      qvalueCutoff = 0.05
    )
    
    if (nrow(ego@result) > 0) {
      go_results <- ego@result[ego@result$p.adjust <= 0.05, ]
      go_results_list[[cluster_id]] <- go_results
      cat("    Found", nrow(go_results), "enriched GO terms\n")
    } else {
      cat("    No significant GO terms\n")
    }
  }, error = function(e) {
    cat("    Error in GO enrichment:", e$message, "\n")
  })
}

# Save GO results
if (length(go_results_list) > 0) {
  go_file <- file.path(output_dir, "go_enrichment_results.csv")
  
  combined_go <- bind_rows(
    lapply(names(go_results_list), function(cluster_id) {
      df <- go_results_list[[cluster_id]]
      df$cluster <- cluster_id
      return(df)
    })
  )
  
  write_csv(combined_go, go_file)
  cat("\nSaved GO results to:", go_file, "\n")
}

################################################################################
# MSigDB Enrichment Analysis
################################################################################

cat("\n>>> Performing MSigDB enrichment analysis...\n")

# Get MSigDB gene sets
msigdb_categories <- c("H", "C2", "C5")  # Hallmark, Curated, GO
msigdb_results_list <- list()

for (cluster_id in cluster_names) {
  cat("  Analyzing", cluster_id, "...\n")
  
  comparison_id <- paste0(cluster_id, "_", comparison_name)
  deg_genes <- unique(c(
    deseq_data$deg_lists[[paste0(comparison_id, "_up")]],
    deseq_data$deg_lists[[paste0(comparison_id, "_down")]]
  ))
  
  if (length(deg_genes) < 5) {
    cat("    Skipping (too few DEGs)\n")
    next
  }
  
  # Enrichment for each MSigDB category
  for (category in msigdb_categories) {
    tryCatch({
      msigdb_sets <- msigdbr(species = "Homo sapiens", category = category)
      
      em <- enricher(
        gene = deg_genes,
        TERM2GENE = msigdb_sets[, c("gs_name", "gene_symbol")],
        pAdjustMethod = "BH",
        pvalueCutoff = 0.01,
        qvalueCutoff = 0.05
      )
      
      if (nrow(em@result) > 0) {
        results <- em@result[em@result$p.adjust <= 0.05, ]
        results$cluster <- cluster_id
        results$category <- category
        msigdb_results_list[[paste0(cluster_id, "_", category)]] <- results
        cat("    ", category, ":", nrow(results), "enriched sets\n")
      }
    }, error = function(e) {
      cat("    Error in", category, "enrichment:", e$message, "\n")
    })
  }
}

# Save MSigDB results
if (length(msigdb_results_list) > 0) {
  msigdb_file <- file.path(output_dir, "msigdb_enrichment_results.csv")
  combined_msigdb <- bind_rows(msigdb_results_list)
  write_csv(combined_msigdb, msigdb_file)
  cat("\nSaved MSigDB results to:", msigdb_file, "\n")
}

################################################################################
# GO Enrichment Heatmaps
################################################################################

if (length(go_results_list) > 0) {
  cat("\n>>> Creating GO enrichment heatmaps...\n")
  
  go_heatmap_file <- file.path(output_dir, "go_enrichment_heatmaps.pdf")
  pdf(go_heatmap_file, width = 18, height = 14)
  
  # Create overlap matrix for each cluster's GO terms
  for (cluster_id in names(go_results_list)) {
    go_results <- go_results_list[[cluster_id]]
    
    if (nrow(go_results) == 0) next
    
    # Get top 20 GO terms
    if (nrow(go_results) > 20) {
      go_results <- go_results[1:20, ]
    }
    
    # Extract genes for each GO term
    go_genes <- str_split(go_results$geneID, "/")
    names(go_genes) <- go_results$ID
    
    # Get DEG lists for this cluster
    comparison_id <- paste0(cluster_id, "_", comparison_name)
    up_genes <- deseq_data$deg_lists[[paste0(comparison_id, "_up")]]
    down_genes <- deseq_data$deg_lists[[paste0(comparison_id, "_down")]]
    
    # Create matrix
    overlap_matrix <- matrix(
      0,
      nrow = nrow(go_results),
      ncol = 2,
      dimnames = list(
        go_results$Description,
        c("Upregulated", "Downregulated")
      )
    )
    
    for (i in seq_len(nrow(go_results))) {
      term_genes <- go_genes[[go_results$ID[i]]]
      overlap_matrix[i, "Upregulated"] <- length(intersect(term_genes, up_genes))
      overlap_matrix[i, "Downregulated"] <- length(intersect(term_genes, down_genes))
    }
    
    # Plot heatmap
    if (nrow(overlap_matrix) > 1) {
      pheatmap(
        overlap_matrix,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "blue"))(250),
        border_color = NA,
        cellwidth = 15,
        cellheight = 10,
        main = paste0(cluster_id, " - GO Terms vs DEGs")
      )
    }
  }
  
  dev.off()
  cat("Saved GO heatmaps to:", go_heatmap_file, "\n")
}

################################################################################
# Calculate Expression Z-scores
################################################################################

cat("\n>>> Calculating expression Z-scores...\n")

# Z-score by pseudobulk (cluster_sample)
vst_matrix <- deseq_data$vst_matrix
zscore_by_sample <- t(apply(vst_matrix, 1, scale))
colnames(zscore_by_sample) <- colnames(vst_matrix)

# Mean Z-score by cluster for control samples
control_samples <- metadata$cluster_sample[metadata$group == "Control"]
zscore_by_cluster_control <- matrix(
  nrow = nrow(vst_matrix),
  ncol = length(cluster_names),
  dimnames = list(rownames(vst_matrix), cluster_names)
)

for (cluster_id in cluster_names) {
  cluster_control_samples <- metadata$cluster_sample[
    metadata$cluster_name == cluster_id & metadata$group == "Control"
  ]
  
  if (length(cluster_control_samples) > 0) {
    zscore_by_cluster_control[, cluster_id] <- rowMeans(
      zscore_by_sample[, cluster_control_samples, drop = FALSE]
    )
  }
}

# Calculate delta Z-scores (C9-ALS - Control)
zscore_delta <- matrix(
  nrow = nrow(vst_matrix),
  ncol = length(cluster_names),
  dimnames = list(rownames(vst_matrix), cluster_names)
)

for (cluster_id in cluster_names) {
  c9_samples <- metadata$cluster_sample[
    metadata$cluster_name == cluster_id & metadata$group == "C9-ALS"
  ]
  control_samples <- metadata$cluster_sample[
    metadata$cluster_name == cluster_id & metadata$group == "Control"
  ]
  
  if (length(c9_samples) > 0 && length(control_samples) > 0) {
    c9_zscore <- rowMeans(zscore_by_sample[, c9_samples, drop = FALSE])
    control_zscore <- rowMeans(zscore_by_sample[, control_samples, drop = FALSE])
    zscore_delta[, cluster_id] <- c9_zscore - control_zscore
  }
}

################################################################################
# DEG Expression Heatmaps
################################################################################

cat("\n>>> Creating DEG expression heatmaps...\n")

deg_heatmap_file <- file.path(output_dir, "deg_expression_heatmaps.pdf")
pdf(deg_heatmap_file, width = 25, height = 60)

# Get all DEGs
comparison_id <- paste0(cluster_names[1], "_", comparison_name)
all_degs <- unique(unlist(deseq_data$deg_lists[grepl(comparison_name, names(deseq_data$deg_lists))]))

if (length(all_degs) > 1) {
  # Delta Z-scores for all DEGs
  max_val <- 0.5 * max(abs(zscore_delta[all_degs, ]), na.rm = TRUE)
  
  pheatmap(
    zscore_delta[all_degs, ],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("magenta", "black", "yellow"))(250),
    breaks = seq(-max_val, max_val, length.out = 251),
    border_color = NA,
    cellwidth = 12,
    cellheight = 10,
    fontsize = 8,
    main = "All DEGs - Expression Change (C9-ALS vs Control, Z-score)"
  )
  
  # Control expression for context
  pheatmap(
    zscore_by_cluster_control[all_degs, ],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = viridis(250),
    border_color = NA,
    cellwidth = 12,
    cellheight = 10,
    fontsize = 8,
    main = "All DEGs - Baseline Expression in Control (Z-score)"
  )
}

dev.off()
cat("Saved DEG expression heatmaps to:", deg_heatmap_file, "\n")

################################################################################
# GO Functional Heatmaps
################################################################################

if (length(go_results_list) > 0) {
  cat("\n>>> Creating GO functional heatmaps...\n")
  
  go_func_file <- file.path(output_dir, "go_functional_heatmaps.pdf")
  pdf(go_func_file, width = 25, height = 60)
  
  for (cluster_id in names(go_results_list)) {
    go_results <- go_results_list[[cluster_id]]
    
    if (nrow(go_results) == 0) next
    if (nrow(go_results) > 20) go_results <- go_results[1:20, ]
    
    go_genes <- str_split(go_results$geneID, "/")
    names(go_genes) <- go_results$ID
    
    for (go_id in names(go_genes)) {
      go_term_genes <- go_genes[[go_id]]
      
      if (length(go_term_genes) < 2) next
      
      max_val <- 0.5 * max(abs(zscore_delta[go_term_genes, ]), na.rm = TRUE)
      
      pheatmap(
        zscore_delta[go_term_genes, ],
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("magenta", "black", "yellow"))(250),
        breaks = seq(-max_val, max_val, length.out = 251),
        border_color = NA,
        cellwidth = 12,
        cellheight = 10,
        fontsize = 10,
        main = paste0(cluster_id, " - ", go_results$Description[go_results$ID == go_id],
                     "\nExpression Change (C9-ALS vs Control)")
      )
    }
  }
  
  dev.off()
  cat("Saved GO functional heatmaps to:", go_func_file, "\n")
}

################################################################################
# Session Info
################################################################################

session_file <- file.path(output_dir, "session_info_d11.txt")
writeLines(capture.output(sessionInfo()), session_file)

################################################################################
# Completion Summary
################################################################################

cat("\n")
cat("========================================================================\n")
cat("D11: Astrocyte Pathway Enrichment Complete\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n")
cat("\nSummary:\n")
cat("  Clusters analyzed:", length(cluster_names), "\n")
cat("  GO enrichment performed:", length(go_results_list) > 0, "\n")
cat("  MSigDB enrichment performed:", length(msigdb_results_list) > 0, "\n")
cat("  Total DEGs:", length(all_degs), "\n")
cat("\nOutput files created:\n")
if (file.exists(barplot_file)) cat("  -", barplot_file, "\n")
if (length(go_results_list) > 0) cat("  - GO enrichment CSV and heatmaps\n")
if (length(msigdb_results_list) > 0) cat("  - MSigDB enrichment CSV and heatmaps\n")
cat("  - DEG expression heatmaps\n")
cat("  - GO functional heatmaps\n")
cat("========================================================================\n\n")
