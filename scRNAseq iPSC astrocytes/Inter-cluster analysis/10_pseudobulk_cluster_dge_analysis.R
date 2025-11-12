#!/usr/bin/env Rscript

# ============================================================================
# Pseudobulk Cluster Differential Gene Expression Analysis
# ============================================================================
# Description: Performs pseudobulk aggregation and DESeq2 analysis comparing
#              each cluster against all others, followed by pathway enrichment
# 
# Author: Analysis Pipeline
# Last Updated: 2025-01-09
# ============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript pseudobulk_cluster_dge_analysis.R <main_directory> <metadata_file> <seurat_object>
  
  Arguments:
    main_directory : Path to main analysis directory for outputs
    metadata_file  : Path to sample metadata CSV file (e.g., results/02_filtered_metadata.csv)
    seurat_object  : Path to integrated Seurat object RDA file (e.g., results/04_integrated_seurat.rda)
  
  Example:
    Rscript pseudobulk_cluster_dge_analysis.R \\
      analysis/ \\
      results/02_filtered_metadata.csv \\
      results/04_integrated_seurat.rda
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
cat("Pseudobulk Cluster Differential Gene Expression Analysis\n")
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
  library(tidyverse)
  library(Signac)
  library(Seurat)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(colorRamps)
  library(viridis)
  library(lmerTest)
  library(pheatmap)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(DOSE)
  library(ReactomePA)
  library(ggplot2)
  library(ggrepel)
  library(rrvgo)
})
cat("All packages loaded successfully.\n\n")

# Create output directory structure
cat("Setting up output directories...\n")
out_dir <- paste0(main_dir, "/pseudobulk_cluster_dge/")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Create subdirectories
subdirs <- c("deseq2_results", "filtered_results", "go_enrichment", 
             "kegg_pathways", "reactome_pathways")
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

# Set identity to clustering resolution
Idents(seur) <- "SCT_snn_res.0.3"
cat("Using clustering resolution: SCT_snn_res.0.3\n")
cat("Number of clusters:", length(unique(Idents(seur))), "\n\n")

# ============================================================================
# PSEUDOBULK AGGREGATION
# ============================================================================

cat("============================================================================\n")
cat("Performing pseudobulk aggregation...\n")
cat("============================================================================\n")

pseudo_seur <- AggregateExpression(
  seur,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = c("sample", "SCT_snn_res.0.3")
)

cat("Pseudobulk aggregation complete\n")
cat("Number of pseudobulk samples:", ncol(pseudo_seur), "\n\n")

# Extract counts matrix
counts <- GetAssayData(pseudo_seur, layer = "counts", assay = "RNA")

# Create metadata table for DESeq2
metadata <- data.frame(row.names = colnames(counts))
metadata$sample <- gsub("-", "_", metadata$sample)
metadata$sample <- gsub("^(.*?)_.*$", "\\1", rownames(metadata))
metadata$cluster <- gsub("^.*?_(.*?)$", "\\1", rownames(metadata))
metadata$treatment <- gr_tab$treatment[match(metadata$sample, gr_tab$sample)]

cat("Metadata prepared for", nrow(metadata), "pseudobulk samples\n\n")

# ============================================================================
# DESEQ2 DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================

cat("============================================================================\n")
cat("Running DESeq2 analysis for each cluster vs all others...\n")
cat("============================================================================\n")

clusters <- unique(metadata$cluster)
results_list <- list()

for (cluster in clusters) {
  cat("Processing cluster", cluster, "...\n")
  
  # Create comparison column
  metadata$is_target_cluster <- ifelse(metadata$cluster == cluster, "target", "other")
  metadata$is_target_cluster <- factor(metadata$is_target_cluster, levels = c("other", "target"))
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = metadata,
    design = ~ is_target_cluster
  )
  
  # Filter low count genes
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds, contrast = c("is_target_cluster", "target", "other"))
  
  # Store results
  results_list[[cluster]] <- res
  
  # Write results to file
  res_df <- as.data.frame(res)
  write.csv(res_df, paste0(out_dir, "deseq2_results/cluster_", cluster, "_results.csv"))
  
  # Create MA plot
  pdf(paste0(out_dir, "deseq2_results/cluster_", cluster, "_MAplot.pdf"))
  plotMA(res, main = paste("Cluster", cluster, "vs Others"))
  dev.off()
  
  # Create volcano plot
  pdf(paste0(out_dir, "deseq2_results/cluster_", cluster, "_volcano.pdf"))
  with(res_df, plot(log2FoldChange, -log10(pvalue), 
                    pch = 20, main = paste("Cluster", cluster, "vs Others"),
                    xlab = "log2 fold change", ylab = "-log10 p-value"))
  with(subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1), 
       points(log2FoldChange, -log10(pvalue), pch = 20, col = "red"))
  dev.off()
}

# Save all results
save(results_list, file = paste0(out_dir, "deseq2_results/all_deseq2_results.rda"))
cat("\nDESeq2 analysis complete for", length(results_list), "clusters\n\n")

# Generate summary table
summary_df <- data.frame(
  Cluster = character(),
  Total_DE_Genes = integer(),
  Upregulated = integer(),
  Downregulated = integer(),
  stringsAsFactors = FALSE
)

for (cluster in names(results_list)) {
  res <- results_list[[cluster]]
  res_df <- as.data.frame(res)
  de_genes <- sum(res_df$padj < 0.05 & !is.na(res_df$padj))
  up_genes <- sum(res_df$padj < 0.05 & !is.na(res_df$padj) & res_df$log2FoldChange > 0)
  down_genes <- sum(res_df$padj < 0.05 & !is.na(res_df$padj) & res_df$log2FoldChange < 0)
  
  summary_df <- rbind(summary_df, data.frame(
    Cluster = cluster,
    Total_DE_Genes = de_genes,
    Upregulated = up_genes,
    Downregulated = down_genes
  ))
}

write.csv(summary_df, paste0(out_dir, "deseq2_results/de_summary_all_clusters.csv"), row.names = FALSE)
cat("Summary:\n")
print(summary_df)
cat("\n")

# ============================================================================
# FILTERING BY FOLD CHANGE
# ============================================================================

cat("============================================================================\n")
cat("Filtering results by fold change threshold (FC > 2)...\n")
cat("============================================================================\n")

filtered_results_list <- list()
log2FC_threshold <- log2(2)

for (cluster in names(results_list)) {
  res_df <- as.data.frame(results_list[[cluster]])
  res_df$gene <- rownames(res_df)
  
  # Filter
  filtered_res <- res_df[!is.na(res_df$padj) & 
                         res_df$padj < 0.05 & 
                         abs(res_df$log2FoldChange) > log2FC_threshold, ]
  
  filtered_results_list[[cluster]] <- filtered_res
  
  # Separate up and down regulated genes
  up_genes <- filtered_res[filtered_res$log2FoldChange > 0, ]
  down_genes <- filtered_res[filtered_res$log2FoldChange < 0, ]
  
  # Get top genes
  top_up_fc <- up_genes[order(-up_genes$log2FoldChange)[1:min(10, nrow(up_genes))], ]
  top_down_fc <- down_genes[order(down_genes$log2FoldChange)[1:min(10, nrow(down_genes))], ]
  top_up_pval <- up_genes[order(up_genes$padj)[1:min(10, nrow(up_genes))], ]
  top_down_pval <- down_genes[order(down_genes$padj)[1:min(10, nrow(down_genes))], ]
  
  top_genes_rows <- unique(c(rownames(top_up_fc), rownames(top_down_fc),
                             rownames(top_up_pval), rownames(top_down_pval)))
  
  # Create regulation status column
  res_df$regulation <- "Not Significant"
  res_df$regulation[res_df$log2FoldChange > log2FC_threshold & 
                    res_df$padj < 0.05 & !is.na(res_df$padj)] <- "Up-regulated"
  res_df$regulation[res_df$log2FoldChange < -log2FC_threshold & 
                    res_df$padj < 0.05 & !is.na(res_df$padj)] <- "Down-regulated"
  
  res_df$to_label <- rownames(res_df) %in% top_genes_rows
  
  # Create volcano plot with labels
  plot_title <- paste("Cluster", as.numeric(cluster) + 1)
  max_y <- max(-log10(res_df$padj[!is.na(res_df$padj) & is.finite(-log10(res_df$padj))]))
  y_upper <- ceiling(max_y)
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = regulation), size = 1.5) +
    scale_color_manual(values = c("Not Significant" = "gray80", 
                                  "Up-regulated" = "red", 
                                  "Down-regulated" = "blue")) +
    geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), 
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), 
               linetype = "dashed", color = "black") +
    geom_text_repel(
      data = subset(res_df, to_label == TRUE),
      aes(label = gene),
      color = "black",
      size = 5.5,
      box.padding = 0.5,
      point.padding = 0.3,
      force = 10,
      segment.size = 0.2,
      segment.color = "grey50",
      max.overlaps = 30
    ) +
    labs(
      title = plot_title,
      x = expression(log[2]~fold~change),
      y = expression(-log[10]~adjusted~p-value),
      color = "Gene regulation"
    ) +
    scale_y_continuous(
      breaks = seq(0, y_upper, by = 10),
      limits = c(0, y_upper),
      expand = c(0, 0)
    ) +
    theme_bw(base_size = 20) +
    theme(
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(
    filename = paste0(out_dir, "filtered_results/cluster_", cluster, "_volcano_labeled.pdf"),
    plot = p,
    width = 12,
    height = 9,
    dpi = 300
  )
  
  # Write files
  write.csv(filtered_res, paste0(out_dir, "filtered_results/cluster_", cluster, "_filtered_fc2.csv"))
  write.csv(up_genes, paste0(out_dir, "filtered_results/cluster_", cluster, "_upregulated_fc2.csv"))
  write.csv(down_genes, paste0(out_dir, "filtered_results/cluster_", cluster, "_downregulated_fc2.csv"))
}

save(filtered_results_list, file = paste0(out_dir, "filtered_results/filtered_results.rda"))
cat("Filtering complete for", length(filtered_results_list), "clusters\n\n")

# Generate filtered summary
filtered_summary_df <- data.frame(
  Cluster = character(),
  Total_DE_Genes = integer(),
  DE_Genes_FC2 = integer(),
  Upregulated_FC2 = integer(),
  Downregulated_FC2 = integer(),
  stringsAsFactors = FALSE
)

for (cluster in names(results_list)) {
  res_df <- as.data.frame(results_list[[cluster]])
  filtered_res <- filtered_results_list[[cluster]]
  
  de_genes <- sum(res_df$padj < 0.05 & !is.na(res_df$padj))
  fc_genes <- nrow(filtered_res)
  up_genes <- sum(filtered_res$log2FoldChange > 0)
  down_genes <- sum(filtered_res$log2FoldChange < 0)
  
  filtered_summary_df <- rbind(filtered_summary_df, data.frame(
    Cluster = cluster,
    Total_DE_Genes = de_genes,
    DE_Genes_FC2 = fc_genes,
    Upregulated_FC2 = up_genes,
    Downregulated_FC2 = down_genes
  ))
}

write.csv(filtered_summary_df, paste0(out_dir, "filtered_results/filtered_summary.csv"), row.names = FALSE)
cat("Filtered summary:\n")
print(filtered_summary_df)
cat("\n")

# ============================================================================
# GO ENRICHMENT ANALYSIS
# ============================================================================

cat("============================================================================\n")
cat("Running GO enrichment analysis...\n")
cat("============================================================================\n")

# Function to run GO analysis
run_go_analysis <- function(gene_list, cluster, direction, output_dir) {
  if (length(gene_list) < 5) {
    message("Too few genes (", length(gene_list), ") for GO analysis in cluster ", 
            cluster, " ", direction)
    return(NULL)
  }
  
  message("Running GO analysis for ", length(gene_list), " ", direction, " genes in cluster ", cluster)
  
  cluster_display <- as.numeric(cluster) + 1
  
  tryCatch({
    # Biological Process
    go_bp <- enrichGO(
      gene = gene_list,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    
    if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
      go_bp <- simplify(go_bp, cutoff = 0.7, by = "p.adjust", select_fun = min)
    }
    
    # Molecular Function
    go_mf <- enrichGO(
      gene = gene_list,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    
    if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
      go_mf <- simplify(go_mf, cutoff = 0.7, by = "p.adjust", select_fun = min)
    }
    
    # Cellular Component
    go_cc <- enrichGO(
      gene = gene_list,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    
    if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
      go_cc <- simplify(go_cc, cutoff = 0.7, by = "p.adjust", select_fun = min)
    }
    
    # Save results
    if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
      write.csv(go_bp@result, paste0(output_dir, "GO_BP_", direction, "_cluster_", cluster, ".csv"))
      
      if (nrow(go_bp@result) >= 5) {
        pdf(paste0(output_dir, "GO_BP_dotplot_", direction, "_cluster_", cluster, ".pdf"), width=12, height=10)
        print(dotplot(go_bp, showCategory=15, title=paste("GO BP - Cluster", cluster_display, direction)) +
                theme(text = element_text(size = 18),
                      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
                      axis.text = element_text(size = 16)))
        dev.off()
        
        pdf(paste0(output_dir, "GO_BP_enrichmap_", direction, "_cluster_", cluster, ".pdf"), width=14, height=12)
        print(emapplot(pairwise_termsim(go_bp), showCategory = 15) +
                theme(text = element_text(size = 18),
                      legend.text = element_text(size = 16)))
        dev.off()
      }
    }
    
    if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
      write.csv(go_mf@result, paste0(output_dir, "GO_MF_", direction, "_cluster_", cluster, ".csv"))
      
      if (nrow(go_mf@result) >= 5) {
        pdf(paste0(output_dir, "GO_MF_dotplot_", direction, "_cluster_", cluster, ".pdf"), width=12, height=10)
        print(dotplot(go_mf, showCategory=15, title=paste("GO MF - Cluster", cluster_display, direction)) +
                theme(text = element_text(size = 18),
                      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
                      axis.text = element_text(size = 16)))
        dev.off()
      }
    }
    
    if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
      write.csv(go_cc@result, paste0(output_dir, "GO_CC_", direction, "_cluster_", cluster, ".csv"))
      
      if (nrow(go_cc@result) >= 5) {
        pdf(paste0(output_dir, "GO_CC_dotplot_", direction, "_cluster_", cluster, ".pdf"), width=12, height=10)
        print(dotplot(go_cc, showCategory=15, title=paste("GO CC - Cluster", cluster_display, direction)) +
                theme(text = element_text(size = 18),
                      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
                      axis.text = element_text(size = 16)))
        dev.off()
      }
    }
    
    return(list(BP = go_bp, MF = go_mf, CC = go_cc))
    
  }, error = function(e) {
    message("Error in GO analysis for cluster ", cluster, " ", direction, ": ", e$message)
    return(NULL)
  })
}

# Run GO analysis for each cluster
go_results <- list()

for (cluster in names(filtered_results_list)) {
  message("Processing GO analysis for cluster ", cluster)
  
  filtered_genes <- filtered_results_list[[cluster]]
  
  if (nrow(filtered_genes) == 0) {
    message("No significant genes for cluster ", cluster)
    next
  }
  
  up_genes <- rownames(filtered_genes[filtered_genes$log2FoldChange > 0, ])
  down_genes <- rownames(filtered_genes[filtered_genes$log2FoldChange < 0, ])
  
  cluster_dir <- paste0(out_dir, "go_enrichment/cluster_", cluster, "/")
  if (!dir.exists(cluster_dir)) {
    dir.create(cluster_dir, recursive = TRUE)
  }
  
  go_results[[paste0(cluster, "_up")]] <- run_go_analysis(up_genes, cluster, "upregulated", cluster_dir)
  go_results[[paste0(cluster, "_down")]] <- run_go_analysis(down_genes, cluster, "downregulated", cluster_dir)
}

save(go_results, file = paste0(out_dir, "go_enrichment/all_go_results.rda"))
cat("\nGO enrichment analysis complete\n\n")

# Create summary plot
extract_top_go <- function(go_result, top_n = 10) {
  if (is.null(go_result) || is.null(go_result$BP) || 
      is.null(go_result$BP@result) || nrow(go_result$BP@result) == 0) 
    return(data.frame())
  
  go_result$BP@result %>%
    arrange(pvalue) %>%
    slice_head(n = top_n)
}

pdf(paste0(out_dir, "go_enrichment/GO_terms_by_cluster.pdf"), width = 14, height = 12)

for (cluster in unique(gsub("_.*", "", names(go_results)))) {
  cluster_display <- as.numeric(cluster) + 1
  
  top_go_up <- extract_top_go(go_results[[paste0(cluster, "_up")]], top_n = 10)
  top_go_down <- extract_top_go(go_results[[paste0(cluster, "_down")]], top_n = 10)
  
  if (nrow(top_go_up) == 0 & nrow(top_go_down) == 0) {
    message("No GO terms found for cluster ", cluster)
    next
  }
  
  if (nrow(top_go_up) > 0) top_go_up$Direction <- "Up"
  if (nrow(top_go_down) > 0) top_go_down$Direction <- "Down"
  
  top_go_combined <- bind_rows(top_go_up, top_go_down)
  
  print(
    ggplot(top_go_combined, aes(x = Direction, y = Description)) +
      geom_point(aes(size = Count, color = -log10(pvalue))) +
      scale_color_gradient(low = "blue", high = "red") +
      theme_bw(base_size = 18) +
      labs(
        x = "Direction",
        y = "GO Term",
        color = expression(-log[10]~p-value),
        size = "Gene Count",
        title = paste("Cluster", cluster_display)
      ) +
      theme(
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()
      )
  )
}

dev.off()

# ============================================================================
# KEGG PATHWAY ANALYSIS
# ============================================================================

cat("============================================================================\n")
cat("Running KEGG pathway analysis...\n")
cat("============================================================================\n")

# Function for KEGG analysis
run_kegg_analysis <- function(gene_list, cluster_name, output_dir) {
  gene_entrez <- bitr(gene_list, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Hs.eg.db")
  
  if (nrow(gene_entrez) == 0) {
    warning(paste("No genes could be converted to ENTREZ IDs for cluster", cluster_name))
    return(NULL)
  }
  
  kegg_result <- enrichKEGG(
    gene = gene_entrez$ENTREZID,
    organism = 'hsa',
    keyType = 'kegg',
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  if (is.null(kegg_result) || nrow(kegg_result) == 0) {
    message(paste("No significant KEGG pathways found for cluster", cluster_name))
    return(NULL)
  }
  
  kegg_df <- as.data.frame(kegg_result)
  write.csv(kegg_df, file = paste0(output_dir, "kegg_pathways_cluster_", cluster_name, ".csv"), 
            row.names = FALSE)
  
  kegg_plot_data <- kegg_df %>%
    mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))) %>%
    arrange(desc(GeneRatio)) %>%
    head(12) %>%
    arrange(GeneRatio) %>%
    mutate(Description = factor(Description, levels = Description))
  
  n_pathways <- nrow(kegg_plot_data)
  plot_height <- max(8, 5 + (n_pathways * 0.4))
  
  p <- ggplot(kegg_plot_data, aes(x = GeneRatio, y = Description)) +
    geom_point(aes(size = Count, color = -log10(pvalue))) +
    scale_color_gradient(low = "blue", high = "red", name = expression(-log[10]~p-value)) +
    scale_size_continuous(name = "Gene Count", range = c(3, 10)) +
    theme_bw(base_size = 16) +
    labs(
      title = paste("KEGG Pathways - Cluster", cluster_name),
      x = "Gene Ratio",
      y = NULL
    ) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      plot.margin = margin(10, 10, 10, 20)
    )
  
  ggsave(
    filename = paste0(output_dir, "kegg_pathways_cluster_", cluster_name, ".pdf"),
    plot = p,
    width = 12,
    height = plot_height,
    dpi = 300
  )
  
  return(kegg_result)
}

kegg_results_list <- list()

for (cluster in names(filtered_results_list)) {
  cluster_display <- as.character(as.numeric(cluster) + 1)
  
  up_genes <- filtered_results_list[[cluster]][filtered_results_list[[cluster]]$log2FoldChange > 0, ]
  
  if (nrow(up_genes) == 0) {
    message(paste("No upregulated genes with FC > 2 found for cluster", cluster_display))
    next
  }
  
  gene_symbols <- rownames(up_genes)
  
  message(paste("Running KEGG analysis on", length(gene_symbols), "upregulated genes for cluster", cluster_display))
  
  kegg_result <- run_kegg_analysis(gene_symbols, cluster_display, paste0(out_dir, "kegg_pathways/"))
  
  if (!is.null(kegg_result)) {
    kegg_results_list[[cluster_display]] <- kegg_result
  }
}

save(kegg_results_list, file = paste0(out_dir, "kegg_pathways/all_kegg_results.rda"))

# Create summary
kegg_summary <- data.frame(
  Cluster = character(),
  Significant_Pathways = integer(),
  Top_Pathway = character(),
  Top_Pathway_PValue = numeric(),
  stringsAsFactors = FALSE
)

for (cluster in names(kegg_results_list)) {
  if (!is.null(kegg_results_list[[cluster]])) {
    kegg_df <- as.data.frame(kegg_results_list[[cluster]])
    if (nrow(kegg_df) > 0) {
      kegg_summary <- rbind(kegg_summary, data.frame(
        Cluster = cluster,
        Significant_Pathways = nrow(kegg_df),
        Top_Pathway = kegg_df[1, "Description"],
        Top_Pathway_PValue = kegg_df[1, "pvalue"]
      ))
    }
  }
}

write.csv(kegg_summary, paste0(out_dir, "kegg_pathways/kegg_summary.csv"), row.names = FALSE)
cat("\nKEGG pathway analysis complete\n\n")

# ============================================================================
# REACTOME PATHWAY ANALYSIS
# ============================================================================

cat("============================================================================\n")
cat("Running Reactome pathway analysis...\n")
cat("============================================================================\n")

# Function for Reactome analysis
run_reactome_analysis <- function(gene_list, cluster_name, output_dir) {
  gene_entrez <- bitr(gene_list, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Hs.eg.db")
  
  if (nrow(gene_entrez) == 0) {
    warning(paste("No genes could be converted to ENTREZ IDs for cluster", cluster_name))
    return(NULL)
  }
  
  reactome_result <- enrichPathway(
    gene = gene_entrez$ENTREZID,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
  )
  
  if (is.null(reactome_result) || nrow(reactome_result) == 0) {
    message(paste("No significant Reactome pathways found for cluster", cluster_name))
    return(NULL)
  }
  
  reactome_df <- as.data.frame(reactome_result)
  write.csv(reactome_df, file = paste0(output_dir, "reactome_pathways_cluster_", cluster_name, ".csv"), 
            row.names = FALSE)
  
  reactome_plot_data <- reactome_df %>%
    mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))) %>%
    arrange(desc(GeneRatio)) %>%
    head(15) %>%
    arrange(GeneRatio) %>%
    mutate(Description = factor(Description, levels = Description))
  
  n_pathways <- nrow(reactome_plot_data)
  plot_height <- max(8, 5 + (n_pathways * 0.35))
  
  p <- ggplot(reactome_plot_data, aes(x = GeneRatio, y = Description)) +
    geom_point(aes(size = Count, color = -log10(pvalue))) +
    scale_color_gradient(low = "blue", high = "red", name = expression(-log[10]~p-value)) +
    scale_size_continuous(name = "Gene Count", range = c(3, 10)) +
    theme_bw(base_size = 16) +
    labs(
      title = paste("Reactome Pathways - Cluster", cluster_name),
      x = "Gene Ratio",
      y = NULL
    ) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 13),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      plot.margin = margin(10, 10, 10, 20)
    )
  
  ggsave(
    filename = paste0(output_dir, "reactome_pathways_cluster_", cluster_name, ".pdf"),
    plot = p,
    width = 14,
    height = plot_height,
    dpi = 300
  )
  
  return(reactome_result)
}

reactome_results_list <- list()

for (cluster in names(filtered_results_list)) {
  cluster_display <- as.character(as.numeric(cluster) + 1)
  
  up_genes <- filtered_results_list[[cluster]][filtered_results_list[[cluster]]$log2FoldChange > 0, ]
  
  if (nrow(up_genes) == 0) {
    message(paste("No upregulated genes with FC > 2 found for cluster", cluster_display))
    next
  }
  
  gene_symbols <- rownames(up_genes)
  
  message(paste("Running Reactome analysis on", length(gene_symbols), "upregulated genes for cluster", cluster_display))
  
  reactome_result <- run_reactome_analysis(gene_symbols, cluster_display, paste0(out_dir, "reactome_pathways/"))
  
  if (!is.null(reactome_result)) {
    reactome_results_list[[cluster_display]] <- reactome_result
  }
}

save(reactome_results_list, file = paste0(out_dir, "reactome_pathways/all_reactome_results.rda"))

# Create summary
reactome_summary <- data.frame(
  Cluster = character(),
  Significant_Pathways = integer(),
  Top_Pathway = character(),
  Top_Pathway_PValue = numeric(),
  Top_Pathway_Genes = character(),
  stringsAsFactors = FALSE
)

for (cluster in names(reactome_results_list)) {
  if (!is.null(reactome_results_list[[cluster]])) {
    reactome_df <- as.data.frame(reactome_results_list[[cluster]])
    if (nrow(reactome_df) > 0) {
      reactome_summary <- rbind(reactome_summary, data.frame(
        Cluster = cluster,
        Significant_Pathways = nrow(reactome_df),
        Top_Pathway = reactome_df[1, "Description"],
        Top_Pathway_PValue = reactome_df[1, "pvalue"],
        Top_Pathway_Genes = reactome_df[1, "geneID"]
      ))
    }
  }
}

write.csv(reactome_summary, paste0(out_dir, "reactome_pathways/reactome_summary.csv"), row.names = FALSE)
cat("\nReactome pathway analysis complete\n\n")

# ============================================================================
# FINAL SUMMARY AND COMPLETION
# ============================================================================

cat("\n============================================================================\n")
cat("ANALYSIS SUMMARY\n")
cat("============================================================================\n")

cat("\nDESeq2 Results:\n")
cat("  Total clusters analyzed:", length(results_list), "\n")
cat("  Total DE genes (any cluster):", sum(summary_df$Total_DE_Genes), "\n")

cat("\nFiltered Results (FC > 2):\n")
cat("  Total DE genes passing filter:", sum(filtered_summary_df$DE_Genes_FC2), "\n")

if (nrow(kegg_summary) > 0) {
  cat("\nKEGG Pathways:\n")
  cat("  Total significant pathways:", sum(kegg_summary$Significant_Pathways), "\n")
}

if (nrow(reactome_summary) > 0) {
  cat("\nReactome Pathways:\n")
  cat("  Total significant pathways:", sum(reactome_summary$Significant_Pathways), "\n")
  cat("\n  Top pathway per cluster:\n")
  for (i in 1:min(5, nrow(reactome_summary))) {
    cat(sprintf("    Cluster %s: %s (p = %.2e)\n", 
                reactome_summary$Cluster[i],
                reactome_summary$Top_Pathway[i],
                reactome_summary$Top_Pathway_PValue[i]))
  }
}

cat("\n============================================================================\n")
cat("OUTPUT FILES LOCATION\n")
cat("============================================================================\n")
cat("Main output directory:", out_dir, "\n")
cat("  - deseq2_results/     : DESeq2 results and plots\n")
cat("  - filtered_results/   : FC-filtered results and enhanced volcanos\n")
cat("  - go_enrichment/      : GO term enrichment by cluster\n")
cat("  - kegg_pathways/      : KEGG pathway enrichment\n")
cat("  - reactome_pathways/  : Reactome pathway enrichment\n")

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
