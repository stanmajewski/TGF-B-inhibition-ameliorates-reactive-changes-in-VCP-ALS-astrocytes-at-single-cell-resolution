#!/usr/bin/env Rscript

# ============================================================================
# SB Treatment Effect on VCP-Associated Gene Expression
# ============================================================================
# Description: Analyzes how SB431542 treatment affects VCP-associated 
#              differential gene expression, identifying genes that are
#              normalized, exacerbated, or unaffected by treatment
# 
# Author: Analysis Pipeline
# Last Updated: 2025-01-09
# ============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript sb_treatment_effect_vcp_analysis.R <main_directory> <metadata_file> <deseq_results_file> <tf_reference_file>
  
  Arguments:
    main_directory     : Path to main analysis directory for outputs
    metadata_file      : Path to sample metadata CSV file (e.g., B02_gr_tab_filtered.csv)
    deseq_results_file : Path to RDA file with DESeq2 results (e.g., bulk_data_with_DESeq_results.rda)
    tf_reference_file  : Path to transcription factor reference CSV (e.g., Transcription_Factors_hg19.csv)
  
  Example:
    Rscript sb_treatment_effect_vcp_analysis.R \\
      /Users/username/scRNAseq/ \\
      B02_gr_tab_filtered.csv \\
      bulk_data_with_DESeq_results.rda \\
      Transcription_Factors_hg19.csv
  ", call. = FALSE)
}

# Set environment and main directory
main_dir <- normalizePath(args[1], mustWork = TRUE)
metadata_file <- normalizePath(args[2], mustWork = TRUE)
deseq_file <- normalizePath(args[3], mustWork = TRUE)
tf_file <- normalizePath(args[4], mustWork = TRUE)

# Record start time
start_time <- Sys.time()

setwd(main_dir)

# Log start time and parameters
cat("\n============================================================================\n")
cat("SB Treatment Effect on VCP-Associated Gene Expression\n")
cat("============================================================================\n")
cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Main directory:", main_dir, "\n")
cat("Metadata file:", metadata_file, "\n")
cat("DESeq results file:", deseq_file, "\n")
cat("TF reference file:", tf_file, "\n")
cat("R version:", R.version.string, "\n")
cat("============================================================================\n\n")

# Load packages
cat("Loading required packages...\n")
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(colorRamps)
  library(viridis)
  library(pheatmap)
  library(clusterProfiler)
  library(DOSE)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(ggrepel)
  library(enrichplot)
  library(gridExtra)
})
cat("All packages loaded successfully.\n\n")

# Create output directory structure
cat("Setting up output directories...\n")
out_dir <- paste0(main_dir, "/sb_treatment_effect_vcp/")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Create subdirectories
subdirs <- c("treatment_analysis", "go_enrichment", "visualizations")
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

load(file = deseq_file)
cat("Loaded DESeq2 results object\n")

# Load marker gene panels
GOI <- list()
t1 <- read_csv(tf_file, show_col_types = FALSE)
GOI$TF <- t1$Symbol
cat("Loaded", length(GOI$TF), "transcription factors\n\n")

# Define group comparisons for each cluster
comp_groups <- list(
  VCP.u_vs_ctrl.u = c("VCP_u", "ctrl_u"),
  VCP.SB_vs_VCP.u = c("VCP_SB", "VCP_u"),
  VCP.SB_vs_ctrl.u = c("VCP_SB", "ctrl_u")
)

# ============================================================================
# DETERMINE CLUSTERS TO ANALYZE
# ============================================================================

cat("============================================================================\n")
cat("Determining available clusters...\n")
cat("============================================================================\n")

# Function to extract available clusters from deseq_results
get_available_clusters <- function(bulk_data) {
  all_keys <- names(bulk_data$deseq_results)
  
  cluster_names <- unique(sapply(all_keys, function(key) {
    parts <- strsplit(key, "_")[[1]]
    return(parts[1])
  }))
  
  return(cluster_names[!cluster_names %in% c("all")])
}

comp_clusters <- get_available_clusters(bulk_data)
cat("Found", length(comp_clusters), "clusters to analyze:", paste(comp_clusters, collapse = ", "), "\n\n")

# ============================================================================
# ANALYZE TREATMENT IMPACT ON VCP-AFFECTED GENES
# ============================================================================

cat("============================================================================\n")
cat("Analyzing treatment impact on VCP-affected genes...\n")
cat("============================================================================\n")

# Function to analyze treatment impact using Log2FC-based thresholds
analyze_treatment_impact <- function(bulk_data, clusters = comp_clusters) {
  treatment_impact <- list()
  
  for (cl in clusters) {
    cat("  * Analyzing treatment impact for cluster:", cl, "\n")
    
    # Get DEGs from VCP vs control comparison
    ctrl_vs_vcp_key <- paste0(cl, "_VCP.u_vs_ctrl.u")
    ctrl_vs_vcp_results <- bulk_data$deseq_results[[ctrl_vs_vcp_key]]
    
    if (is.null(ctrl_vs_vcp_results)) {
      cat("    - No results found for", ctrl_vs_vcp_key, "\n")
      next
    }
    
    # Get genes up and down in VCP vs control
    vcp_up_genes <- bulk_data$DEGs[[paste0(ctrl_vs_vcp_key, "_down")]]  # UP in VCP
    vcp_down_genes <- bulk_data$DEGs[[paste0(ctrl_vs_vcp_key, "_up")]]  # DOWN in VCP
    
    cat("    - Found", length(vcp_up_genes), "genes up-regulated in VCP\n")
    cat("    - Found", length(vcp_down_genes), "genes down-regulated in VCP\n")
    
    # Get results from VCP+SB vs control comparison
    ctrl_vs_vcp_sb_key <- paste0(cl, "_VCP.SB_vs_ctrl.u")
    ctrl_vs_vcp_sb_results <- bulk_data$deseq_results[[ctrl_vs_vcp_sb_key]]
    
    if (is.null(ctrl_vs_vcp_sb_results)) {
      cat("    - No results found for", ctrl_vs_vcp_sb_key, "\n")
      next
    }
    
    # Analyze VCP upregulated genes
    if (length(vcp_up_genes) > 0) {
      normalized_genes <- c()
      exacerbated_genes <- c()
      unaffected_genes <- c()
      
      for (gene in vcp_up_genes) {
        if (gene %in% rownames(ctrl_vs_vcp_sb_results)) {
          l2fc_vcp <- ctrl_vs_vcp_results[gene, "log2FoldChange"]
          l2fc_sb <- ctrl_vs_vcp_sb_results[gene, "log2FoldChange"]
          
          # Classification based on Log2FC change
          l2fc_change <- abs(l2fc_sb) - abs(l2fc_vcp)
          
          if (l2fc_change <= -0.5) {
            normalized_genes <- c(normalized_genes, gene)
          } else if (l2fc_change >= 0.5) {
            exacerbated_genes <- c(exacerbated_genes, gene)
          } else {
            unaffected_genes <- c(unaffected_genes, gene)
          }
        }
      }
      
      cat("    - VCP up-regulated genes normalized by treatment:", length(normalized_genes), "\n")
      cat("    - VCP up-regulated genes exacerbated by treatment:", length(exacerbated_genes), "\n")
      cat("    - VCP up-regulated genes unaffected by treatment:", length(unaffected_genes), "\n")
      
      treatment_impact[[paste0(cl, "_vcp_up_normalized_by_treatment")]] <- normalized_genes
      treatment_impact[[paste0(cl, "_vcp_up_exacerbated_by_treatment")]] <- exacerbated_genes
      treatment_impact[[paste0(cl, "_vcp_up_unaffected_by_treatment")]] <- unaffected_genes
    }
    
    # Analyze VCP downregulated genes
    if (length(vcp_down_genes) > 0) {
      normalized_genes <- c()
      exacerbated_genes <- c()
      unaffected_genes <- c()
      
      for (gene in vcp_down_genes) {
        if (gene %in% rownames(ctrl_vs_vcp_sb_results)) {
          l2fc_vcp <- ctrl_vs_vcp_results[gene, "log2FoldChange"]
          l2fc_sb <- ctrl_vs_vcp_sb_results[gene, "log2FoldChange"]
          
          # Classification based on Log2FC change
          l2fc_change <- abs(l2fc_sb) - abs(l2fc_vcp)
          
          if (l2fc_change <= -0.5) {
            normalized_genes <- c(normalized_genes, gene)
          } else if (l2fc_change >= 0.5) {
            exacerbated_genes <- c(exacerbated_genes, gene)
          } else {
            unaffected_genes <- c(unaffected_genes, gene)
          }
        }
      }
      
      cat("    - VCP down-regulated genes normalized by treatment:", length(normalized_genes), "\n")
      cat("    - VCP down-regulated genes exacerbated by treatment:", length(exacerbated_genes), "\n")
      cat("    - VCP down-regulated genes unaffected by treatment:", length(unaffected_genes), "\n")
      
      treatment_impact[[paste0(cl, "_vcp_down_normalized_by_treatment")]] <- normalized_genes
      treatment_impact[[paste0(cl, "_vcp_down_exacerbated_by_treatment")]] <- exacerbated_genes
      treatment_impact[[paste0(cl, "_vcp_down_unaffected_by_treatment")]] <- unaffected_genes
    }
  }
  
  return(treatment_impact)
}

# Run the analysis
treatment_impact <- analyze_treatment_impact(bulk_data)
bulk_data$treatment_impact <- treatment_impact

cat("\nTreatment impact analysis complete\n\n")

# ============================================================================
# GENERATE SUMMARY STATISTICS
# ============================================================================

cat("Generating summary statistics...\n")

treatment_summary <- data.frame(
  cluster = character(),
  vcp_up_total = integer(),
  vcp_up_normalized = integer(),
  vcp_up_exacerbated = integer(),
  vcp_up_unaffected = integer(),
  vcp_down_total = integer(),
  vcp_down_normalized = integer(),
  vcp_down_exacerbated = integer(),
  vcp_down_unaffected = integer(),
  stringsAsFactors = FALSE
)

for (cl in comp_clusters) {
  vcp_up_total <- length(bulk_data$DEGs[[paste0(cl, "_VCP.u_vs_ctrl.u_down")]])
  vcp_down_total <- length(bulk_data$DEGs[[paste0(cl, "_VCP.u_vs_ctrl.u_up")]])
  
  vcp_up_normalized <- length(treatment_impact[[paste0(cl, "_vcp_up_normalized_by_treatment")]])
  vcp_up_exacerbated <- length(treatment_impact[[paste0(cl, "_vcp_up_exacerbated_by_treatment")]])
  vcp_up_unaffected <- length(treatment_impact[[paste0(cl, "_vcp_up_unaffected_by_treatment")]])
  
  vcp_down_normalized <- length(treatment_impact[[paste0(cl, "_vcp_down_normalized_by_treatment")]])
  vcp_down_exacerbated <- length(treatment_impact[[paste0(cl, "_vcp_down_exacerbated_by_treatment")]])
  vcp_down_unaffected <- length(treatment_impact[[paste0(cl, "_vcp_down_unaffected_by_treatment")]])
  
  treatment_summary <- rbind(treatment_summary, data.frame(
    cluster = cl,
    vcp_up_total = vcp_up_total,
    vcp_up_normalized = vcp_up_normalized,
    vcp_up_exacerbated = vcp_up_exacerbated,
    vcp_up_unaffected = vcp_up_unaffected,
    vcp_down_total = vcp_down_total,
    vcp_down_normalized = vcp_down_normalized,
    vcp_down_exacerbated = vcp_down_exacerbated,
    vcp_down_unaffected = vcp_down_unaffected
  ))
}

bulk_data$treatment_summary <- treatment_summary

# Save summary
write.csv(treatment_summary, 
          paste0(out_dir, "treatment_analysis/treatment_summary.csv"), 
          row.names = FALSE)

cat("Summary statistics saved\n\n")

# ============================================================================
# VISUALIZATION 1: STACKED BAR CHART BY CLUSTER (PERCENTAGES)
# ============================================================================

cat("============================================================================\n")
cat("Creating treatment effect visualizations...\n")
cat("============================================================================\n")

# Calculate percentages
cluster_summary_long <- treatment_summary %>%
  mutate(
    total_vcp_genes = vcp_up_total + vcp_down_total,
    total_normalized = vcp_up_normalized + vcp_down_normalized,
    total_exacerbated = vcp_up_exacerbated + vcp_down_exacerbated,
    total_unaffected = vcp_up_unaffected + vcp_down_unaffected
  ) %>%
  select(cluster, total_normalized, total_unaffected, total_exacerbated, total_vcp_genes) %>%
  pivot_longer(cols = c(total_normalized, total_unaffected, total_exacerbated),
               names_to = "effect",
               values_to = "count")

cluster_summary_pct <- cluster_summary_long %>%
  group_by(cluster) %>%
  mutate(percentage = (count / total_vcp_genes) * 100) %>%
  ungroup() %>%
  mutate(effect = case_when(
    effect == "total_normalized" ~ "Normalized",
    effect == "total_unaffected" ~ "Unaffected",
    effect == "total_exacerbated" ~ "Exacerbated"
  ))

cluster_summary_pct$effect <- factor(cluster_summary_pct$effect,
                                     levels = c("Normalized", "Unaffected", "Exacerbated"))

# Create stacked bar plot
p1 <- ggplot(cluster_summary_pct, aes(x = cluster, y = percentage, fill = effect)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c(
    "Normalized" = "#0072B2",
    "Unaffected" = "#7F8C8D",
    "Exacerbated" = "#D55E00"
  )) +
  labs(
    title = "Treatment Impact on VCP-affected Genes by Cluster",
    subtitle = "Percentage of genes normalized, unaffected, or exacerbated by SB treatment",
    x = "Cluster",
    y = "Percentage of Genes (%)",
    fill = "Treatment Effect"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 16, color = "gray40"),
    legend.position = "bottom",
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22, face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    axis.title = element_text(size = 22, face = "bold"),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt")
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25))

print(p1)

ggsave(
  paste0(out_dir, "visualizations/treatment_impact_barplot.pdf"),
  plot = p1,
  width = 12,
  height = 8
)

cat("Stacked bar chart created\n")

# ============================================================================
# VISUALIZATION 2: LOG2FC SCATTER PLOTS PER CLUSTER
# ============================================================================

cat("Creating Log2FC comparison scatter plots...\n")

for (cl in treatment_summary$cluster) {
  ctrl_vs_vcp_key <- paste0(cl, "_VCP.u_vs_ctrl.u")
  ctrl_vs_vcp_sb_key <- paste0(cl, "_VCP.SB_vs_ctrl.u")
  
  vcp_results <- bulk_data$deseq_results[[ctrl_vs_vcp_key]]
  vcp_sb_results <- bulk_data$deseq_results[[ctrl_vs_vcp_sb_key]]
  
  if (!is.null(vcp_results) && !is.null(vcp_sb_results)) {
    vcp_up_genes <- bulk_data$DEGs[[paste0(ctrl_vs_vcp_key, "_down")]]
    vcp_down_genes <- bulk_data$DEGs[[paste0(ctrl_vs_vcp_key, "_up")]]
    all_deg_genes <- c(vcp_up_genes, vcp_down_genes)
    
    if (length(all_deg_genes) > 0) {
      common_genes <- intersect(all_deg_genes, rownames(vcp_sb_results))
      
      cluster_data <- data.frame(
        gene = common_genes,
        log2FC_VCP = vcp_results[common_genes, "log2FoldChange"],
        log2FC_VCP_SB = vcp_sb_results[common_genes, "log2FoldChange"]
      )
      
      cluster_data$l2fc_change <- abs(cluster_data$log2FC_VCP_SB) - abs(cluster_data$log2FC_VCP)
      cluster_data$effect <- ifelse(cluster_data$l2fc_change <= -0.5, "Normalized",
                                    ifelse(cluster_data$l2fc_change >= 0.5, "Exacerbated",
                                           "Unaffected"))
      
      cluster_data$effect <- factor(cluster_data$effect,
                                    levels = c("Normalized", "Unaffected", "Exacerbated"))
      
      # Calculate axis limits
      quantile_95 <- quantile(abs(c(cluster_data$log2FC_VCP, cluster_data$log2FC_VCP_SB)),
                              probs = 0.95, na.rm = TRUE)
      axis_limit <- ceiling(quantile_95)
      
      # Create scatter plot
      p_cluster <- ggplot(cluster_data, aes(x = log2FC_VCP, y = log2FC_VCP_SB, color = effect)) +
        geom_point(alpha = 0.6, size = 3, shape = 16) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 1) +
        geom_hline(yintercept = 0, linetype = "solid", color = "gray60", linewidth = 0.5) +
        geom_vline(xintercept = 0, linetype = "solid", color = "gray60", linewidth = 0.5) +
        scale_color_manual(values = c(
          "Normalized" = "#0072B2",
          "Unaffected" = "#7F8C8D",
          "Exacerbated" = "#D55E00"
        )) +
        labs(
          title = paste0("Cluster ", cl),
          x = "Log2FC (VCP_U vs CTRL)",
          y = "Log2FC (VCP_SB vs CTRL)",
          color = "Treatment Effect"
        ) +
        theme_minimal(base_size = 22) +
        theme(
          plot.title = element_text(face = "bold", size = 26),
          legend.position = "bottom",
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 22, face = "bold"),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 22),
          axis.title = element_text(size = 22, face = "bold")
        ) +
        coord_fixed(ratio = 1, xlim = c(-axis_limit, axis_limit), ylim = c(-axis_limit, axis_limit))
      
      print(p_cluster)
      
      ggsave(
        paste0(out_dir, "visualizations/cluster_", cl, "_log2fc_comparison.pdf"),
        plot = p_cluster,
        width = 10,
        height = 9
      )
    }
  }
}

cat("Log2FC scatter plots created\n\n")

# ============================================================================
# GO ENRICHMENT ANALYSIS ON NORMALIZED GENES
# ============================================================================

cat("============================================================================\n")
cat("Running GO enrichment analysis on normalized genes...\n")
cat("============================================================================\n")

# Function to convert gene symbols to Entrez IDs
convert_to_entrez <- function(gene_symbols) {
  entrez_ids <- bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  return(entrez_ids$ENTREZID)
}

# Initialize storage for GO results
go_results <- list()

# Analyze each cluster
for (cl in comp_clusters) {
  cat("  * Running GO analysis for cluster", cl, "\n")
  
  # Get normalized genes
  vcp_up_normalized <- treatment_impact[[paste0(cl, "_vcp_up_normalized_by_treatment")]]
  vcp_down_normalized <- treatment_impact[[paste0(cl, "_vcp_down_normalized_by_treatment")]]
  
  # Analyze UP-regulated genes normalized by treatment
  if (length(vcp_up_normalized) >= 5) {
    tryCatch({
      entrez_up <- convert_to_entrez(vcp_up_normalized)
      
      if (length(entrez_up) >= 5) {
        # GO Biological Process
        go_bp_up <- enrichGO(
          gene = entrez_up,
          OrgDb = org.Hs.eg.db,
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2,
          readable = TRUE
        )
        
        # GO Molecular Function
        go_mf_up <- enrichGO(
          gene = entrez_up,
          OrgDb = org.Hs.eg.db,
          ont = "MF",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2,
          readable = TRUE
        )
        
        # GO Cellular Component
        go_cc_up <- enrichGO(
          gene = entrez_up,
          OrgDb = org.Hs.eg.db,
          ont = "CC",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2,
          readable = TRUE
        )
        
        go_results[[paste0(cl, "_up_normalized_BP")]] <- go_bp_up
        go_results[[paste0(cl, "_up_normalized_MF")]] <- go_mf_up
        go_results[[paste0(cl, "_up_normalized_CC")]] <- go_cc_up
        
        cat("    - Completed GO analysis for UP-regulated normalized genes\n")
      }
    }, error = function(e) {
      cat("    - Error in UP-regulated analysis:", e$message, "\n")
    })
  }
  
  # Analyze DOWN-regulated genes normalized by treatment
  if (length(vcp_down_normalized) >= 5) {
    tryCatch({
      entrez_down <- convert_to_entrez(vcp_down_normalized)
      
      if (length(entrez_down) >= 5) {
        # GO Biological Process
        go_bp_down <- enrichGO(
          gene = entrez_down,
          OrgDb = org.Hs.eg.db,
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2,
          readable = TRUE
        )
        
        # GO Molecular Function
        go_mf_down <- enrichGO(
          gene = entrez_down,
          OrgDb = org.Hs.eg.db,
          ont = "MF",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2,
          readable = TRUE
        )
        
        # GO Cellular Component
        go_cc_down <- enrichGO(
          gene = entrez_down,
          OrgDb = org.Hs.eg.db,
          ont = "CC",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2,
          readable = TRUE
        )
        
        go_results[[paste0(cl, "_down_normalized_BP")]] <- go_bp_down
        go_results[[paste0(cl, "_down_normalized_MF")]] <- go_mf_down
        go_results[[paste0(cl, "_down_normalized_CC")]] <- go_cc_down
        
        cat("    - Completed GO analysis for DOWN-regulated normalized genes\n")
      }
    }, error = function(e) {
      cat("    - Error in DOWN-regulated analysis:", e$message, "\n")
    })
  }
}

# Store results in bulk_data
bulk_data$go_normalized_genes <- go_results

cat("\nGO enrichment analysis complete\n\n")

# ============================================================================
# GO VISUALIZATION 1: DOT PLOTS FOR EACH CLUSTER
# ============================================================================

cat("Creating GO dot plots...\n")

for (cl in comp_clusters) {
  # UP-regulated normalized genes - Biological Process
  go_key_up_bp <- paste0(cl, "_up_normalized_BP")
  if (!is.null(go_results[[go_key_up_bp]])) {
    if (nrow(go_results[[go_key_up_bp]]) > 0) {
      p <- dotplot(
        go_results[[go_key_up_bp]],
        showCategory = 15,
        title = paste0("Cluster ", cl, ": GO-BP\nNormalized UP-regulated Genes")
      ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(face = "bold", size = 16),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12)
        )
      
      print(p)
      
      ggsave(
        paste0(out_dir, "go_enrichment/GO_BP_dotplot_cluster_", cl, "_up_normalized.pdf"),
        plot = p,
        width = 10,
        height = 8
      )
    }
  }
  
  # DOWN-regulated normalized genes - Biological Process
  go_key_down_bp <- paste0(cl, "_down_normalized_BP")
  if (!is.null(go_results[[go_key_down_bp]])) {
    if (nrow(go_results[[go_key_down_bp]]) > 0) {
      p <- dotplot(
        go_results[[go_key_down_bp]],
        showCategory = 15,
        title = paste0("Cluster ", cl, ": GO-BP\nNormalized DOWN-regulated Genes")
      ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(face = "bold", size = 16),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12)
        )
      
      print(p)
      
      ggsave(
        paste0(out_dir, "go_enrichment/GO_BP_dotplot_cluster_", cl, "_down_normalized.pdf"),
        plot = p,
        width = 10,
        height = 8
      )
    }
  }
}

cat("GO dot plots created\n")

# ============================================================================
# GO VISUALIZATION 2: BAR PLOTS
# ============================================================================

cat("Creating GO bar plots...\n")

for (cl in comp_clusters) {
  # UP-regulated normalized genes
  go_key_up_bp <- paste0(cl, "_up_normalized_BP")
  if (!is.null(go_results[[go_key_up_bp]])) {
    if (nrow(go_results[[go_key_up_bp]]) > 0) {
      p <- barplot(
        go_results[[go_key_up_bp]],
        showCategory = 15,
        title = paste0("Cluster ", cl, ": GO-BP\nNormalized UP-regulated Genes")
      ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(face = "bold", size = 16),
          axis.text.y = element_text(size = 11)
        )
      
      print(p)
      
      ggsave(
        paste0(out_dir, "go_enrichment/GO_BP_barplot_cluster_", cl, "_up_normalized.pdf"),
        plot = p,
        width = 10,
        height = 8
      )
    }
  }
  
  # DOWN-regulated normalized genes
  go_key_down_bp <- paste0(cl, "_down_normalized_BP")
  if (!is.null(go_results[[go_key_down_bp]])) {
    if (nrow(go_results[[go_key_down_bp]]) > 0) {
      p <- barplot(
        go_results[[go_key_down_bp]],
        showCategory = 15,
        title = paste0("Cluster ", cl, ": GO-BP\nNormalized DOWN-regulated Genes")
      ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(face = "bold", size = 16),
          axis.text.y = element_text(size = 11)
        )
      
      print(p)
      
      ggsave(
        paste0(out_dir, "go_enrichment/GO_BP_barplot_cluster_", cl, "_down_normalized.pdf"),
        plot = p,
        width = 10,
        height = 8
      )
    }
  }
}

cat("GO bar plots created\n")

# ============================================================================
# GO VISUALIZATION 3: COMPARISON PLOTS
# ============================================================================

cat("Creating GO comparison plots...\n")

for (cl in comp_clusters) {
  go_key_up <- paste0(cl, "_up_normalized_BP")
  go_key_down <- paste0(cl, "_down_normalized_BP")
  
  has_up <- !is.null(go_results[[go_key_up]]) && nrow(go_results[[go_key_up]]) > 0
  has_down <- !is.null(go_results[[go_key_down]]) && nrow(go_results[[go_key_down]]) > 0
  
  if (has_up && has_down) {
    up_df <- as.data.frame(go_results[[go_key_up]])[1:10, ]
    down_df <- as.data.frame(go_results[[go_key_down]])[1:10, ]
    
    up_df$Direction <- "UP-regulated\n(normalized)"
    down_df$Direction <- "DOWN-regulated\n(normalized)"
    
    combined_df <- rbind(up_df, down_df)
    combined_df$neglog10p <- -log10(combined_df$p.adjust)
    
    p <- ggplot(combined_df, aes(x = reorder(Description, neglog10p),
                                 y = neglog10p,
                                 fill = Direction)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      facet_wrap(~Direction, scales = "free_y") +
      scale_fill_manual(values = c(
        "UP-regulated\n(normalized)" = "#3498DB",
        "DOWN-regulated\n(normalized)" = "#E67E22"
      )) +
      labs(
        title = paste0("Cluster ", cl, ": GO-BP Enrichment Comparison"),
        subtitle = "Top enriched terms in normalized genes",
        x = "GO Term",
        y = "-log10(adjusted p-value)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 13, color = "gray40"),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "none",
        axis.text.y = element_text(size = 10)
      )
    
    print(p)
    
    ggsave(
      paste0(out_dir, "go_enrichment/GO_BP_comparison_cluster_", cl, ".pdf"),
      plot = p,
      width = 14,
      height = 10
    )
  }
}

cat("GO comparison plots created\n")

# ============================================================================
# GO VISUALIZATION 4: SUMMARY HEATMAP ACROSS CLUSTERS
# ============================================================================

cat("Creating GO summary heatmap...\n")

summary_go_terms <- list()

for (cl in comp_clusters) {
  # UP-regulated
  go_key_up <- paste0(cl, "_up_normalized_BP")
  if (!is.null(go_results[[go_key_up]]) && nrow(go_results[[go_key_up]]) > 0) {
    top_terms <- as.data.frame(go_results[[go_key_up]])[1:5, c("Description", "p.adjust")]
    top_terms$Cluster <- cl
    top_terms$Direction <- "UP"
    summary_go_terms[[paste0(cl, "_up")]] <- top_terms
  }
  
  # DOWN-regulated
  go_key_down <- paste0(cl, "_down_normalized_BP")
  if (!is.null(go_results[[go_key_down]]) && nrow(go_results[[go_key_down]]) > 0) {
    top_terms <- as.data.frame(go_results[[go_key_down]])[1:5, c("Description", "p.adjust")]
    top_terms$Cluster <- cl
    top_terms$Direction <- "DOWN"
    summary_go_terms[[paste0(cl, "_down")]] <- top_terms
  }
}

if (length(summary_go_terms) > 0) {
  all_terms <- do.call(rbind, summary_go_terms)
  all_terms$neglog10p <- -log10(all_terms$p.adjust)
  
  term_counts <- table(all_terms$Description)
  common_terms <- names(term_counts[term_counts >= 2])
  
  if (length(common_terms) > 0) {
    filtered_terms <- all_terms[all_terms$Description %in% common_terms, ]
    
    p <- ggplot(filtered_terms, aes(x = Cluster, y = Description, fill = neglog10p)) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradient2(
        low = "white",
        mid = "yellow",
        high = "red",
        midpoint = 2,
        name = "-log10(p.adj)"
      ) +
      facet_wrap(~Direction, scales = "free_x") +
      labs(
        title = "Common GO Terms Across Clusters",
        subtitle = "Normalized genes (terms appearing in 2+ clusters)",
        x = "Cluster",
        y = "GO Term"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 13, color = "gray40"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "right"
      )
    
    print(p)
    
    ggsave(
      paste0(out_dir, "go_enrichment/GO_summary_heatmap_normalized.pdf"),
      plot = p,
      width = 12,
      height = 10
    )
  }
}

cat("GO summary heatmap created\n\n")

# ============================================================================
# EXPORT GO RESULTS TO CSV FILES
# ============================================================================

cat("Exporting GO results to CSV...\n")

for (result_name in names(go_results)) {
  if (!is.null(go_results[[result_name]]) && nrow(go_results[[result_name]]) > 0) {
    write.csv(
      as.data.frame(go_results[[result_name]]),
      file = paste0(out_dir, "go_enrichment/", result_name, ".csv"),
      row.names = FALSE
    )
  }
}

cat("GO results exported\n\n")

# ============================================================================
# SAVE UPDATED BULK DATA
# ============================================================================

cat("Saving updated analysis object...\n")
save(bulk_data, file = paste0(out_dir, "bulk_data_with_treatment_analysis.rda"))
cat("Analysis object saved\n\n")

# ============================================================================
# FINAL SUMMARY AND COMPLETION
# ============================================================================

cat("\n============================================================================\n")
cat("ANALYSIS SUMMARY\n")
cat("============================================================================\n")

cat("\nClusters Analyzed:", length(comp_clusters), "\n")
cat("Clusters:", paste(comp_clusters, collapse = ", "), "\n\n")

cat("Treatment Effect Summary:\n")
print(treatment_summary)

cat("\n")

total_normalized <- sum(treatment_summary$vcp_up_normalized) + sum(treatment_summary$vcp_down_normalized)
total_exacerbated <- sum(treatment_summary$vcp_up_exacerbated) + sum(treatment_summary$vcp_down_exacerbated)
total_unaffected <- sum(treatment_summary$vcp_up_unaffected) + sum(treatment_summary$vcp_down_unaffected)
total_vcp_genes <- sum(treatment_summary$vcp_up_total) + sum(treatment_summary$vcp_down_total)

cat("Overall Treatment Effect (All Clusters):\n")
cat("  Total VCP-affected genes:", total_vcp_genes, "\n")
cat("  Normalized by treatment:", total_normalized, 
    sprintf("(%.1f%%)", (total_normalized/total_vcp_genes)*100), "\n")
cat("  Exacerbated by treatment:", total_exacerbated,
    sprintf("(%.1f%%)", (total_exacerbated/total_vcp_genes)*100), "\n")
cat("  Unaffected by treatment:", total_unaffected,
    sprintf("(%.1f%%)", (total_unaffected/total_vcp_genes)*100), "\n\n")

cat("GO Enrichment Results:\n")
cat("  Total GO analyses:", length(go_results), "\n")
cat("  Analyses with significant terms:", 
    sum(sapply(go_results, function(x) !is.null(x) && nrow(x) > 0)), "\n")

cat("\n============================================================================\n")
cat("OUTPUT FILES LOCATION\n")
cat("============================================================================\n")
cat("Main output directory:", out_dir, "\n")
cat("  - treatment_analysis/   : Summary statistics and treatment impact data\n")
cat("  - visualizations/        : Bar plots and scatter plots\n")
cat("  - go_enrichment/         : GO enrichment results and plots\n")

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
