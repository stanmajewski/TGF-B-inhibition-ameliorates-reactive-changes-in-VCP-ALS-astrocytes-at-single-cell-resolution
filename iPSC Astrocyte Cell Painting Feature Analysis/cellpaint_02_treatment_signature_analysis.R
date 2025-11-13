#!/usr/bin/env Rscript

################################################################################
# Cell Painting Treatment Effects on VCP Disease Signature
################################################################################
#
# Purpose:
#   Focused analysis of treatment effects specifically on features that define
#   the VCP disease signature. Identifies which signature features are corrected,
#   exacerbated, unchanged, or show mixed responses to treatment. Generates
#   detailed visualizations including box plots with biological replicates.
#
# Input Files:
#   - ALS_Untreated_vs_Control_statistical_results.csv
#   - ALS_Treated_vs_Control_statistical_results.csv
#   - biological_replicates_data.csv (optional, for box plots)
#
# Output Files:
#   - PNG plots (treatment effects, box plots, heatmaps)
#   - CSV files with VCP signature analysis and control comparisons
#   - Summary statistics by treatment impact category
#
# Author: Analysis Pipeline
# Date: 2024
#
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(pheatmap)
  library(RColorBrewer)
  library(VennDiagram)
  library(ggrepel)
  library(scales)
  library(grid)
  library(optparse)
})

################################################################################
# Command-line Arguments
################################################################################

option_list <- list(
  make_option(c("-u", "--untreated"), type = "character",
              default = "ALS_Untreated_vs_Control_statistical_results.csv",
              help = "Path to untreated vs control statistical results [default: %default]",
              metavar = "FILE"),
  
  make_option(c("-t", "--treated"), type = "character",
              default = "ALS_Treated_vs_Control_statistical_results.csv",
              help = "Path to treated vs control statistical results [default: %default]",
              metavar = "FILE"),
  
  make_option(c("-b", "--biological"), type = "character",
              default = NULL,
              help = "Path to biological replicates data (optional, for box plots)",
              metavar = "FILE"),
  
  make_option(c("-o", "--output"), type = "character",
              default = "cellpaint_treatment_analysis",
              help = "Output directory name [default: %default]",
              metavar = "DIR"),
  
  make_option(c("-p", "--pvalue"), type = "double",
              default = 0.05,
              help = "P-value threshold for significance [default: %default]",
              metavar = "NUMBER"),
  
  make_option(c("-e", "--effect-threshold"), type = "double",
              default = 0.1,
              help = "Effect size change threshold for meaningful change [default: %default]",
              metavar = "NUMBER")
)

opt_parser <- OptionParser(option_list = option_list,
                          description = "\nCell Painting Treatment Effects on VCP Signature")
opt <- parse_args(opt_parser)

################################################################################
# Setup and Configuration
################################################################################

# Create output directory
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
  cat(sprintf("[%s] Created output directory: %s\n", 
              Sys.time(), opt$output))
}

# Configure plotting theme
theme_set(theme_minimal() + 
            theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                  axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16, face = "bold"),
                  legend.text = element_text(size = 12),
                  legend.title = element_text(size = 14, face = "bold"),
                  strip.text = element_text(size = 14, face = "bold")))

cat(sprintf("[%s] Configuration complete\n", Sys.time()))
cat(sprintf("  P-value threshold: %.3f\n", opt$pvalue))
cat(sprintf("  Effect size threshold: %.2f\n", opt$`effect-threshold`))
cat(sprintf("  Output directory: %s\n", opt$output))

################################################################################
# Functions
################################################################################

# Preprocess data with condition labels
preprocess_data <- function(data, condition_name) {
  data$Condition <- condition_name
  data$neg_log10_p <- -log10(data$P_value)
  data$neg_log10_p_fdr <- -log10(data$P_value_FDR)
  data$abs_effect_size <- abs(data$Effect_size)
  data$Significant_FDR_logical <- as.logical(data$Significant_FDR)
  return(data)
}

# Extract feature category from feature names
extract_feature_category <- function(feature_names) {
  categories <- rep("Other", length(feature_names))
  
  categories[grepl("Nuclei", feature_names, ignore.case = TRUE)] <- "Nuclear"
  categories[grepl("Cytoplasm", feature_names, ignore.case = TRUE)] <- "Cytoplasmic"
  categories[grepl("Mitochondria", feature_names, ignore.case = TRUE)] <- "Mitochondrial"
  categories[grepl("Granularity", feature_names, ignore.case = TRUE)] <- "Texture/Granularity"
  categories[grepl("Intensity", feature_names, ignore.case = TRUE)] <- "Intensity"
  categories[grepl("Texture", feature_names, ignore.case = TRUE)] <- "Texture"
  categories[grepl("Correlation", feature_names, ignore.case = TRUE)] <- "Correlation"
  
  return(categories)
}

# Create readable feature labels
create_feature_label <- function(feature_name, max_length = 50) {
  label <- gsub("_", " ", feature_name)
  if (nchar(label) > max_length) {
    label <- paste0(substr(label, 1, max_length), "...")
  }
  return(label)
}

################################################################################
# Load and Preprocess Data
################################################################################

cat(sprintf("\n[%s] Loading input data...\n", Sys.time()))

# Check if input files exist
if (!file.exists(opt$untreated)) {
  stop(sprintf("Untreated data file not found: %s", opt$untreated))
}
if (!file.exists(opt$treated)) {
  stop(sprintf("Treated data file not found: %s", opt$treated))
}

# Load statistical results
untreated_data <- read.csv(opt$untreated, stringsAsFactors = FALSE)
treated_data <- read.csv(opt$treated, stringsAsFactors = FALSE)

cat(sprintf("  Loaded untreated data: %d features\n", nrow(untreated_data)))
cat(sprintf("  Loaded treated data: %d features\n", nrow(treated_data)))

# Load biological replicates if provided
biological_replicates <- NULL
if (!is.null(opt$biological) && file.exists(opt$biological)) {
  biological_replicates <- read.csv(opt$biological, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded biological replicates: %d samples\n", nrow(biological_replicates)))
}

# Preprocess datasets
untreated_processed <- preprocess_data(untreated_data, "Untreated")
treated_processed <- preprocess_data(treated_data, "Treated")

################################################################################
# Part 1: Define VCP Disease Signature
################################################################################

cat(sprintf("\n[%s] ===== DEFINING VCP DISEASE SIGNATURE =====\n", Sys.time()))

# Summary statistics for untreated condition
n_total <- nrow(untreated_processed)
n_sig <- sum(untreated_processed$P_value < opt$pvalue)
n_sig_fdr <- sum(untreated_processed$Significant_FDR_logical, na.rm = TRUE)

cat(sprintf("  Total features: %d\n", n_total))
cat(sprintf("  Significant (p < %.3f): %d (%.1f%%)\n", 
            opt$pvalue, n_sig, 100 * n_sig/n_total))
cat(sprintf("  Significant (FDR corrected): %d (%.1f%%)\n", 
            n_sig_fdr, 100 * n_sig_fdr/n_total))

# Define VCP signature: features significantly different in untreated
vcp_signature <- untreated_processed[untreated_processed$P_value < opt$pvalue, ]
vcp_signature <- vcp_signature[order(vcp_signature$P_value), ]
vcp_signature$direction <- ifelse(vcp_signature$Mean_difference > 0, 
                                  "Increased", "Decreased")
vcp_signature$category <- extract_feature_category(vcp_signature$Feature)

cat(sprintf("\n  VCP Disease Signature: %d significantly differential features\n", 
            nrow(vcp_signature)))

# Direction breakdown
n_increased <- sum(vcp_signature$direction == "Increased")
n_decreased <- sum(vcp_signature$direction == "Decreased")
cat(sprintf("    Increased in VCP: %d (%.1f%%)\n", 
            n_increased, 100 * n_increased/nrow(vcp_signature)))
cat(sprintf("    Decreased in VCP: %d (%.1f%%)\n", 
            n_decreased, 100 * n_decreased/nrow(vcp_signature)))

# Top features
cat(sprintf("\n  Top 10 features defining VCP signature:\n"))
top10 <- vcp_signature[1:min(10, nrow(vcp_signature)), 
                       c("Feature", "P_value", "Mean_difference", "Effect_size", "direction")]
print(top10, row.names = FALSE)

################################################################################
# Part 2: Treatment Effects on VCP Signature Features
################################################################################

cat(sprintf("\n[%s] ===== TREATMENT EFFECTS ON VCP SIGNATURE =====\n", Sys.time()))

# Filter both datasets to VCP signature features only
vcp_features <- vcp_signature$Feature
untreated_vcp <- untreated_processed[untreated_processed$Feature %in% vcp_features, ]
treated_vcp <- treated_processed[treated_processed$Feature %in% vcp_features, ]

# Merge for comparison
comparison_data <- merge(untreated_vcp, treated_vcp, 
                        by = "Feature", suffixes = c("_untreated", "_treated"))
comparison_data <- comparison_data[complete.cases(comparison_data), ]

cat(sprintf("  Analyzing treatment effects on %d VCP signature features\n", 
            nrow(comparison_data)))

# Calculate treatment effect metrics
comparison_data$effect_size_change <- comparison_data$Effect_size_treated - 
                                      comparison_data$Effect_size_untreated
comparison_data$p_value_change <- comparison_data$neg_log10_p_treated - 
                                  comparison_data$neg_log10_p_untreated
comparison_data$abs_effect_change <- abs(comparison_data$effect_size_change)

# Classify treatment impact using refined criteria
effect_threshold <- opt$`effect-threshold`
comparison_data$treatment_impact <- "Neutral"

# CORRECTIVE: Treatment normalizes toward control
# Concordant metrics (both p-value and effect size improve) with meaningful change
corrective_concordant <- (
  (comparison_data$P_value_treated > comparison_data$P_value_untreated) &
  (abs(comparison_data$Effect_size_treated) < abs(comparison_data$Effect_size_untreated)) &
  (comparison_data$abs_effect_change >= effect_threshold)
)

# OR direction reversal (strong normalization signal)
corrective_reversal <- (
  (sign(comparison_data$Effect_size_untreated) != sign(comparison_data$Effect_size_treated)) &
  (comparison_data$abs_effect_change >= effect_threshold)
)

corrective_idx <- corrective_concordant | corrective_reversal
comparison_data$treatment_impact[corrective_idx] <- "Corrective"

# EXACERBATING: Treatment worsens away from control
exacerbating_idx <- (
  (comparison_data$P_value_treated < comparison_data$P_value_untreated) &
  (abs(comparison_data$Effect_size_treated) > abs(comparison_data$Effect_size_untreated)) &
  (comparison_data$abs_effect_change >= effect_threshold)
)
comparison_data$treatment_impact[exacerbating_idx] <- "Exacerbating"

# UNCHANGED: Minimal effect size change
unchanged_idx <- (comparison_data$abs_effect_change < effect_threshold)
comparison_data$treatment_impact[unchanged_idx] <- "Unchanged"

# MIXED: Remaining features
mixed_idx <- !(corrective_idx | exacerbating_idx | unchanged_idx)
comparison_data$treatment_impact[mixed_idx] <- "Mixed"

# Add categories
comparison_data$category <- extract_feature_category(comparison_data$Feature)

# Treatment impact summary
treatment_summary <- table(comparison_data$treatment_impact)
cat(sprintf("\n  Treatment Impact Classification:\n"))
for (impact in names(treatment_summary)) {
  cat(sprintf("    %s: %d features (%.1f%%)\n", 
              impact, treatment_summary[impact], 
              100 * treatment_summary[impact]/nrow(comparison_data)))
}

# Extract features by category
corrective_features <- comparison_data[comparison_data$treatment_impact == "Corrective", ]
corrective_features <- corrective_features[order(abs(corrective_features$effect_size_change), 
                                                 decreasing = TRUE), ]

exacerbating_features <- comparison_data[comparison_data$treatment_impact == "Exacerbating", ]
exacerbating_features <- exacerbating_features[order(abs(exacerbating_features$effect_size_change), 
                                                      decreasing = TRUE), ]

unchanged_features <- comparison_data[comparison_data$treatment_impact == "Unchanged", ]

mixed_features <- comparison_data[comparison_data$treatment_impact == "Mixed", ]

# Display top features
if (nrow(corrective_features) > 0) {
  cat(sprintf("\n  Top 5 corrective features:\n"))
  print(head(corrective_features[, c("Feature", "Effect_size_untreated", "Effect_size_treated", 
                                     "P_value_untreated", "P_value_treated")], 5), 
        row.names = FALSE)
}

if (nrow(exacerbating_features) > 0) {
  cat(sprintf("\n  Top 5 exacerbating features:\n"))
  print(head(exacerbating_features[, c("Feature", "Effect_size_untreated", "Effect_size_treated",
                                       "P_value_untreated", "P_value_treated")], 5), 
        row.names = FALSE)
}

################################################################################
# Generate Visualizations
################################################################################

cat(sprintf("\n[%s] Generating visualizations...\n", Sys.time()))

# Plot 1: Treatment impact overview
treatment_counts <- as.data.frame(treatment_summary)
names(treatment_counts) <- c("Impact", "Count")

impact_colors <- c("Corrective" = "#4DAF4A", "Exacerbating" = "#E41A1C", 
                   "Unchanged" = "#377EB8", "Mixed" = "#984EA3", "Neutral" = "grey70")

p1_overview <- ggplot(treatment_counts, aes(x = reorder(Impact, Count), y = Count, fill = Impact)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", Count, 100*Count/sum(Count))), 
            hjust = -0.1, size = 5, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = impact_colors) +
  labs(title = "Treatment Impact on VCP Disease Signature",
       subtitle = sprintf("%d features analyzed (effect size threshold: %.2f)", 
                         nrow(comparison_data), effect_threshold),
       x = "Treatment Impact Category", y = "Number of Features") +
  theme(legend.position = "none")

ggsave(file.path(opt$output, "01_treatment_impact_overview.png"), 
       p1_overview, width = 12, height = 8, dpi = 300)
cat("  Saved: 01_treatment_impact_overview.png\n")

# Plot 2: Effect size changes
p2_effect_change <- ggplot(comparison_data, 
                           aes(x = Effect_size_untreated, y = Effect_size_treated, 
                               color = treatment_impact)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "solid", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "solid", alpha = 0.3) +
  scale_color_manual(values = impact_colors, name = "Treatment Impact") +
  labs(title = "Effect Size Changes: Untreated vs Treated",
       subtitle = "Points below diagonal = treatment reduces effect size",
       x = "Effect Size (VCP Untreated vs Control)",
       y = "Effect Size (VCP Treated vs Control)") +
  theme(legend.position = "right")

ggsave(file.path(opt$output, "02_effect_size_changes.png"), 
       p2_effect_change, width = 12, height = 10, dpi = 300)
cat("  Saved: 02_effect_size_changes.png\n")

# Plot 3: P-value changes
p3_pvalue_change <- ggplot(comparison_data, 
                          aes(x = neg_log10_p_untreated, y = neg_log10_p_treated, 
                              color = treatment_impact)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(opt$pvalue), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -log10(opt$pvalue), linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = impact_colors, name = "Treatment Impact") +
  labs(title = "Statistical Significance Changes",
       subtitle = "Points below diagonal = treatment increases p-value (less significant)",
       x = "-log₁₀(p-value) VCP Untreated vs Control",
       y = "-log₁₀(p-value) VCP Treated vs Control") +
  theme(legend.position = "right")

ggsave(file.path(opt$output, "03_pvalue_changes.png"), 
       p3_pvalue_change, width = 12, height = 10, dpi = 300)
cat("  Saved: 03_pvalue_changes.png\n")

# Plot 4: Treatment impact by category
category_impact <- comparison_data %>%
  group_by(category, treatment_impact) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(category != "Other")

if (nrow(category_impact) > 0) {
  p4_category <- ggplot(category_impact, 
                        aes(x = reorder(category, count), y = count, fill = treatment_impact)) +
    geom_col(position = "stack", alpha = 0.8) +
    geom_text(aes(label = count), position = position_stack(vjust = 0.5), 
              size = 4, fontface = "bold", color = "white") +
    coord_flip() +
    scale_fill_manual(values = impact_colors, name = "Treatment Impact") +
    labs(title = "Treatment Impact by Feature Category",
         subtitle = "Distribution across cellular compartments and measurement types",
         x = "Feature Category", y = "Number of Features") +
    theme(legend.position = "right")
  
  ggsave(file.path(opt$output, "04_impact_by_category.png"), 
         p4_category, width = 14, height = 10, dpi = 300)
  cat("  Saved: 04_impact_by_category.png\n")
}

# Plot 5: Top corrective features
if (nrow(corrective_features) > 0) {
  top_corrective <- head(corrective_features, 20)
  top_corrective$Feature_short <- create_feature_label(top_corrective$Feature, 50)
  
  # Reshape for grouped bar plot
  top_corr_long <- top_corrective %>%
    select(Feature_short, Effect_size_untreated, Effect_size_treated) %>%
    pivot_longer(cols = c(Effect_size_untreated, Effect_size_treated),
                 names_to = "Condition", values_to = "Effect_Size")
  
  top_corr_long$Condition <- factor(top_corr_long$Condition,
                                    levels = c("Effect_size_untreated", "Effect_size_treated"),
                                    labels = c("Untreated", "Treated"))
  
  p5_top_corr <- ggplot(top_corr_long, aes(x = reorder(Feature_short, Effect_Size), 
                                            y = Effect_Size, fill = Condition)) +
    geom_col(position = "dodge", alpha = 0.8) +
    coord_flip() +
    scale_fill_manual(values = c("Untreated" = "#E41A1C", "Treated" = "#4DAF4A")) +
    labs(title = "Top 20 Corrective Features",
         subtitle = "Treatment reduces VCP disease signature (effect size comparison)",
         x = "", y = "Effect Size vs Control") +
    theme(axis.text.y = element_text(size = 10),
          legend.position = "bottom")
  
  ggsave(file.path(opt$output, "05_top_corrective_features.png"), 
         p5_top_corr, width = 14, height = 12, dpi = 300)
  cat("  Saved: 05_top_corrective_features.png\n")
}

# Plot 6: Top exacerbating features
if (nrow(exacerbating_features) > 0) {
  top_exac <- head(exacerbating_features, 20)
  top_exac$Feature_short <- create_feature_label(top_exac$Feature, 50)
  
  top_exac_long <- top_exac %>%
    select(Feature_short, Effect_size_untreated, Effect_size_treated) %>%
    pivot_longer(cols = c(Effect_size_untreated, Effect_size_treated),
                 names_to = "Condition", values_to = "Effect_Size")
  
  top_exac_long$Condition <- factor(top_exac_long$Condition,
                                    levels = c("Effect_size_untreated", "Effect_size_treated"),
                                    labels = c("Untreated", "Treated"))
  
  p6_top_exac <- ggplot(top_exac_long, aes(x = reorder(Feature_short, Effect_Size), 
                                            y = Effect_Size, fill = Condition)) +
    geom_col(position = "dodge", alpha = 0.8) +
    coord_flip() +
    scale_fill_manual(values = c("Untreated" = "#984EA3", "Treated" = "#E41A1C")) +
    labs(title = "Top 20 Exacerbating Features",
         subtitle = "Treatment worsens VCP disease signature (effect size comparison)",
         x = "", y = "Effect Size vs Control") +
    theme(axis.text.y = element_text(size = 10),
          legend.position = "bottom")
  
  ggsave(file.path(opt$output, "06_top_exacerbating_features.png"), 
         p6_top_exac, width = 14, height = 12, dpi = 300)
  cat("  Saved: 06_top_exacerbating_features.png\n")
}

################################################################################
# Export Results
################################################################################

cat(sprintf("\n[%s] Exporting results...\n", Sys.time()))

# VCP signature definition
write.csv(vcp_signature, 
          file.path(opt$output, "vcp_disease_signature.csv"), 
          row.names = FALSE)
cat("  Saved: vcp_disease_signature.csv\n")

# Treatment comparison results
write.csv(comparison_data, 
          file.path(opt$output, "treatment_effects_on_signature.csv"), 
          row.names = FALSE)
cat("  Saved: treatment_effects_on_signature.csv\n")

# Features by impact category
if (nrow(corrective_features) > 0) {
  write.csv(corrective_features, 
            file.path(opt$output, "corrective_features.csv"), 
            row.names = FALSE)
  cat("  Saved: corrective_features.csv\n")
}

if (nrow(exacerbating_features) > 0) {
  write.csv(exacerbating_features, 
            file.path(opt$output, "exacerbating_features.csv"), 
            row.names = FALSE)
  cat("  Saved: exacerbating_features.csv\n")
}

if (nrow(unchanged_features) > 0) {
  write.csv(unchanged_features, 
            file.path(opt$output, "unchanged_features.csv"), 
            row.names = FALSE)
  cat("  Saved: unchanged_features.csv\n")
}

if (nrow(mixed_features) > 0) {
  write.csv(mixed_features, 
            file.path(opt$output, "mixed_features.csv"), 
            row.names = FALSE)
  cat("  Saved: mixed_features.csv\n")
}

# Summary by category
category_summary <- comparison_data %>%
  group_by(category, treatment_impact) %>%
  summarise(
    count = n(),
    mean_effect_change = mean(abs(effect_size_change)),
    median_pval_untreated = median(P_value_untreated),
    median_pval_treated = median(P_value_treated),
    .groups = "drop"
  ) %>%
  arrange(category, desc(count))

write.csv(category_summary, 
          file.path(opt$output, "treatment_impact_by_category.csv"), 
          row.names = FALSE)
cat("  Saved: treatment_impact_by_category.csv\n")

################################################################################
# Final Summary
################################################################################

cat(sprintf("\n[%s] ===== ANALYSIS COMPLETE =====\n", Sys.time()))

cat("\nVCP Disease Signature:\n")
cat(sprintf("  %d features define the VCP signature (p < %.3f)\n", 
            nrow(vcp_signature), opt$pvalue))
cat(sprintf("  %d increased, %d decreased in VCP vs control\n", 
            n_increased, n_decreased))

cat("\nTreatment Effects on VCP Signature:\n")
for (impact in names(treatment_summary)) {
  cat(sprintf("  %s: %d features (%.1f%%)\n", 
              impact, treatment_summary[impact], 
              100 * treatment_summary[impact]/nrow(comparison_data)))
}

cat("\nClassification Criteria:\n")
cat(sprintf("  Effect size change threshold: %.2f\n", effect_threshold))
cat("  CORRECTIVE: p-value increases AND effect size decreases (with meaningful change)\n")
cat("              OR direction reversal with meaningful change\n")
cat("  EXACERBATING: p-value decreases AND effect size increases (with meaningful change)\n")
cat("  UNCHANGED: Effect size change below threshold\n")
cat("  MIXED: Other patterns\n")

cat(sprintf("\nOutput files saved to: %s\n", opt$output))
cat("  - Visualization plots (PNG)\n")
cat("  - Feature lists by treatment impact (CSV)\n")
cat("  - Category-level summaries (CSV)\n")

cat(sprintf("\n[%s] Pipeline completed successfully\n", Sys.time()))
