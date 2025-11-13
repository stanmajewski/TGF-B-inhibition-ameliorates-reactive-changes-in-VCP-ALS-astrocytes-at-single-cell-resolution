#!/usr/bin/env Rscript

################################################################################
# Cell Painting Category-Specific Analysis
################################################################################
#
# Purpose:
#   Focused analysis of Cell Painting features organized by CellProfiler
#   measurement categories (AreaShape, Granularity, Intensity, Texture, etc.).
#   Examines how different feature types respond to VCP disease and treatment,
#   generating category-specific visualizations and statistics.
#
# Input Files:
#   - ALS_Untreated_vs_Control_statistical_results.csv
#   - ALS_Treated_vs_Control_statistical_results.csv
#
# Output Files:
#   - Category-specific volcano plots and heatmaps
#   - Treatment impact analysis by category
#   - CSV files with category-level statistics
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
  
  make_option(c("-o", "--output"), type = "character",
              default = "cellpaint_category_analysis",
              help = "Output directory name [default: %default]",
              metavar = "DIR"),
  
  make_option(c("-p", "--pvalue"), type = "double",
              default = 0.05,
              help = "P-value threshold for significance [default: %default]",
              metavar = "NUMBER"),
  
  make_option(c("-n", "--top-features"), type = "integer",
              default = 3,
              help = "Number of top features per category to display [default: %default]",
              metavar = "NUMBER")
)

opt_parser <- OptionParser(option_list = option_list,
                          description = "\nCell Painting Category-Specific Analysis")
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
            theme(plot.title = element_text(hjust = 0.5, size = 23, face = "bold"),
                  axis.text = element_text(size = 22),
                  axis.title = element_text(size = 22, face = "bold"),
                  legend.text = element_text(size = 20),
                  legend.title = element_text(size = 22, face = "bold"),
                  strip.text = element_text(size = 20, face = "bold")))

cat(sprintf("[%s] Configuration complete\n", Sys.time()))
cat(sprintf("  P-value threshold: %.3f\n", opt$pvalue))
cat(sprintf("  Top features per category: %d\n", opt$`top-features`))
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

# Extract CellProfiler measurement category from feature names
extract_feature_category <- function(feature_names) {
  categories <- rep("Other", length(feature_names))
  
  # Define categories based on CellProfiler measurement types
  # Order matters - more specific patterns first
  categories[grepl("_AreaShape_", feature_names, ignore.case = TRUE)] <- "AreaShape"
  categories[grepl("_Granularity_", feature_names, ignore.case = TRUE)] <- "Granularity"
  categories[grepl("_Intensity_", feature_names, ignore.case = TRUE)] <- "Intensity"
  categories[grepl("_RadialDistribution_", feature_names, ignore.case = TRUE)] <- "RadialDistribution"
  categories[grepl("_Correlation_", feature_names, ignore.case = TRUE)] <- "Correlation"
  categories[grepl("_ObjectSkeleton_", feature_names, ignore.case = TRUE)] <- "ObjectSkeleton"
  categories[grepl("_Texture_", feature_names, ignore.case = TRUE)] <- "Texture"
  
  # Handle cases without underscores (only if not already categorized)
  categories[categories == "Other" & grepl("AreaShape", feature_names, ignore.case = TRUE)] <- "AreaShape"
  categories[categories == "Other" & grepl("Granularity", feature_names, ignore.case = TRUE)] <- "Granularity"
  categories[categories == "Other" & grepl("Intensity", feature_names, ignore.case = TRUE)] <- "Intensity"
  categories[categories == "Other" & grepl("RadialDistribution", feature_names, ignore.case = TRUE)] <- "RadialDistribution"
  categories[categories == "Other" & grepl("Correlation", feature_names, ignore.case = TRUE)] <- "Correlation"
  categories[categories == "Other" & grepl("ObjectSkeleton", feature_names, ignore.case = TRUE)] <- "ObjectSkeleton"
  categories[categories == "Other" & grepl("Texture", feature_names, ignore.case = TRUE)] <- "Texture"
  
  return(categories)
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

# Load datasets
untreated_data <- read.csv(opt$untreated, stringsAsFactors = FALSE)
treated_data <- read.csv(opt$treated, stringsAsFactors = FALSE)

cat(sprintf("  Loaded untreated data: %d features\n", nrow(untreated_data)))
cat(sprintf("  Loaded treated data: %d features\n", nrow(treated_data)))

# Preprocess datasets
untreated_processed <- preprocess_data(untreated_data, "Untreated")
treated_processed <- preprocess_data(treated_data, "Treated")

# Add categories to all features
untreated_processed$category <- extract_feature_category(untreated_processed$Feature)
treated_processed$category <- extract_feature_category(treated_processed$Feature)

################################################################################
# Part 1: VCP-ALS Untreated vs Control by Category
################################################################################

cat(sprintf("\n[%s] ===== VCP-ALS UNTREATED VS CONTROL BY CATEGORY =====\n", Sys.time()))

# Overall statistics
n_total <- nrow(untreated_processed)
n_sig <- sum(untreated_processed$P_value < opt$pvalue)

cat(sprintf("  Total features: %d\n", n_total))
cat(sprintf("  Significant (p < %.3f): %d (%.1f%%)\n", 
            opt$pvalue, n_sig, 100 * n_sig/n_total))

# Extract significant features
untreated_significant <- untreated_processed[untreated_processed$P_value < opt$pvalue, ]
untreated_significant <- untreated_significant[order(untreated_significant$P_value), ]
untreated_significant$direction <- ifelse(untreated_significant$Mean_difference > 0, 
                                         "Increased", "Decreased")

# Category distribution
cat("\n  Feature Category Distribution (All Features):\n")
all_category_counts <- table(untreated_processed$category)
for (cat_name in names(all_category_counts)) {
  cat(sprintf("    %s: %d\n", cat_name, all_category_counts[cat_name]))
}

cat("\n  Feature Category Distribution (Significant Features):\n")
sig_category_counts <- table(untreated_significant$category)
for (cat_name in names(sig_category_counts)) {
  pct <- 100 * sig_category_counts[cat_name] / all_category_counts[cat_name]
  cat(sprintf("    %s: %d (%.1f%% of total %s features)\n", 
              cat_name, sig_category_counts[cat_name], pct, cat_name))
}

# Display sample categorization
cat("\n  Sample Feature Categorization:\n")
sample_features <- head(untreated_significant[, c("Feature", "category")], 10)
print(sample_features, row.names = FALSE)

# Top features by category
cat("\n  Top Features per Category:\n")
for (cat_name in names(sig_category_counts)) {
  cat_features <- untreated_significant[untreated_significant$category == cat_name, ]
  if (nrow(cat_features) > 0) {
    cat(sprintf("\n    %s (top %d):\n", cat_name, min(opt$`top-features`, nrow(cat_features))))
    top_cat <- head(cat_features[, c("Feature", "P_value", "Effect_size", "direction")], 
                   opt$`top-features`)
    print(top_cat, row.names = FALSE)
  }
}

################################################################################
# Part 2: Treatment Effects by Category
################################################################################

cat(sprintf("\n[%s] ===== TREATMENT EFFECTS BY CATEGORY =====\n", Sys.time()))

# Merge datasets
comparison_data <- merge(untreated_processed, treated_processed, 
                        by = "Feature", suffixes = c("_untreated", "_treated"))
comparison_data <- comparison_data[complete.cases(comparison_data), ]

cat(sprintf("  Merged dataset: %d features\n", nrow(comparison_data)))

# Calculate treatment metrics
comparison_data$effect_size_change <- comparison_data$Effect_size_treated - 
                                      comparison_data$Effect_size_untreated
comparison_data$p_value_change <- comparison_data$neg_log10_p_treated - 
                                  comparison_data$neg_log10_p_untreated

# Classify treatment impact
comparison_data$treatment_impact <- "Neutral"

# Corrective
corrective_idx <- (comparison_data$P_value_untreated < opt$pvalue & 
                   comparison_data$P_value_treated >= comparison_data$P_value_untreated)
comparison_data$treatment_impact[corrective_idx] <- "Corrective"

# Exacerbating
exacerbating_idx <- (comparison_data$P_value_untreated >= opt$pvalue & 
                     comparison_data$P_value_treated < opt$pvalue)
comparison_data$treatment_impact[exacerbating_idx] <- "Exacerbating"

# Both significant, treated worse
both_sig_worse <- (comparison_data$P_value_untreated < opt$pvalue & 
                   comparison_data$P_value_treated < opt$pvalue & 
                   sign(comparison_data$Effect_size_untreated) == sign(comparison_data$Effect_size_treated) &
                   abs(comparison_data$Effect_size_treated) > abs(comparison_data$Effect_size_untreated))
comparison_data$treatment_impact[both_sig_worse] <- "Exacerbating"

# Both significant, treated better
both_sig_better <- (comparison_data$P_value_untreated < opt$pvalue & 
                    comparison_data$P_value_treated < opt$pvalue & 
                    sign(comparison_data$Effect_size_untreated) == sign(comparison_data$Effect_size_treated) &
                    abs(comparison_data$Effect_size_treated) < abs(comparison_data$Effect_size_untreated))
comparison_data$treatment_impact[both_sig_better] <- "Corrective"

# Use category from untreated analysis
comparison_data$category <- comparison_data$category_untreated

# Treatment summary by category
treatment_by_category <- comparison_data %>%
  filter(treatment_impact %in% c("Corrective", "Exacerbating")) %>%
  group_by(category, treatment_impact) %>%
  summarise(
    count = n(),
    mean_effect_change = mean(abs(effect_size_change)),
    .groups = "drop"
  ) %>%
  arrange(category, treatment_impact)

cat("\n  Treatment Effects by Category:\n")
print(treatment_by_category, row.names = FALSE)

# Extract corrective and exacerbating features
corrective_features <- comparison_data[comparison_data$treatment_impact == "Corrective", ]
corrective_features <- corrective_features[order(abs(corrective_features$effect_size_change), 
                                                 decreasing = TRUE), ]

exacerbating_features <- comparison_data[comparison_data$treatment_impact == "Exacerbating", ]
exacerbating_features <- exacerbating_features[order(abs(exacerbating_features$effect_size_change), 
                                                      decreasing = TRUE), ]

################################################################################
# Generate Visualizations
################################################################################

cat(sprintf("\n[%s] Generating visualizations...\n", Sys.time()))

# Plot 1: Volcano plot colored by category
category_colors <- RColorBrewer::brewer.pal(
  min(length(unique(untreated_processed$category)), 8), "Set2")
names(category_colors) <- unique(untreated_processed$category)[1:length(category_colors)]

p1_volcano_category <- ggplot(untreated_processed, 
                              aes(x = Mean_difference, y = neg_log10_p, color = category)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_hline(yintercept = -log10(opt$pvalue), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
  scale_color_manual(values = category_colors, name = "Feature Category") +
  labs(title = "Differential Features by Category",
       subtitle = "VCP-ALS Untreated vs Control - colored by measurement type",
       x = "Mean Difference (Effect Size)",
       y = "-log₁₀(p-value)") +
  theme(legend.position = "right")

ggsave(file.path(opt$output, "01_volcano_plot_by_category.png"), 
       p1_volcano_category, width = 16, height = 10, dpi = 300)
cat("  Saved: 01_volcano_plot_by_category.png\n")

# Plot 2: Category distribution
category_dist_data <- untreated_significant %>%
  group_by(category, direction) %>%
  summarise(count = n(), .groups = "drop")

p2_category_dist <- ggplot(category_dist_data, 
                           aes(x = reorder(category, count), y = count, fill = direction)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(aes(label = count), position = position_dodge(width = 0.9), 
            hjust = -0.2, size = 6, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("Increased" = "firebrick", "Decreased" = "steelblue"),
                   name = "Direction") +
  labs(title = "Significant Features by Category",
       subtitle = sprintf("Distribution of up/down-regulated features (p < %.3f)", opt$pvalue),
       x = "Feature Category", y = "Number of Significant Features") +
  theme(legend.position = "bottom")

ggsave(file.path(opt$output, "02_category_distribution.png"), 
       p2_category_dist, width = 12, height = 8, dpi = 300)
cat("  Saved: 02_category_distribution.png\n")

# Plot 3: Top features per category (faceted)
top_per_category <- untreated_significant %>%
  filter(category != "Other") %>%
  group_by(category) %>%
  slice_head(n = opt$`top-features`) %>%
  ungroup()

if (nrow(top_per_category) > 0) {
  top_per_category$Feature_short <- substr(gsub("_", " ", top_per_category$Feature), 1, 40)
  
  p3_top_features <- ggplot(top_per_category, 
                            aes(x = reorder(Feature_short, abs(Effect_size)), 
                                y = abs(Effect_size), fill = direction)) +
    geom_col(alpha = 0.8) +
    coord_flip() +
    facet_wrap(~ category, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("Increased" = "firebrick", "Decreased" = "steelblue")) +
    labs(title = sprintf("Top %d Features per Category", opt$`top-features`),
         subtitle = "Highest effect size features in each measurement category",
         x = "", y = "Absolute Effect Size",
         fill = "Direction") +
    theme(axis.text.y = element_text(size = 10),
          strip.text = element_text(size = 12, face = "bold"))
  
  ggsave(file.path(opt$output, "03_top_features_per_category.png"), 
         p3_top_features, width = 16, height = 14, dpi = 300)
  cat("  Saved: 03_top_features_per_category.png\n")
}

# Plot 4: Treatment impact by category
if (nrow(treatment_by_category) > 0) {
  p4_treatment_category <- ggplot(treatment_by_category, 
                                  aes(x = reorder(category, count), y = count, 
                                      fill = treatment_impact)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(aes(label = count), position = position_dodge(width = 0.9), 
              hjust = -0.2, size = 6) +
    coord_flip() +
    scale_fill_manual(values = c("Corrective" = "green", "Exacerbating" = "red"),
                     name = "Treatment Impact") +
    labs(title = "Treatment Impact by Feature Category",
         subtitle = "How treatment affects different measurement types",
         x = "Feature Category", y = "Number of Features") +
    theme(legend.position = "bottom")
  
  ggsave(file.path(opt$output, "04_treatment_impact_by_category.png"), 
         p4_treatment_category, width = 12, height = 8, dpi = 300)
  cat("  Saved: 04_treatment_impact_by_category.png\n")
}

# Plot 5: Effect size distribution by category
effect_by_category <- untreated_significant %>%
  filter(category != "Other")

if (nrow(effect_by_category) > 0) {
  p5_effect_dist <- ggplot(effect_by_category, 
                           aes(x = category, y = abs(Effect_size), fill = category)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    scale_fill_manual(values = category_colors) +
    labs(title = "Effect Size Distribution by Category",
         subtitle = "Range of effect sizes across measurement types",
         x = "Feature Category", y = "Absolute Effect Size") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  ggsave(file.path(opt$output, "05_effect_size_by_category.png"), 
         p5_effect_dist, width = 12, height = 8, dpi = 300)
  cat("  Saved: 05_effect_size_by_category.png\n")
}

# Plot 6: Corrective features heatmap by category
if (nrow(corrective_features) > 0) {
  top_corrective_per_cat <- corrective_features %>%
    filter(category != "Other") %>%
    group_by(category) %>%
    arrange(desc(abs(effect_size_change))) %>%
    slice_head(n = opt$`top-features`) %>%
    ungroup()
  
  if (nrow(top_corrective_per_cat) > 0) {
    heatmap_data <- top_corrective_per_cat %>%
      select(Feature, Effect_size_untreated, Effect_size_treated, category)
    
    heatmap_data$Feature_Clean <- substr(gsub("_", " ", heatmap_data$Feature), 1, 50)
    
    heatmap_long <- heatmap_data %>%
      select(Feature_Clean, Effect_size_untreated, Effect_size_treated, category) %>%
      pivot_longer(cols = c(Effect_size_untreated, Effect_size_treated),
                   names_to = "Condition", values_to = "EffectSize")
    
    heatmap_long$Condition <- factor(heatmap_long$Condition,
                                     levels = c("Effect_size_untreated", "Effect_size_treated"),
                                     labels = c("Untreated vs Control", "Treated vs Control"))
    
    p6_heatmap <- ggplot(heatmap_long, 
                         aes(x = Condition, y = reorder(Feature_Clean, EffectSize), 
                             fill = EffectSize)) +
      geom_tile(color = "black", size = 0.5) +
      geom_text(aes(label = round(EffectSize, 2)), color = "black", size = 4, fontface = "bold") +
      scale_fill_gradient2(low = "#0571b0", mid = "white", high = "#ca0020", 
                          midpoint = 0, name = "Effect Size") +
      facet_wrap(~ category, scales = "free_y", ncol = 2) +
      labs(title = "Top Corrective Features by Category",
           subtitle = "Effect size changes from untreated to treated conditions",
           x = "Comparison", y = "Features") +
      theme(axis.text.y = element_text(size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            strip.text = element_text(size = 12, face = "bold"))
    
    ggsave(file.path(opt$output, "06_corrective_features_by_category.png"), 
           p6_heatmap, width = 14, height = 12, dpi = 300)
    cat("  Saved: 06_corrective_features_by_category.png\n")
  }
}

# Plot 7: Faceted volcano plots by category
p7_volcano_faceted <- ggplot(untreated_processed[untreated_processed$category != "Other", ], 
                             aes(x = Mean_difference, y = neg_log10_p)) +
  geom_point(alpha = 0.5, size = 2, aes(color = P_value < opt$pvalue)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"), 
                    name = "Significant", labels = c("No", "Yes")) +
  geom_hline(yintercept = -log10(opt$pvalue), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
  facet_wrap(~ category, scales = "free", ncol = 3) +
  labs(title = "Volcano Plots by Feature Category",
       subtitle = "Distribution of effects across different measurement types",
       x = "Mean Difference (Effect Size)",
       y = "-log₁₀(p-value)") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom")

ggsave(file.path(opt$output, "07_volcano_plots_faceted.png"), 
       p7_volcano_faceted, width = 16, height = 12, dpi = 300)
cat("  Saved: 07_volcano_plots_faceted.png\n")

################################################################################
# Export Results
################################################################################

cat(sprintf("\n[%s] Exporting results...\n", Sys.time()))

# Category-specific summary
category_summary <- untreated_significant %>%
  group_by(category) %>%
  summarise(
    Total_Significant = n(),
    Upregulated = sum(direction == "Increased"),
    Downregulated = sum(direction == "Decreased"),
    Mean_Effect_Size = mean(abs(Effect_size)),
    Median_P_Value = median(P_value),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Significant))

write.csv(category_summary, 
          file.path(opt$output, "category_summary_statistics.csv"), 
          row.names = FALSE)
cat("  Saved: category_summary_statistics.csv\n")

# Treatment effects by category
write.csv(treatment_by_category, 
          file.path(opt$output, "treatment_effects_by_category.csv"), 
          row.names = FALSE)
cat("  Saved: treatment_effects_by_category.csv\n")

# Full comparison data with categories
write.csv(comparison_data, 
          file.path(opt$output, "treatment_comparison_with_categories.csv"), 
          row.names = FALSE)
cat("  Saved: treatment_comparison_with_categories.csv\n")

# Corrective features
if (nrow(corrective_features) > 0) {
  write.csv(corrective_features, 
            file.path(opt$output, "corrective_features_all.csv"), 
            row.names = FALSE)
  cat("  Saved: corrective_features_all.csv\n")
}

# Exacerbating features
if (nrow(exacerbating_features) > 0) {
  write.csv(exacerbating_features, 
            file.path(opt$output, "exacerbating_features_all.csv"), 
            row.names = FALSE)
  cat("  Saved: exacerbating_features_all.csv\n")
}

################################################################################
# Final Summary
################################################################################

cat(sprintf("\n[%s] ===== ANALYSIS COMPLETE =====\n", Sys.time()))

cat("\nCategory-Specific Summary:\n")
print(category_summary, row.names = FALSE)

cat("\nKey Findings by Category:\n")
for (cat_name in category_summary$category) {
  cat_data <- category_summary[category_summary$category == cat_name, ]
  cat(sprintf("  %s:\n", cat_name))
  cat(sprintf("    Total significant: %d\n", cat_data$Total_Significant))
  cat(sprintf("    Upregulated: %d\n", cat_data$Upregulated))
  cat(sprintf("    Downregulated: %d\n", cat_data$Downregulated))
  cat(sprintf("    Mean effect size: %.3f\n", cat_data$Mean_Effect_Size))
}

cat(sprintf("\nOutput files saved to: %s\n", opt$output))
cat("  - 7 visualization plots (PNG format)\n")
cat("  - Category-specific statistical summaries (CSV)\n")
cat("  - Full feature lists with categories (CSV)\n")

cat(sprintf("\n[%s] Pipeline completed successfully\n", Sys.time()))
