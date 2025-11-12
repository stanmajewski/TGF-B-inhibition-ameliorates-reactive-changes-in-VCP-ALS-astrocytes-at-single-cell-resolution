message("\n\n##########################################################################\n",
        "# Start 09: Cell cycle scoring and analysis ", Sys.time(),
        "\n##########################################################################\n",
        "\n   Cell cycle phase assignment and distribution analysis",
        "\n##########################################################################\n\n")

#set environment and main directory
main_dir = getwd()
setwd(main_dir)

# Open packages necessary for analysis.
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

#specify script/output index as prefix for file names
script_ind = "09_"

#specify output directory
out_dir = paste0(main_dir, "/results/")

#load dataset
load(file = paste0(out_dir, "04_seur_integr_labelled.rda")) 

###########################################################
# Cell Cycle Scoring
###########################################################

message("\n\n          *** Performing cell cycle scoring... ", Sys.time(),"\n\n")

# Load Seurat's cell cycle gene lists
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

# Score cell cycle phases
seur <- CellCycleScoring(seur, s.features = s_genes, g2m.features = g2m_genes)

# Save updated object with cell cycle scores
save(seur, file = paste0(out_dir, script_ind, "seur_with_cell_cycle.rda"))

###########################################################
# UMAP visualization by cell cycle phase
###########################################################

message("\n\n          *** Creating cell cycle UMAP visualization... ", Sys.time(),"\n\n")

# Create basic UMAP plot
pdf(file.path(out_dir, paste0(script_ind, "umap_cell_cycle_phase.pdf")), width = 10, height = 8)

print(
  DimPlot(seur, group.by = "Phase", label = FALSE, raster = FALSE) +
    labs(title = "Cell Cycle Phase") +
    scale_color_manual(
      values = c("G1" = "#8E44AD",    # Rich Purple
                 "S" = "#E67E22",     # Warm Orange
                 "G2M" = "#3498DB")   # Cool Teal
    ) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 16),
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 16),
      legend.key = element_rect(fill = "white", color = NA)
    )
)

dev.off()

# Create high-resolution PNG version
png(file.path(out_dir, paste0(script_ind, "umap_cell_cycle_phase_highres.png")), 
    width = 12, height = 10, units = "in", res = 300)

print(
  DimPlot(seur, group.by = "Phase", label = FALSE, raster = FALSE) +
    labs(title = "Cell Cycle Phase") +
    scale_color_manual(
      values = c("G1" = "#8E44AD",    # Rich Purple
                 "S" = "#E67E22",     # Warm Orange
                 "G2M" = "#3498DB")   # Cool Teal
    ) +
    theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 22),
      panel.background = element_blank(),
      axis.line = element_line(color = "black", linewidth = 1),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 20),
      legend.key = element_rect(fill = "white", color = NA),
      legend.position = "right",
      plot.margin = margin(20, 20, 20, 20)
    )
)

dev.off()

###########################################################
# Cell Cycle Phase Proportions by Group
###########################################################

message("\n\n          *** Analyzing cell cycle distribution by group... ", Sys.time(),"\n\n")

# Extract phase and group information from metadata
cell_cycle_data <- as.data.frame(seur@meta.data) %>%
  dplyr::select(Phase, group)

# Calculate counts and proportions
phase_proportions <- cell_cycle_data %>%
  group_by(group, Phase) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(group) %>%
  mutate(proportion = count / sum(count),
         percentage = round(proportion * 100, 1))

# Create a summary table
phase_summary <- phase_proportions %>%
  pivot_wider(
    id_cols = group,
    names_from = Phase,
    values_from = c(count, percentage),
    names_glue = "{Phase}_{.value}"
  )

# Save summary tables
write.csv(phase_summary, 
          file.path(out_dir, paste0(script_ind, "cell_cycle_summary_by_group.csv")), 
          row.names = FALSE)
write.csv(phase_proportions, 
          file.path(out_dir, paste0(script_ind, "cell_cycle_proportions_by_group.csv")), 
          row.names = FALSE)

# Print the summary table
print(phase_summary)

###########################################################
# Stacked bar plots of cell cycle distribution
###########################################################

message("\n\n          *** Creating distribution plots... ", Sys.time(),"\n\n")

pdf(file.path(out_dir, paste0(script_ind, "cell_cycle_distribution_by_group.pdf")), 
    width = 10, height = 8)

# Plot 1: Stacked bar plot without percentages
print(
  ggplot(phase_proportions, aes(x = group, y = proportion, fill = Phase)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(
      values = c("G1" = "#8E44AD",    # Rich Purple
                 "S" = "#E67E22",     # Warm Orange
                 "G2M" = "#3498DB")   # Cool Teal
    ) +
    labs(
      title = "Cell Cycle Phase Distribution by Group",
      x = "Group",
      y = "Proportion",
      fill = "Cell Cycle Phase"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 18),
      axis.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 20, face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      legend.position = "right"
    ) +
    scale_y_continuous(labels = scales::percent_format())
)

# Plot 2: Stacked bar plot with percentages displayed on bars
print(
  ggplot(phase_proportions, aes(x = group, y = proportion, fill = Phase)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    geom_text(
      aes(label = paste0(percentage, "%")),
      position = position_stack(vjust = 0.5),
      size = 5, color = "white", fontface = "bold"
    ) +
    scale_fill_manual(
      values = c("G1" = "#8E44AD",    # Rich Purple
                 "S" = "#E67E22",     # Warm Orange
                 "G2M" = "#3498DB")   # Cool Teal
    ) +
    labs(
      title = "Cell Cycle Phase Distribution by Group (with Percentages)",
      x = "Group",
      y = "Proportion",
      fill = "Cell Cycle Phase"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 18),
      axis.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 20, face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      legend.position = "right"
    ) +
    scale_y_continuous(labels = scales::percent_format())
)

dev.off()

###########################################################
# Statistical analysis
###########################################################

message("\n\n          *** Running statistical tests... ", Sys.time(),"\n\n")

# Chi-square test to determine if cell cycle distribution differs between groups
chisq_result <- chisq.test(table(cell_cycle_data$group, cell_cycle_data$Phase))

# Save chi-square results
sink(file.path(out_dir, paste0(script_ind, "cell_cycle_chisq_test.txt")))
cat("Chi-square test for cell cycle distribution across groups:\n")
cat("============================================================\n\n")
print(chisq_result)
cat("\n")
if(chisq_result$p.value < 0.05) {
  cat("Interpretation: The cell cycle distribution differs significantly between groups (p < 0.05)\n")
} else {
  cat("Interpretation: No significant difference in cell cycle distribution between groups (p >= 0.05)\n")
}
sink()

message("\n\nChi-square test results:")
print(chisq_result)

message("\n\n##########################################################################\n",
        "# Completed 09 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")
