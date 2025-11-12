message("\n\n##########################################################################\n",
        "# Start 04: Clustering tests and basic cell population analysis ", Sys.time(),
        "\n##########################################################################\n",
        "\n   clustering with selected cluster resolution ",
        "\n   differential abundance analysis on cluster level", 
        "\n##########################################################################\n\n")

#set environment and main directory
main_dir = getwd()
setwd(main_dir)

# Open packages necessary for analysis.
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
library(sccomp)
library(cmdstanr)


#specify script/output index as prefix for file names
script_ind = "04_"

#specify output directory
out_dir = paste0(main_dir, "/results/")

#load group and file info
gr_tab = read_csv("results/02_gr_tab_filtered.csv")

#load dataset
load(file = paste0(out_dir, "results/04_seur_integr_labelled.rda")) 

#get marker gene panels
GOI = list()
t1 = read_csv("reference_data/cell_type_markers.csv")
GOI$cell_type_markers = t1$gene[t1$level %in% c("cell_types", "neuronal_lineage")]
GOI$subtype_markers = t1$gene[t1$level %in% c("NPC_patterning","cortical_layers", "I_clusterneuron_subtypes", "Astrocyte_subtypes")]

#define clustering resolution for downstream analyses (for full dataset)
clust_resol = 0.3

#load cluster labels
clust_tab = read_csv(paste0(out_dir, "results/03_cluster_assignment_updated.csv"))


####################################
#Functions
####################################

#custom colour palette for variable values defined in vector v
pal = function(v){
  v2 = length(unique(v))
  if (v2 == 2){
    p2 = c("grey20", "dodgerblue")
  } else if (v2 ==3){
    p2 = c("dodgerblue", "grey20", "orange")
  } else if (v2<6){
    p2 = matlab.like(6)[1:v2]
  } else {
    p2 = matlab.like(v2)
  }
  return(p2)
}

#distinct scale (larger number of colours, colour vector gets shuffled)
pal_dist = function(v){
  v2 = length(unique(v))
  if (v2 < 6){p2 = matlab.like(6)[1:v2]} else {
    p2 = matlab.like(v2)
    set.seed(12)
    p2 = sample(p2)
    }
  return(p2)
}



###########################################################
# save cluster/cell class labels to dataset, plot with cluster/cell_class names
###########################################################

DefaultAssay(seur) = "SCT"
Idents(seur) <- "SCT_snn_res.0.3"

# Add 1 to all cluster numbers in the SCT_snn_res.0.3 column
seur@meta.data$SCT_snn_res.0.3 <- as.factor(as.numeric(as.character(seur@meta.data$SCT_snn_res.0.3)) + 1)

# Set the identities to use this modified column
Idents(seur) <- "SCT_snn_res.0.3"

#add clusters with default resolution
seur$seurat_clusters = seur@meta.data[,paste0("SCT_snn_res.", clust_resol)]

seur$cluster_name = clust_tab$cluster_name[match(seur$seurat_clusters, clust_tab$cluster)]
seur$cell_type = clust_tab$cell_type[match(seur$seurat_clusters, clust_tab$cluster)]
seur$cell_class = clust_tab$cell_class[match(seur$seurat_clusters, clust_tab$cluster)]

save(seur, file = paste0(out_dir,script_ind, "seur_integr_filtered.rda")) 


### plot whole dataset as UMAP with named clusters/cell classes

#define grouping variables
gr = unique(gr_tab$group)
samples = unique(gr_tab$sample)
clusters = clust_tab$cluster
cluster_names = clust_tab$cluster_name
cell_types = unique(clust_tab$cell_type)
cell_classes = unique(clust_tab$cell_class)


#subsample seurat object for plotting for large datasets (else plots become too large (pdf with hundreds of MB))

seur0 = seur

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>100000){
  set.seed(1234)
  seur = seur[, sample(cells, size =100000, replace=F)]
} 


#umap plots for cell cluster and cell class (labelled)
pl = list()

pl[["group"]] = DimPlot(seur, group.by = "group", shuffle = TRUE,
                             label = FALSE, reduction = "umap", 
                             pt.size = 0.01)+
  scale_color_manual(limits = gr, values = pal(gr))+
  labs(title = "UMAP by group")

pl[["cluster"]] = DimPlot(seur, group.by = "seurat_clusters", shuffle = TRUE, repel = TRUE, 
                                    label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = clusters, values = pal(clusters))+
  NoLegend()+labs(title = "clusters")

pl[["cluster_dist"]] = DimPlot(seur, group.by = "seurat_clusters", shuffle = TRUE, repel = TRUE, 
                          label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = clusters, values = pal_dist(clusters))+
  NoLegend()+labs(title = "clusters")

pl[["cluster_name_umap"]] = DimPlot(seur, group.by = "cluster_name", shuffle = TRUE, repel = TRUE, 
                                    label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cluster_names, values = pal(cluster_names))+
  NoLegend()+labs(title = "cluster_names")

pl[["cluster_name_umap_dist"]] = DimPlot(seur, group.by = "cluster_name", shuffle = TRUE, repel = TRUE, 
                                    label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cluster_names, values = pal_dist(cluster_names))+
  NoLegend()+labs(title = "cluster_names")

pl[["cell_type_umap"]] = DimPlot(seur, group.by = "cell_type", shuffle = TRUE, repel = TRUE, 
                                    label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cell_types, values = pal(cell_types))+
  NoLegend()+labs(title = "cell_types")

pl[["cell_type_umap_dist"]] = DimPlot(seur, group.by = "cell_type", shuffle = TRUE, repel = TRUE, 
                                  label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cell_types, values = pal_dist(cell_types))+
  NoLegend()+labs(title = "cell_types")

pl[["cell_class_umap"]] = DimPlot(seur, group.by = "cell_class", shuffle = TRUE, repel = TRUE, 
                                  label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cell_classes, values = pal(cell_classes))+
  NoLegend()+labs(title = "cell_classes")

pl[["cell_class_umap_dist"]] = DimPlot(seur, group.by = "cell_class", shuffle = TRUE, repel = TRUE, 
                                label = TRUE, reduction = "umap", pt.size = 0.01)+
  scale_color_manual(limits = cell_classes, values = pal_dist(cell_classes))+
  NoLegend()+labs(title = "cell_classes")

pdf(file = paste0(out_dir,script_ind, "UMAP_clusters_labelled.pdf"), width = 6, height = 5)
lapply(pl, function(x){x})
dev.off()


#labelled dotplot for marker gene expression by cluster (split by cell type markers vs neuron subtype markers)

seur = seur0

DefaultAssay(seur) = "SCT"

p1 = DotPlot(seur, features = intersect(GOI$cell_type_markers, rownames(seur)), 
             group.by = "cluster_name", scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_discrete(limits = cluster_names)

p2 = DotPlot(seur, features = intersect(GOI$subtype_markers, rownames(seur)), 
             group.by = "cluster_name", scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_y_discrete(limits = cluster_names)

pdf(file = paste0(out_dir,script_ind, "Cell_markers_dotplot_clusters_labelled.pdf"), 
    width = 10, height = 8)
plot(p1)
plot(p2)
dev.off()



#######################################################################################
### add plots split by group vs repl/batch_seq/dev_PCW/stage colored by cluster_name
#######################################################################################

seur = seur0

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>100000){
  set.seed(1234)
  seur = seur[, sample(cells, size =100000, replace=F)]
} 


t1 = FetchData(seur, vars = c("umap_1", "umap_2", "cluster_name","group", "sample"))

# Convert cluster_name to factor
t1$cluster_name <- as.factor(t1$cluster_name)

#plot group vs sample

p1 = ggplot(data = t1, aes(x = umap_1, y = umap_2, color = cluster_name))+geom_point(size = 0.05, alpha = 0.5)+
  scale_color_manual(limits = cluster_names, values = pal_dist(cluster_names))+
  ggplot2::facet_wrap(facets = factor(t1$sample, levels = samples), nrow = 2)+
  theme_bw()+
  labs(title = "UMAP by group vs sample")+RestoreLegend()

pdf(file = paste0(out_dir, script_ind, "UMAP_split_by_group_vs_sample.pdf"), width = 10, height = 5)
plot(p1)
dev.off()




######################################################
# plots for marker expression on UMAP
######################################################

#subsample seurat object for large datasets (else plots become too large (pdf with hundreds of MB))

seur = seur0

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>10000){
  set.seed(1234)
  seur = seur[, sample(cells, size =10000, replace=F)]
} 

DefaultAssay(seur) = "SCT"

pl_umap_expr = FeaturePlot(seur, features = c(GOI$cell_type_markers, GOI$subtype_markers), order = TRUE,
                           reduction = "umap", pt.size = 0.01, ncol = 9)

pdf(file = paste0(out_dir, script_ind, "marker_expr_UMAP_integr_dataset.pdf"), width = 30, height = 30)
plot(pl_umap_expr)
dev.off()


#############################################################################
# cluster quantification and differential abundance analysis
#############################################################################

seur = seur0

meta = seur@meta.data

# create table with one row for each cluster for each sample

stat_tab = tibble(cluster = unlist(lapply(cluster_names, rep, length.out = length(samples))), 
                  sample = rep(samples, length(cluster_names)))
stat_tab$group = gr_tab$group[match(stat_tab$sample, gr_tab$sample)] #add group assignment

#count cells per cluster per sample, add to stat_tab
t1 = meta %>% group_by(cluster_name, sample) %>% summarize(N_cells = n())
stat_tab$N_cells = t1$N_cells[match( paste0(stat_tab$cluster,stat_tab$sample), paste0(t1$cluster_name,t1$sample) )]
stat_tab$N_cells[is.na(stat_tab$N_cells)] = 0

#add total cells per sample and per cluster, fraction of cluster, fraction of sample
N_sample = stat_tab %>% group_by(sample) %>% summarize(N = sum(N_cells))
stat_tab$N_sample = N_sample$N[match(stat_tab$sample, N_sample$sample)]
stat_tab$fract_sample = stat_tab$N_cells / stat_tab$N_sample
N_cluster = stat_tab %>% group_by(cluster) %>% summarize(N = sum(N_cells))
stat_tab$N_cluster = N_cluster$N[match(stat_tab$cluster, N_cluster$cluster)]
stat_tab$fract_cluster = stat_tab$N_cells / stat_tab$N_cluster

#filtered_table <- `N_cluster` %>% filter(cluster %in% 1:15)

write_csv(stat_tab, file = paste0(out_dir,script_ind,"cell_abundance_by_sample_cluster.csv"))


### crossbar-dotplot quantification of cluster contribution fraction of sample

t1 = stat_tab
t2 = t1 %>% group_by(cluster, group) %>% 
  summarise(mean_fract = mean(fract_sample), sd_fract = sd(fract_sample))

p1 = ggplot()+
  geom_col(data = t2, aes(x = cluster, y = mean_fract, color = group), fill = "grey90", 
           position = position_dodge(), width = 0.5, lwd = 0.3)+
  geom_errorbar(data = t2, aes(x = cluster,
                               ymin = mean_fract-sd_fract,
                               y = mean_fract,
                               ymax = mean_fract+sd_fract,
                               color = group),
                position = position_dodge(width = 0.5), width = 0.3, lwd = 0.2)+
  geom_point(data = t1, aes(x = cluster, y = fract_sample, color = group), 
             position = position_dodge(width = 0.5), size = 0.5, stroke = 0.3)+
  geom_hline(yintercept = 0)+
  scale_x_discrete(limits = cluster_names)+
  scale_color_manual(limits = gr, values = pal(gr))+
  scale_fill_manual(limits = gr, values = pal(gr))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file = paste0(out_dir,script_ind,"cell_abundance_by_sample_cluster.pdf"), 
    width = 12, height = 7)
plot(p1)
dev.off()


### heatmap quantification of cluster contribution fraction of sample

heatmap_vir = function(m1, main = "", cluster_rows = FALSE, cluster_cols = FALSE){
  p1 = pheatmap(m1, show_rownames=TRUE, cluster_rows = cluster_rows,
                cluster_cols = cluster_cols, show_colnames = TRUE, 
                clustering_distance_rows = "euclidean",
                clustering_method = "ward.D2",
                treeheight_row = 50,
                color = viridis_pal(option = "magma")(250),
                breaks = seq(0, max(m1), length.out = 251),
                border_color = NA, fontsize = 10,
                cellwidth = 10, cellheight = 10,
                main = main
  )
}

t1 = stat_tab

m1 = matrix(nrow = length(cluster_names), ncol = length(samples), dimnames = list(cluster_names, samples))

for (i in samples){
  t2 = t1[t1$sample==i,]
  m1[,i] = t2$fract_sample[match(rownames(m1),t2$cluster)]
}
m1[is.na(m1)] = 0

pdf(file = paste0(out_dir,script_ind,"cell_abundance_by_sample_cluster_heatmap.pdf"), 
    width = 10, height = 8)
heatmap_vir(m1, main = "cell_abundance_by_area_and_sample")
heatmap_vir(m1, main = "cell_abundance_by_area_and_sample", 
            cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()



### statistical analysis t-tests by cluster

t1 = stat_tab

t2 = tibble(cluster = cluster_names, p = NA, padj = NA)

for (i in cluster_names){
  t4 = t1[t1$cluster == i,]
  t5 = pairwise.t.test(t4$fract_sample, t4$group, p.adjust.method = "none")
  t2$p[t2$cluster == i] = t5$p.value[1,1]
}

t2$padj = p.adjust(t2$p, method = "BH")

write_csv(t2, paste0(out_dir,script_ind,"cell_abundance_by_sample_cluster_t-Test.csv"))

#Setting ctrl_u as reference

seur@meta.data$group <- factor(seur@meta.data$group, levels = c("ctrl_u", "ctrl_SB", "VCP_u", "VCP_SB"))


###########################################################
# sccomp differential cell cluster abundance anlaysis
###########################################################

# Set identities to the existing resolution 0.3 clustering
Idents(seur) <- "SCT_snn_res.0.3"

# Set cluster names properly - since clusters are already 1-16, just format with leading zeros
seur$cluster_name <- sprintf("%02d", as.numeric(Idents(seur)))

print("Final cluster names:")
print(sort(unique(seur$cluster_name)))

# Now continue with sccomp analysis
sccomp_result = 
  seur |>
  sccomp_estimate( 
    formula_composition = ~ group, 
    .sample =  sample, 
    .cell_group = cluster_name, 
    bimodal_mean_variability_association = TRUE,
    cores = 8 
  ) |> 
  sccomp_test()

write_csv(sccomp_result, paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group.csv"))

pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_estimates.pdf"), 
    width = 10, height = 12)
sccomp_result |> 
  plot_1D_intervals() +
  theme(
    text = element_text(size = 20),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    plot.title = element_text(size = 22)
  )
dev.off()

# Create simplified boxplot with significance coloring
pdf(file = paste0(out_dir,script_ind,"sccomp_cell_abundance_by_cluster_group_boxplot.pdf"), 
    width = 20, height = 14)

# Extract proportion data from Seurat object
prop_data <- seur@meta.data %>%
  dplyr::filter(group %in% c("ctrl_u", "ctrl_SB", "VCP_u", "VCP_SB")) %>%  # Only include known groups
  group_by(sample, group, cluster_name) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample, group) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  mutate(
    cluster_name = factor(cluster_name, levels = sprintf("%02d", 1:16)),
    group = factor(group, levels = c("ctrl_u", "ctrl_SB", "VCP_u", "VCP_SB"))
  )

# Process significance data
sig_data <- sccomp_result %>%
  dplyr::filter(parameter != "(Intercept)") %>%
  dplyr::select(cluster_name, parameter, c_FDR) %>%
  mutate(
    significant = c_FDR < 0.05,
    sig_group = case_when(
      parameter == "groupctrl_SB" ~ "ctrl_SB",
      parameter == "groupVCP_SB" ~ "VCP_SB", 
      parameter == "groupVCP_u" ~ "VCP_u",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(significant & !is.na(sig_group)) %>%
  dplyr::select(cluster_name, sig_group) %>%
  mutate(is_significant = TRUE)

# Add significance info to proportion data
prop_data_final <- prop_data %>%
  left_join(sig_data, by = c("cluster_name", "group" = "sig_group")) %>%
  mutate(
    is_significant = ifelse(is.na(is_significant), FALSE, TRUE),
    fill_group = ifelse(is_significant, as.character(group), "Not significant")
  )

# Define colors - only for the 4 known groups
colors <- c(
  "ctrl_u" = "#66C2A5",     
  "ctrl_SB" = "#FC8D62",
  "VCP_u" = "#8DA0CB", 
  "VCP_SB" = "#E78AC3",
  "Not significant" = "white"
)

# Create the plot
prop_data_final %>%
  ggplot(aes(x = group, y = proportion)) +
  geom_boxplot(aes(fill = fill_group), 
               linewidth = 1.2, outlier.size = 3, width = 0.6) +
  geom_point(aes(color = group), size = 3, alpha = 0.6, 
             position = position_jitter(width = 0.15)) +
  facet_wrap(~ cluster_name, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = colors, name = "Significance") +
  scale_color_manual(values = colors[1:4], guide = "none") +  # Remove color guide, only show fill
  scale_x_discrete(limits = c("ctrl_u", "ctrl_SB", "VCP_u", "VCP_SB")) +
  labs(x = "Biological condition", 
       y = "Cell-group proportion",
       title = "Cell Abundance by Cluster and Group",
       subtitle = "Colored boxes indicate significant differences compared to ctrl_u") +
  theme_bw() +
  theme(
    text = element_text(size = 24),
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 26, face = "bold"),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 24, face = "bold"),
    strip.text = element_text(size = 24, face = "bold"),
    plot.title = element_text(size = 28, face = "bold"),
    plot.subtitle = element_text(size = 20),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20),
    legend.position = "bottom",
    plot.margin = margin(25, 25, 25, 25),
    panel.grid.major = element_line(linewidth = 0.8),
    panel.grid.minor = element_line(linewidth = 0.4),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) +
  guides(fill = guide_legend(nrow = 1))

dev.off()

# Export comprehensive statistical information for ALL sccomp comparisons
cat("\n=== Exporting comprehensive sccomp statistical results ===\n")

# 1. Complete sccomp results with additional formatting
complete_sccomp_results <- sccomp_result %>%
  mutate(
    comparison_type = case_when(
      parameter == "groupctrl_SB" ~ "ctrl_SB vs ctrl_u",
      parameter == "groupVCP_u" ~ "VCP_u vs ctrl_u", 
      parameter == "groupVCP_SB" ~ "VCP_SB vs ctrl_u",
      parameter == "(Intercept)" ~ "Intercept",
      TRUE ~ parameter
    ),
    significant_FDR = c_FDR < 0.05,
    effect_direction = case_when(
      c_effect > 0 ~ "increased",
      c_effect < 0 ~ "decreased", 
      TRUE ~ "no_change"
    ),
    effect_magnitude = case_when(
      abs(c_effect) < 0.5 ~ "small",
      abs(c_effect) < 1.0 ~ "medium",
      abs(c_effect) < 2.0 ~ "large",
      TRUE ~ "very_large"
    )
  ) %>%
  arrange(cluster_name, parameter)

# 2. Generate comprehensive descriptive statistics for all groups and clusters
descriptive_stats_all <- prop_data %>%
  group_by(cluster_name, group) %>%
  summarise(
    n_samples = n(),
    mean_proportion = mean(proportion),
    median_proportion = median(proportion),
    sd_proportion = sd(proportion),
    min_proportion = min(proportion),
    max_proportion = max(proportion),
    q25_proportion = quantile(proportion, 0.25),
    q75_proportion = quantile(proportion, 0.75),
    iqr_proportion = IQR(proportion),
    .groups = "drop"
  ) %>%
  arrange(cluster_name, group)

# 3. Create pairwise comparison summary
pairwise_comparisons <- complete_sccomp_results %>%
  dplyr::filter(parameter != "(Intercept)") %>%
  dplyr::select(cluster_name, comparison_type, c_effect, c_lower, c_upper, 
                c_FDR, significant_FDR, effect_direction, effect_magnitude) %>%
  arrange(cluster_name, comparison_type)

# 4. Summary of significant results only
significant_results_summary <- pairwise_comparisons %>%
  dplyr::filter(significant_FDR == TRUE) %>%
  arrange(c_FDR)

# 5. Cluster-wise summary (how many significant comparisons per cluster)
cluster_significance_summary <- pairwise_comparisons %>%
  group_by(cluster_name) %>%
  summarise(
    total_comparisons = n(),
    significant_comparisons = sum(significant_FDR),
    prop_significant = significant_comparisons / total_comparisons,
    min_FDR = min(c_FDR),
    max_effect_size = max(abs(c_effect)),
    .groups = "drop"
  ) %>%
  arrange(desc(significant_comparisons), min_FDR)

# 6. Comparison-wise summary (how many clusters significant per comparison type)
comparison_significance_summary <- pairwise_comparisons %>%
  group_by(comparison_type) %>%
  summarise(
    total_clusters = n(),
    significant_clusters = sum(significant_FDR),
    prop_significant = significant_clusters / total_clusters,
    mean_effect_size = mean(abs(c_effect)),
    median_FDR = median(c_FDR),
    .groups = "drop"
  ) %>%
  arrange(desc(significant_clusters))

# 7. Individual sample data for transparency
individual_sample_data <- prop_data %>%
  dplyr::select(cluster_name, group, sample, proportion) %>%
  arrange(cluster_name, group, sample)

# 8. Wide format descriptive stats for easy comparison
descriptive_stats_wide <- descriptive_stats_all %>%
  pivot_wider(
    names_from = group,
    values_from = c(mean_proportion, median_proportion, sd_proportion, n_samples),
    names_sep = "_"
  ) %>%
  arrange(cluster_name)

# Save all statistical files
write_csv(complete_sccomp_results, paste0(out_dir, script_ind, "sccomp_complete_results_formatted.csv"))
write_csv(descriptive_stats_all, paste0(out_dir, script_ind, "sccomp_descriptive_stats_all_groups.csv"))
write_csv(descriptive_stats_wide, paste0(out_dir, script_ind, "sccomp_descriptive_stats_wide_format.csv"))
write_csv(pairwise_comparisons, paste0(out_dir, script_ind, "sccomp_pairwise_comparisons_summary.csv"))
write_csv(significant_results_summary, paste0(out_dir, script_ind, "sccomp_significant_results_only.csv"))
write_csv(cluster_significance_summary, paste0(out_dir, script_ind, "sccomp_cluster_significance_summary.csv"))
write_csv(comparison_significance_summary, paste0(out_dir, script_ind, "sccomp_comparison_significance_summary.csv"))
write_csv(individual_sample_data, paste0(out_dir, script_ind, "sccomp_individual_sample_proportions.csv"))

# Generate comprehensive analysis report
report_lines <- c(
  "=== SCCOMP DIFFERENTIAL ABUNDANCE ANALYSIS REPORT ===",
  paste("Generated on:", Sys.time()),
  "",
  "ANALYSIS OVERVIEW:",
  paste("- Total clusters analyzed:", length(unique(complete_sccomp_results$cluster_name))),
  paste("- Groups compared:", paste(unique(prop_data$group), collapse = ", ")),
  paste("- Total pairwise comparisons:", nrow(pairwise_comparisons)),
  paste("- Significant comparisons (FDR < 0.05):", nrow(significant_results_summary)),
  paste("- Overall significance rate:", round(nrow(significant_results_summary)/nrow(pairwise_comparisons)*100, 1), "%"),
  "",
  "MOST SIGNIFICANT RESULTS (Top 10):",
  paste("Rank | Cluster | Comparison | Effect | FDR | Direction")
)

if(nrow(significant_results_summary) > 0) {
  top_results <- significant_results_summary %>% 
    slice_head(n = 10) %>%
    mutate(rank = row_number())
  
  for(i in 1:nrow(top_results)) {
    report_lines <- c(report_lines,
                      paste(sprintf("%2d", top_results$rank[i]), "|", 
                            top_results$cluster_name[i], "|",
                            top_results$comparison_type[i], "|",
                            sprintf("%.3f", top_results$c_effect[i]), "|",
                            sprintf("%.2e", top_results$c_FDR[i]), "|",
                            top_results$effect_direction[i])
    )
  }
} else {
  report_lines <- c(report_lines, "No significant results found.")
}

report_lines <- c(report_lines,
                  "",
                  "CLUSTER SUMMARY (Most affected clusters):",
                  "Cluster | Sig.Comparisons | Total | Min.FDR | Max.Effect"
)

top_clusters <- cluster_significance_summary %>% slice_head(n = 10)
for(i in 1:nrow(top_clusters)) {
  report_lines <- c(report_lines,
                    paste(top_clusters$cluster_name[i], "|",
                          top_clusters$significant_comparisons[i], "|",
                          top_clusters$total_comparisons[i], "|", 
                          sprintf("%.2e", top_clusters$min_FDR[i]), "|",
                          sprintf("%.3f", top_clusters$max_effect_size[i]))
  )
}

report_lines <- c(report_lines,
                  "",
                  "COMPARISON SUMMARY:",
                  "Comparison | Sig.Clusters | Total | Mean.Effect | Median.FDR"
)

for(i in 1:nrow(comparison_significance_summary)) {
  report_lines <- c(report_lines,
                    paste(comparison_significance_summary$comparison_type[i], "|",
                          comparison_significance_summary$significant_clusters[i], "|",
                          comparison_significance_summary$total_clusters[i], "|",
                          sprintf("%.3f", comparison_significance_summary$mean_effect_size[i]), "|",
                          sprintf("%.2e", comparison_significance_summary$median_FDR[i]))
  )
}

report_lines <- c(report_lines,
                  "",
                  "FILES GENERATED:",
                  "1. sccomp_complete_results_formatted.csv - Complete sccomp results with formatting",
                  "2. sccomp_descriptive_stats_all_groups.csv - Descriptive statistics by cluster and group", 
                  "3. sccomp_descriptive_stats_wide_format.csv - Wide format descriptive statistics",
                  "4. sccomp_pairwise_comparisons_summary.csv - Summary of all pairwise comparisons",
                  "5. sccomp_significant_results_only.csv - Only significant results",
                  "6. sccomp_cluster_significance_summary.csv - Per-cluster significance summary",
                  "7. sccomp_comparison_significance_summary.csv - Per-comparison type summary", 
                  "8. sccomp_individual_sample_proportions.csv - Individual sample data",
                  "9. sccomp_comprehensive_analysis_report.txt - This report"
)

# Save the comprehensive report
writeLines(report_lines, paste0(out_dir, script_ind, "sccomp_comprehensive_analysis_report.txt"))

# Print summary to console
cat("Statistical export complete! Files saved:\n")
cat("1. Complete results:", paste0(out_dir, script_ind, "sccomp_complete_results_formatted.csv"), "\n")
cat("2. Descriptive stats:", paste0(out_dir, script_ind, "sccomp_descriptive_stats_all_groups.csv"), "\n") 
cat("3. Wide format stats:", paste0(out_dir, script_ind, "sccomp_descriptive_stats_wide_format.csv"), "\n")
cat("4. Pairwise comparisons:", paste0(out_dir, script_ind, "sccomp_pairwise_comparisons_summary.csv"), "\n")
cat("5. Significant only:", paste0(out_dir, script_ind, "sccomp_significant_results_only.csv"), "\n")
cat("6. Cluster summary:", paste0(out_dir, script_ind, "sccomp_cluster_significance_summary.csv"), "\n")
cat("7. Comparison summary:", paste0(out_dir, script_ind, "sccomp_comparison_significance_summary.csv"), "\n")
cat("8. Individual samples:", paste0(out_dir, script_ind, "sccomp_individual_sample_proportions.csv"), "\n")
cat("9. Analysis report:", paste0(out_dir, script_ind, "sccomp_comprehensive_analysis_report.txt"), "\n")

cat("\nQuick overview:\n")
cat("- Total comparisons:", nrow(pairwise_comparisons), "\n")
cat("- Significant comparisons:", nrow(significant_results_summary), "\n")
cat("- Most affected cluster:", cluster_significance_summary$cluster_name[1], 
    "(", cluster_significance_summary$significant_comparisons[1], " significant comparisons)\n")


######## Create a new plot showing only significant clusters (excluding 02, 12, 14)
pdf(file = paste0(out_dir,script_ind,"sccomp_significant_clusters_only.pdf"), 
    width = 20, height = 8)  # Wider and shorter

# Extract proportion data from Seurat object with complete combinations
prop_data_raw <- seur@meta.data %>%
  dplyr::filter(group %in% c("ctrl_u", "ctrl_SB", "VCP_u", "VCP_SB")) %>%
  group_by(sample, group, cluster_name) %>%
  summarise(count = n(), .groups = "drop")

# Get the actual sample-group relationships
sample_groups <- seur@meta.data %>% 
  dplyr::select(sample, group) %>% 
  distinct() %>%
  dplyr::filter(group %in% c("ctrl_u", "ctrl_SB", "VCP_u", "VCP_SB"))

# Create complete combinations only for valid sample-group pairs
all_combinations <- sample_groups %>%
  crossing(cluster_name = sprintf("%02d", 1:16))

# Join with actual counts and fill missing with 0
prop_data <- all_combinations %>%
  left_join(prop_data_raw, by = c("sample", "group", "cluster_name")) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  group_by(sample, group) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  mutate(
    cluster_name = factor(cluster_name, levels = sprintf("%02d", 1:16)),
    group = factor(group, levels = c("ctrl_u", "ctrl_SB", "VCP_u", "VCP_SB"))
  )

# Process significance data and identify significant clusters
sig_data <- sccomp_result %>%
  dplyr::filter(parameter != "(Intercept)") %>%
  dplyr::select(cluster_name, parameter, c_FDR) %>%
  mutate(
    significant = c_FDR < 0.05,
    sig_group = case_when(
      parameter == "groupctrl_SB" ~ "ctrl_SB",
      parameter == "groupVCP_SB" ~ "VCP_SB", 
      parameter == "groupVCP_u" ~ "VCP_u",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(significant & !is.na(sig_group))

# Get list of clusters with ANY significant difference, excluding specified clusters
significant_clusters <- unique(sig_data$cluster_name)
excluded_clusters <- c("02", "12", "14")
significant_clusters_filtered <- significant_clusters[!significant_clusters %in% excluded_clusters]

print("Clusters with significant differences (after exclusions):")
print(significant_clusters_filtered)

# Filter data to only include significant clusters (minus excluded ones)
prop_data_sig_only <- prop_data %>%
  dplyr::filter(cluster_name %in% significant_clusters_filtered)

# Add significance info
sig_info <- sig_data %>%
  dplyr::filter(cluster_name %in% significant_clusters_filtered) %>%
  dplyr::select(cluster_name, sig_group) %>%
  mutate(is_significant = TRUE)

prop_data_final <- prop_data_sig_only %>%
  left_join(sig_info, by = c("cluster_name", "group" = "sig_group")) %>%
  mutate(
    is_significant = ifelse(is.na(is_significant), FALSE, TRUE),
    fill_group = ifelse(is_significant, as.character(group), "Not significant")
  )

# Define colors with better contrast
colors <- c(
  "ctrl_u" = "#66C2A5",     
  "ctrl_SB" = "#FC8D62",
  "VCP_u" = "#8DA0CB", 
  "VCP_SB" = "#E78AC3",
  "Not significant" = "white"
)

# Create the enhanced plot
prop_data_final %>%
  ggplot(aes(x = group, y = proportion)) +
  geom_boxplot(aes(fill = fill_group), 
               linewidth = 1.2, outlier.size = 3, width = 0.7, alpha = 0.8) +
  geom_point(aes(fill = group), size = 3.5, alpha = 0.9, 
             position = position_jitter(width = 0.2, seed = 42),
             shape = 21, color = "black", stroke = 1.2) +
  facet_wrap(~ cluster_name, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = colors, name = "Significance") +
  scale_x_discrete(limits = c("ctrl_u", "ctrl_SB", "VCP_u", "VCP_SB")) +
  labs(x = "Biological condition", 
       y = "Proportion") +
  theme_minimal() +
  theme(
    text = element_text(size = 22, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    axis.text.x = element_text(size = 20, color = "black", face = "bold", angle = 45, hjust = 1),
    axis.title = element_text(size = 24, face = "bold", color = "black"),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.title.x = element_text(margin = margin(t = 15)),
    legend.text = element_text(size = 20, color = "black"),
    legend.title = element_text(size = 22, face = "bold", color = "black"),
    strip.text = element_text(size = 22, face = "bold", color = "black",
                              margin = margin(5, 5, 5, 5)),
    legend.position = "bottom",
    plot.margin = margin(20, 20, 20, 20),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(nrow = 1, 
                             override.aes = list(color = "black", stroke = 1)))

dev.off()

# Generate and save comprehensive statistical information for selected clusters comparison


# Create a plot comparing only clusters 01, 05, 07, 10 between ctrl_u and VCP_u
pdf(file = paste0(out_dir,script_ind,"sccomp_selected_clusters_ctrl_VCP_comparison.pdf"), 
    width = 12, height = 8)

# Define the clusters of interest (with proper formatting)
clusters_of_interest <- c("01", "05", "07", "10")

# Filter data for only ctrl_u and VCP_u conditions and selected clusters
prop_data_selected <- prop_data %>%
  dplyr::filter(
    group %in% c("ctrl_u", "VCP_u"),
    cluster_name %in% clusters_of_interest
  ) %>%
  mutate(
    cluster_name = factor(cluster_name, levels = clusters_of_interest),
    group = factor(group, levels = c("ctrl_u", "VCP_u"))
  )

# Get significance data for these specific clusters and comparison
sig_data_selected <- sccomp_result %>%
  dplyr::filter(
    parameter == "groupVCP_u",  # This compares VCP_u vs ctrl_u (reference)
    cluster_name %in% clusters_of_interest
  ) %>%
  dplyr::select(cluster_name, c_FDR, c_effect) %>%
  mutate(
    significant = c_FDR < 0.05,
    direction = ifelse(c_effect > 0, "increased", "decreased")
  )

# Add significance info to proportion data
prop_data_selected_final <- prop_data_selected %>%
  left_join(sig_data_selected, by = "cluster_name") %>%
  mutate(
    is_significant = ifelse(is.na(significant), FALSE, significant),
    fill_group = case_when(
      is_significant & group == "VCP_u" ~ paste0("VCP_u (", direction, ")"),
      group == "ctrl_u" ~ "ctrl_u",
      TRUE ~ "VCP_u (not significant)"
    )
  )

# Define colors for this specific comparison
comparison_colors <- c(
  "ctrl_u" = "white",
  "VCP_u (increased)" = "#8DA0CB",
  "VCP_u (decreased)" = "#8DA0CB", 
  "VCP_u (not significant)" = "lightgrey"
)

# Create the comparison plot
comparison_plot <- prop_data_selected_final %>%
  ggplot(aes(x = group, y = proportion)) +
  geom_boxplot(aes(fill = fill_group), 
               linewidth = 1.2, outlier.size = 3, width = 0.6, alpha = 0.8) +
  geom_point(aes(fill = fill_group), size = 3.5, alpha = 0.9, 
             position = position_jitter(width = 0.15, seed = 42),
             shape = 21, color = "black", stroke = 1.2) +
  facet_wrap(~ paste("Cluster", cluster_name), scales = "free_y", ncol = 2) +
  scale_fill_manual(values = comparison_colors, name = "Group & Significance") +
  scale_x_discrete(labels = c("ctrl_u" = "Control", "VCP_u" = "VCP")) +
  labs(
    x = "Condition", 
    y = "Proportion",
    title = "Cell Abundance Comparison: Control vs VCP",
    subtitle = "Selected clusters (01, 05, 07, 10) - Untreated conditions only"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 18, color = "black", face = "bold"),
    axis.title = element_text(size = 22, face = "bold", color = "black"),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.title.x = element_text(margin = margin(t = 15)),
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 18, face = "bold", color = "black"),
    strip.text = element_text(size = 18, face = "bold", color = "black",
                              margin = margin(8, 8, 8, 8)),
    plot.title = element_text(size = 24, face = "bold", color = "black", hjust = 0.5),
    plot.subtitle = element_text(size = 16, color = "black", hjust = 0.5),
    legend.position = "bottom",
    plot.margin = margin(20, 20, 20, 20),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(nrow = 2, 
                             override.aes = list(color = "black", stroke = 1)))

print(comparison_plot)
dev.off()

# Generate and save summary statistics
summary_stats <- prop_data_selected %>%
  group_by(cluster_name, group) %>%
  summarise(
    mean_prop = mean(proportion),
    median_prop = median(proportion),
    sd_prop = sd(proportion),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(cluster_name, group)

# Generate significance results
significance_summary <- sig_data_selected %>%
  dplyr::select(cluster_name, c_effect, c_FDR, significant, direction) %>%
  arrange(cluster_name)

# Combine both summaries into one comprehensive table
combined_summary <- summary_stats %>%
  left_join(
    significance_summary %>% dplyr::select(cluster_name, c_effect, c_FDR, significant, direction), 
    by = "cluster_name"
  ) %>%
  mutate(
    # Add significance info only for VCP_u rows
    c_effect = ifelse(group == "VCP_u", c_effect, NA),
    c_FDR = ifelse(group == "VCP_u", c_FDR, NA),
    significant = ifelse(group == "VCP_u", significant, NA),
    direction = ifelse(group == "VCP_u", direction, NA)
  )

# Save the combined summary to CSV
write_csv(combined_summary, paste0(out_dir, script_ind, "selected_clusters_ctrl_VCP_summary_stats.csv"))

# Also save just the significance results separately
write_csv(significance_summary, paste0(out_dir, script_ind, "selected_clusters_ctrl_VCP_significance.csv"))

# Print to console as well
cat("\n=== Summary Statistics for Selected Clusters (ctrl_u vs VCP_u) ===\n")
print(summary_stats)

cat("\n=== Significance Results (VCP_u vs ctrl_u) ===\n")
print(significance_summary)

cat("\n=== Files saved ===\n")
cat("Combined summary:", paste0(out_dir, script_ind, "selected_clusters_ctrl_VCP_summary_stats.csv"), "\n")
cat("Significance only:", paste0(out_dir, script_ind, "selected_clusters_ctrl_VCP_significance.csv"), "\n")




message("\n\n##########################################################################\n",
        "# Completed 04 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


