message("\n\n##########################################################################\n",
        "# Start 03: Clustering tests and basic cell population analysis ", Sys.time(),
        "\n##########################################################################\n",
        "\n   \n", 
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
library(colorRamps)
library(GenomicRanges)

#specify script/output index as prefix for file names
script_ind = "03_"

#specify output directory
out_dir = paste0(main_dir, "/results/")

#load group and file info
gr_tab = read_csv("results/02_gr_tab_filtered.csv")

#load raw dataset
load(file = paste0(main_dir,"seur_integr_pat.rda")) 

#get marker gene panels
GOI = list()
t1 = read_csv("reference_data/cell_type_markers.csv")
GOI$cell_type_markers = t1$gene[t1$level %in% c("cell_types", "neuronal_lineage")]
GOI$subtype_markers = t1$gene[t1$level %in% c("NPC_patterning","cortical_layers", "Interneuron_subtypes")]


###########################################################
#functions
###########################################################

#custom colour palette for variable values defined in vector v
pal = function(v){
  v2 = length(unique(v))
  if (v2 == 2){
    p2 = c("grey", "blue")
  } else if (v2 ==3){
    p2 = c("blue", "grey", "orange")
  } else if (v2<6){
    p2 = matlab.like(6)[1:v2]
  } else {
    p2 = matlab.like(v2)
  }
  return(p2)
}


#############################################
#Plot UMAP and marker expression (test different clustering resolutions)
#############################################

set.seed(1234)

#define resolutions for clustering tests
test_res = c(0.1, 0.2, 0.25, 0.3)


#subsample seurat object for plotting for large datasets (else plots become too large (pdf with hundreds of MB))

seur0 = seur

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>100000){
  set.seed(1234)
  seur = seur[, sample(cells, size =100000, replace=F)]
} 


#create lists for plot objects
pl_umap = list()
pl_dotplots_cell_type_markers = list()
pl_dotplots_subtype_markers = list()
pl_umap_expr = list()


### generate umap plot by group and sample


#define grouping variables and colors for plotting
gr = unique(gr_tab$group)
samples = unique(gr_tab$sample)


pl_umap[["group"]] = DimPlot(seur, group.by = "group", shuffle = TRUE,
                             label = FALSE, reduction = "umap", 
                             pt.size = 0.01)+
  scale_color_manual(limits = gr, values = pal(gr))+
  labs(title = "UMAP by group")

pl_umap[["sample"]] = DimPlot(seur, group.by = "sample", shuffle = TRUE,
                              label = FALSE, reduction = "umap", 
                              pt.size = 0.01)+
  scale_color_manual(limits = samples, values = pal(samples))+
  labs(title = "UMAP by sample")




### create plots for clusters and marker dot plot for different resolutions (each test_res)
j = test_res[1]
for (j in test_res){
  
  set.seed(1234)
  
  seur0 <- FindClusters(
    object = seur0,
    graph.name = "SCT_snn",
    algorithm = 1,
    resolution = j,
    verbose = FALSE
  )
  
  # Convert cluster labels from 0-based to 1-based
  seur0$seurat_clusters <- as.factor(as.numeric(as.character(seur0$seurat_clusters)) + 1)
  
  #subset for plotting
  
  seur = seur0
  
  cells = colnames(seur@assays$SCT@scale.data)
  
  if (length(cells)>100000){
    set.seed(1234)
    seur = seur[, sample(cells, size =100000, replace=F)]
  } 
  
  cl = unique(seur$seurat_clusters)
  pl_umap[[paste0("resol_", j)]] = DimPlot(seur, group.by = "seurat_clusters", shuffle = TRUE,
                                           label = TRUE, reduction = "umap", pt.size = 0.01)+
    scale_color_manual(limits = cl, values = pal(cl))+
    NoLegend()+labs(title = paste0("clusters resol_", j))
  
  
  #create dotplot with gene expression of marker genes
  
  DefaultAssay(seur0) = "SCT"
  
  pl_dotplots_cell_type_markers[[paste0("resol_", j)]] = 
    DotPlot(seur0, features = intersect(GOI$cell_type_markers, rownames(seur0)), 
            scale.by = "size") +RotatedAxis()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste0("resol_", j))
  
  pl_dotplots_subtype_markers[[paste0("resol_", j)]] = 
    DotPlot(seur0, features = intersect(GOI$subtype_markers, rownames(seur0)), 
            scale.by = "size") +RotatedAxis()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste0("resol_", j))
  
}



seur = seur0

#save dataset with cluster assignment for different cluster resolutions
save(seur, file = paste0(out_dir, script_ind, "seur_integr_filtered.rda")) 


### save umap clustering and dotplots

pdf(file = paste0(out_dir, script_ind, "UMAP_clustering_test_integr_dataset.pdf"), width = 6, height = 5)
lapply(pl_umap, function(x){x})
dev.off()

pdf(file = paste0(out_dir, script_ind, "Dotplots_clustering_test_integr_dataset_cell_type_markers.pdf"), width = 12, height = 10)
lapply(pl_dotplots_cell_type_markers, function(x){x})
dev.off()

pdf(file = paste0(out_dir, script_ind, "Dotplots_clustering_test_integr_dataset_subtype_markers.pdf"), width = 12, height = 10)
lapply(pl_dotplots_subtype_markers, function(x){x})
dev.off()



############################################################
# UMAP plots split by group vs repl/batch_seq
############################################################

seur = seur0

cells = colnames(seur@assays$SCT@scale.data)

if (length(cells)>100000){
  set.seed(1234)
  seur = seur[, sample(cells, size =100000, replace=F)]
} 


t1 = FetchData(seur, vars = c("umap_1", "umap_2", "group", "sample"))


#plot group vs sample

p1 = ggplot(data = t1, aes(x = umap_1, y = umap_2))+geom_point(size = 0.05, alpha = 0.5)+
  ggplot2::facet_wrap(facets = factor(t1$sample, levels = samples), nrow = 2)+
  theme_bw()+
  labs(title = "UMAP by group vs sample")+RestoreLegend()

pdf(file = paste0(out_dir, script_ind, "UMAP_split_by_group_vs_sample.pdf"), width = 10, height = 2.5)
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

pdf(file = paste0(out_dir, script_ind, "marker_expr_UMAP_integr_dataset.pdf"), width = 30, height = 15)
plot(pl_umap_expr)
dev.off()


###################################################################################
# save table for cluster assignment/labelling (with highest resolution) 
###################################################################################

seur = seur0

v1 = unique(seur$seurat_clusters)
v2 = v1[order(v1)]
t1 = tibble(cluster = v2, cluster_name = paste0("s",v2), cell_type = "all", cell_class = "all")

if (file.exists(paste0(out_dir, script_ind, "cluster_assignment.csv"))){
  write_csv(t1, paste0(out_dir, script_ind, "cluster_assignment_new.csv"))
} else {write_csv(t1, paste0(out_dir, script_ind, "cluster_assignment.csv"))}



message("\n\n##########################################################################\n",
        "# Completed 03 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")




