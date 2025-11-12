message("\n\n##########################################################################\n",
        "# Start 07: DEG characterisation ", Sys.time(),
        "\n##########################################################################\n",
        "\n   plot DEG stats per cluster, GO/MSigDB analysis per cluster up vs downreg genes",
        "\n##########################################################################\n\n")

#set environment and main directory
main_dir = getwd()
setwd(main_dir)


# Open packages necessary for analysis.
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

#specify script/output index as prefix for file names
script_ind = "07_"

#specify output directory
out_dir = paste0(main_dir, "/results/")

#load group and file info
gr_tab = read_csv("results/02_gr_tab_filtered.csv")

#load DEseq2 dataset, get genes for module clustering analysis
load(file = paste0(out_dir,"results/06_bulk_data_with_DESeq_results.rda")) 


#get marker gene panels
GOI = list()
t1 = read_csv(paste0(main_dir,"reference_data/transcription_factors.csv"))
GOI$TF = t1$Symbol


#define group comparisons for each cluster
comp_groups = list(VCP.u_vs_ctrl.u = c("VCP_u","ctrl_u"),
                   VCP.SB_vs_VCP.u = c("VCP_SB", "VCP_u"),
                   VCP.SB_vs_ctrl.u = c("VCP_SB", "ctrl_u"))




###########################################################
# functions
###########################################################

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


###function: plot bulk gene expression heatmap with annotations
pl_meta = bulk_data$meta

bulkdata_heatmap = function(pl_mat, pl_meta, x_col, pl_genes = NULL, 
                            meta_annot_cols = NULL,
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                            pl_lims = NULL,  cellwidth = 15, cellheight = 10, 
                            fontsize = 10, title = "Z-score vst-norm gene expression"){
  
  if (is.null(pl_genes)){pl_genes = rownames(pl_mat)}
  
  pl_mat = pl_mat[match(pl_genes, rownames(pl_mat), nomatch = 0),]
  pl_meta = pl_meta[match(colnames(pl_mat), pl_meta[[x_col]]),]
  
  #define annotation bars
  
  if (!is.null(meta_annot_cols)){
    
    annot_col = data.frame(row.names = pl_meta[[x_col]]) #if not defined cbind converts factor values to factor levels
    
    for (col1 in meta_annot_cols){
      v1 = pl_meta[match(colnames(pl_mat), pl_meta[[x_col]]),][[col1]]
      v1 = factor(v1, levels = unique(pl_meta[[col1]]))
      annot_col = as.data.frame(cbind(annot_col, v1))
    }
    colnames(annot_col) = meta_annot_cols
    rownames(annot_col) = colnames(pl_mat)
    
    annot_colors = lapply(meta_annot_cols, function(x){
      v1 = pal(unique(annot_col[[x]]) )
      names(v1) = levels(annot_col[[x]])
      return(v1)
    })
    names(annot_colors) = meta_annot_cols
    
  }else{
    annot_col = NULL
    annot_colors = NULL
  }
  
  # define plot limits
  
  if (is.null(pl_lims)){pl_lims = c(-0.7*max(abs(na.omit(pl_mat))), 0.7*max(abs(na.omit(pl_mat))))}
  
  #create plot
  
  p1 = pheatmap(pl_mat, cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                color = color,
                breaks = seq(pl_lims[1], pl_lims[2], length.out = length(color)+1),
                annotation_col = annot_col, annotation_colors = annot_colors,
                border_color = NA, cellwidth = cellwidth, cellheight = cellheight, 
                fontsize = fontsize, main = title
  )
  
  return(p1)
  
}




#################################################
# Plot number of DEG by cluster by comparison
#################################################

l1 = bulk_data$DEGs
meta = bulk_data$meta

cluster_names = unique(meta$cluster_name)

t1 = tibble(comp = names(l1), comp_group = NA, cluster_name = NA, up_down = NA,
            N_genes = lengths(l1))

for (cl in cluster_names){
  t1$cluster_name[grepl(cl, t1$comp)] = cl
}

for (comp_group in names(comp_groups)){
  t1$comp_group[grepl(comp_group, t1$comp)] = comp_group
}

t1$up_down[grepl("up", t1$comp)] = "up"
t1$up_down[grepl("down", t1$comp)] = "down"

DEGs_by_comp_cl = t1


pl = list()

for (comp_group in names(comp_groups)){
  
  t1 = DEGs_by_comp_cl[DEGs_by_comp_cl$comp_group == comp_group,]
  
  pl[[comp_group]] = ggplot(t1, aes(x = cluster_name, group = up_down))+
    geom_col(aes(y = N_genes, fill = up_down), width = 0.7, position = position_dodge(width = 0.7))+
    scale_x_discrete(limits = cluster_names)+
    scale_fill_manual(limits = c("down", "up"), values = c("blue", "red"))+
    scale_color_manual(limits = c("down", "up"), values = c("blue", "red"))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    labs(title = paste0(comp_groups[[comp_group]][1], " vs ", comp_groups[[comp_group]][2], 
                        " - N DEGs by cluster"))
}


pdf(file = paste0(out_dir,script_ind,"DEGs_numbers_by_comp_and_cluster_deg_characterization.pdf"), 
    width = 5, height = 4)
{
  lapply(pl, function(x){x})
}
dev.off()




#################################################
### GO-BP over-representation analysis of DEGs by cluster combined => mapping to clusters (for each comp_group)
#################################################

message("\n\n          *** GO analysis DEGs vs clusters ... ", Sys.time(),"\n\n")

### extract GO terms

for (comp_group in names(comp_groups)){
  
  message("               Extracting GO terms for ", comp_group)
  
  l1 = bulk_data$DEGs[grepl(comp_group, names(bulk_data$DEGs))]
  comp_DEGs = unique(unlist(l1))
  
  ego = enrichGO(gene         = comp_DEGs,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
  
  
  GO_results_tab = ego@result[ego@result$p.adjust<=0.05,]
  bulk_data$GO_results$full[[comp_group]] = GO_results_tab
  
  # simpify GO results (usually helpful, sometimes removes also many intersting terms) 
  ego_simplified =  dropGO(ego, level = c(1,2)) #drops vast majority of terms in some cases
  ego_simplified = simplify(ego_simplified) #drops vast majority of terms in some cases
  
  GO_results_simplified_tab = ego_simplified@result
  bulk_data$GO_results$simplified[[comp_group]] = GO_results_simplified_tab
  
}


### identify overlap of GO genes with cluster DEGs (full GO terms)

for (comp_group in names(comp_groups)){
  
  comp_DEGs_list = bulk_data$DEGs[grepl(comp_group, names(bulk_data$DEGs))]
  
  GO_results_tab = bulk_data$GO_results$full[[comp_group]]
  v1 = GO_results_tab$geneID
  l1 = str_split(v1, "/")
  names(l1) = GO_results_tab$ID
  
  GO_gene_list = l1
  
  for (comp in names(comp_DEGs_list)){
    
    GO_results_tab[[comp]] = 0
    
    comp_DEGs = comp_DEGs_list[[comp]]
    
    for (go in names(GO_gene_list)){
      
      go_genes = GO_gene_list[[go]]
      GO_DEGs_overlap = length(intersect(comp_DEGs, go_genes))
      GO_results_tab[go,comp] = GO_DEGs_overlap
    }
    
  }
  
  bulk_data$GO_results$full[[comp_group]] = GO_results_tab
  
}



### identify overlap of GO genes with cluster DEGs (simplified GO terms)

for (comp_group in names(comp_groups)){
  
  comp_DEGs_list = bulk_data$DEGs[grepl(comp_group, names(bulk_data$DEGs))]
  
  GO_results_tab = bulk_data$GO_results$simplified[[comp_group]]
  v1 = GO_results_tab$geneID
  l1 = str_split(v1, "/")
  names(l1) = GO_results_tab$ID
  
  GO_gene_list = l1
  
  for (comp in names(comp_DEGs_list)){
    
    GO_results_tab[[comp]] = 0
    
    comp_DEGs = comp_DEGs_list[[comp]]
    
    for (go in names(GO_gene_list)){
      
      go_genes = GO_gene_list[[go]]
      GO_DEGs_overlap = length(intersect(comp_DEGs, go_genes))
      GO_results_tab[go,comp] = GO_DEGs_overlap
    }
    
  }
  
  bulk_data$GO_results$simplified[[comp_group]] = GO_results_tab
  
}


### save overlap of GO genes with cluster DEGs

GO_dir = paste0(out_dir,script_ind,"GO_results/")
if (!dir.exists(GO_dir)){dir.create(GO_dir, recursive = TRUE)}

for (comp_group in names(comp_groups)){
  
  GO_results_tab = bulk_data$GO_results$full[[comp_group]]
  
  write_csv(GO_results_tab, 
            paste0(GO_dir, "GO_results_full_",comp_group,".csv"))
  
  GO_results_tab = bulk_data$GO_results$simplified[[comp_group]]
  
  write_csv(GO_results_tab, 
            paste0(GO_dir, "GO_results_simplified_",comp_group,".csv"))
  
}





### plot number of genes of GO term reg in each cluster (full GO terms for each comp_group)

pdf(file = paste0(out_dir,script_ind, "GO_terms_top20_vs_clusters_by_comp_group.pdf"), 
    width = 18, height = 14)
{
  
  for (comp_group in names(comp_groups)){
    
    t1 = bulk_data$GO_results$full[[comp_group]]
    t2 = t1[,grepl(comp_group, names(t1))] #get columns with overlap gene numbers
    m1 = as.matrix(t2)
    rownames(m1) = paste0(t1$Description, " (", t1$ID, ")")
    
    if(nrow(m1)>20){m1 = m1[1:20,]}
    
    top_go_mat = m1
    
    pheatmap(top_go_mat, show_rownames=TRUE, show_colnames = TRUE,
             cluster_rows = TRUE, cluster_cols = FALSE,  
             clustering_distance_rows = "euclidean",
             clustering_method = "ward.D2",
             treeheight_row = 10, treeheight_col = 10,
             color = colorRampPalette(c("white", "blue"))(250),
             breaks = seq(0, max(top_go_mat), length.out = 251),
             border_color = NA, fontsize = 10,
             cellwidth = 10, cellheight = 10,
             main = paste0(comp_group, "Top20 GO terms vs DEGs by cluster")
    )
    
  }
}

dev.off()
####

pdf(file = paste0(out_dir,script_ind, "GO_terms_top20_vs_clusters_by_comp_group.pdf"), 
    width = 18, height = 14)
{
  
  for (comp_group in names(comp_groups)){
    
    t1 = bulk_data$GO_results$full[[comp_group]]
    t2 = t1[,grepl(comp_group, names(t1))]
    
    # Skip if no matching columns or empty data
    if(ncol(t2) == 0 || nrow(t2) == 0) {
      cat("Skipping", comp_group, "- no data\n")
      next
    }
    
    m1 = as.matrix(t2)
    rownames(m1) = paste0(t1$Description, " (", t1$ID, ")")
    
    if(nrow(m1)>20){m1 = m1[1:20,]}
    
    top_go_mat = m1
    
    # Skip if matrix is too small for meaningful visualization
    if(nrow(top_go_mat) < 1 || ncol(top_go_mat) < 1) {
      cat("Skipping", comp_group, "- matrix too small\n")
      next
    }
    
    pheatmap(top_go_mat, show_rownames=TRUE, show_colnames = TRUE,
             cluster_rows = (nrow(top_go_mat) >= 2), 
             cluster_cols = FALSE,  
             clustering_distance_rows = "euclidean",
             clustering_method = "ward.D2",
             treeheight_row = 10, treeheight_col = 10,
             color = colorRampPalette(c("white", "blue"))(250),
             breaks = seq(0, max(top_go_mat), length.out = 251),
             border_color = NA, fontsize = 10,
             cellwidth = 10, cellheight = 10,
             main = paste0(comp_group, " Top20 GO terms vs DEGs by cluster")
    )
    
  }
}
dev.off()

####

### plot number of genes of GO term reg in each cluster (simplified GO terms for each comp_group)

pdf(file = paste0(out_dir,script_ind, "GO_terms_simplified_top20_vs_clusters_by_comp_group.pdf"), 
    width = 18, height = 14)
{
  
  for (comp_group in names(comp_groups)){
    
    t1 = bulk_data$GO_results$simplified[[comp_group]]
    t2 = t1[,grepl(comp_group, names(t1))] #get columns with overlap gene numbers
    m1 = as.matrix(t2)
    rownames(m1) = paste0(t1$Description, " (", t1$ID, ")")
    
    if(nrow(m1)>20){m1 = m1[1:20,]}
    
    top_go_mat = m1
    
    if (nrow(top_go_mat)>1){
      
      pheatmap(top_go_mat, show_rownames=TRUE, show_colnames = TRUE,
               cluster_rows = TRUE, cluster_cols = FALSE,  
               clustering_distance_rows = "euclidean",
               clustering_method = "ward.D2",
               treeheight_row = 10, treeheight_col = 10,
               color = colorRampPalette(c("white", "blue"))(250),
               breaks = seq(0, max(top_go_mat), length.out = 251),
               border_color = NA, fontsize = 10,
               cellwidth = 10, cellheight = 10,
               main = paste0(comp_group, " - Top20 GO terms vs DEGs by cluster")
      )
    }
  }
}

dev.off()



#save bulk dataset with GO analysis

save(bulk_data, file = paste0(out_dir, script_ind, "bulk_data_w_expr_z_scores_deg_characterization.rda"))




#################################################
# MSigDB over-representation analysis of DEGs combined => mapping to clusters
#################################################

message("\n\n          *** MSigDB overrepresentation analysis DEGs vs clusters... ", Sys.time(),"\n\n")

### get MSigDB collections of interest

msigdb_tab = msigdbr(species = "human")

msigdb_tab_sel = msigdb_tab[msigdb_tab$gs_cat %in% c("H", "C8")|
                              msigdb_tab$gs_subcat %in% c("CGP", "CP:REACTOME", "TFT:GTRD",
                                                          "GO:BP", "GO:CC", "GO:MF", "HPO"),]

msigdbr_t2g = msigdb_tab_sel %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()


### get DEG enrichment of MSigDB sets

for (comp_group in names(comp_groups)){
  
  message("               Extracting MSigDB set overlap for ", comp_group)
  
  l1 = bulk_data$DEGs[grepl(comp_group, names(bulk_data$DEGs))]
  comp_DEGs = unique(unlist(l1))
  
  enrich_res = enricher(gene = comp_DEGs, TERM2GENE = msigdbr_t2g)
  t2 = as.data.frame(enrich_res)
  t2 = t2[t2$p.adjust <= 0.05,]
  t2$gs_cat = msigdb_tab$gs_cat[match(t2$ID, msigdb_tab$gs_name)]
  t2$gs_subcat = msigdb_tab$gs_subcat[match(t2$ID, msigdb_tab$gs_name)]
  t2 = t2[order(t2$gs_cat, t2$gs_subcat),]
  
  msigdb_results_tab = t2
  
  bulk_data$msigdb_results[[comp_group]] = msigdb_results_tab
  
}


### identify overlap of MSigDB set genes with cluster DEGs

for (comp_group in names(comp_groups)){
  
  comp_DEGs_list = bulk_data$DEGs[grepl(comp_group, names(bulk_data$DEGs))]
  
  msigdb_results_tab = bulk_data$msigdb_results[[comp_group]]
  v1 = msigdb_results_tab$geneID
  l1 = str_split(v1, "/")
  names(l1) = msigdb_results_tab$ID
  
  GO_gene_list = l1
  
  for (comp in names(comp_DEGs_list)){
    
    msigdb_results_tab[[comp]] = 0
    
    comp_DEGs = comp_DEGs_list[[comp]]
    
    for (go in names(GO_gene_list)){
      
      go_genes = GO_gene_list[[go]]
      GO_DEGs_overlap = length(intersect(comp_DEGs, go_genes))
      msigdb_results_tab[go,comp] = GO_DEGs_overlap
    }
    
  }
  
  bulk_data$msigdb_results[[comp_group]] = msigdb_results_tab
  
}



### save overlap of MSigDB genes with cluster DEGs

GO_dir = paste0(out_dir,script_ind,"GO_results/")
if (!dir.exists(GO_dir)){dir.create(GO_dir, recursive = TRUE)}

for (comp_group in names(comp_groups)){
  
  GO_results_tab = bulk_data$msigdb_results[[comp_group]]
  
  write_csv(GO_results_tab, 
            paste0(GO_dir, "MSigDB_results_",comp_group,".csv"))
  
}



### plot number of genes of MSigDB set reg in each cluster (for each category and comp_group)

pdf(file = paste0(out_dir,script_ind, "MSigDB_terms_top20_vs_clusters_by_cat_comp_group.pdf"), 
    width = 18, height = 14)
{
  for (comp_group in names(comp_groups)){
    
    t1 = bulk_data$msigdb_results[[comp_group]]
    t1$gs_cat_comb = paste0(t1$gs_cat, "_", t1$gs_subcat)
    
    for (cat in unique(t1$gs_cat_comb)){
      
      t2 = t1[t1$gs_cat_comb == cat,]
      
      t3 = t2[,grepl(comp_group, names(t2))] #get columns with overlap gene numbers
      m1 = as.matrix(t3)
      rownames(m1) = paste0(t2$Description)
      
      if(nrow(m1)>20){m1 = m1[1:20,]}
      
      top_go_mat = m1
      
      if (nrow(top_go_mat)>1){
        
        pheatmap(top_go_mat, show_rownames=TRUE, show_colnames = TRUE,
                 cluster_rows = TRUE, cluster_cols = FALSE,  
                 clustering_distance_rows = "euclidean",
                 clustering_method = "ward.D2",
                 treeheight_row = 10, treeheight_col = 10,
                 color = colorRampPalette(c("white", "blue"))(250),
                 breaks = seq(0, max(top_go_mat), length.out = 251),
                 border_color = NA, fontsize = 10,
                 cellwidth = 10, cellheight = 10,
                 main = paste0(comp_group, " - ", cat," - Top20 MSigDB sets vs DEGs by cluster")
        )
      }
    }
  }
}

dev.off()


#######################################
# calculate per gene Z-scores, mean Z-scores per cluster_name for CON group 
#        and group difference (DELTA) for each comparison, add to bulk_data
#######################################

message("\n\n          *** Calculate gene expression Z-scores... ", Sys.time(),"\n\n")

#calculate Z-score per gene by pseudobulk (cluster_sample)
cluster_sample_mat = t(apply(bulk_data$vst_mat, 1, scale))
colnames(cluster_sample_mat) = colnames(bulk_data$vst_mat)

bulk_data$gene_Z_scores$by_cluster_sample = cluster_sample_mat


#calculate mean Z-score per gene by per cluster for control samples

cluster_gene_mat = matrix(nrow = nrow(bulk_data$vst_mat),
                          ncol = length(cluster_names), 
                          dimnames = list(rownames(bulk_data$vst_mat),
                                          cluster_names))

for (cl in cluster_names){
  cl_bulks = meta$cluster_sample[meta$cluster_name == cl &
                                   meta$group == "PN"] 
  
  cluster_gene_mat[,cl] = apply(cluster_sample_mat[,cl_bulks],1, mean)
  
}

bulk_data$gene_Z_scores[["PN"]] = cluster_gene_mat




#calculate mean Z-score per gene by per cluster for each group comparison

for (comp_group in names(comp_groups)){
  
  message("Extracting cluster expression Z-score ", comp_group)
  
  cluster_gene_mat1 = cluster_gene_mat
  cluster_gene_mat2 = cluster_gene_mat
  
  for (cl in cluster_names){
    
    cl_bulks1 = meta$cluster_sample[meta$cluster_name == cl &
                                      meta$group == comp_groups[[comp_group]][1] ]
    cl_bulks2 = meta$cluster_sample[meta$cluster_name == cl &
                                      meta$group == comp_groups[[comp_group]][2] ]
    
    cluster_gene_mat1[,cl] = apply(cluster_sample_mat[,cl_bulks1],1, mean)
    cluster_gene_mat2[,cl] = apply(cluster_sample_mat[,cl_bulks2],1, mean)
    
  }
  
  bulk_data$gene_Z_scores[[comp_group]] = cluster_gene_mat1-cluster_gene_mat2
  
}

####

for (comp_group in names(comp_groups)){
  
  message("Extracting cluster expression Z-score ", comp_group)
  
  cluster_gene_mat1 = cluster_gene_mat
  cluster_gene_mat2 = cluster_gene_mat
  
  for (cl in cluster_names){
    
    cl_bulks1 = meta$cluster_sample[meta$cluster_name == cl &
                                      meta$group == comp_groups[[comp_group]][1] ]
    cl_bulks2 = meta$cluster_sample[meta$cluster_name == cl &
                                      meta$group == comp_groups[[comp_group]][2] ]
    
    # Add "X" prefix to match column names in the matrix
    cl_bulks1 = paste0("X", cl_bulks1)
    cl_bulks2 = paste0("X", cl_bulks2)
    
    # Filter for existing columns (safety check)
    cl_bulks1 = cl_bulks1[cl_bulks1 %in% colnames(cluster_sample_mat)]
    cl_bulks2 = cl_bulks2[cl_bulks2 %in% colnames(cluster_sample_mat)]
    
    if(length(cl_bulks1) > 0) {
      if(length(cl_bulks1) == 1) {
        cluster_gene_mat1[,cl] = cluster_sample_mat[,cl_bulks1]
      } else {
        cluster_gene_mat1[,cl] = apply(cluster_sample_mat[,cl_bulks1, drop=FALSE], 1, mean)
      }
    } else {
      warning(paste("No samples found for cluster", cl, "group", comp_groups[[comp_group]][1]))
    }
    
    if(length(cl_bulks2) > 0) {
      if(length(cl_bulks2) == 1) {
        cluster_gene_mat2[,cl] = cluster_sample_mat[,cl_bulks2]
      } else {
        cluster_gene_mat2[,cl] = apply(cluster_sample_mat[,cl_bulks2, drop=FALSE], 1, mean)
      }
    } else {
      warning(paste("No samples found for cluster", cl, "group", comp_groups[[comp_group]][2]))
    }
  }
  
  bulk_data$gene_Z_scores[[comp_group]] = cluster_gene_mat1 - cluster_gene_mat2
}

####


save(bulk_data, file = paste0(out_dir, script_ind, "bulk_data_w_expr_z_scores_deg_characterization.rda"))



#########################################
#plot DEGs by comp_group
#########################################

meta = bulk_data$meta


message("\n\n          *** Plotting DEGs cluster expression Z-score... ", Sys.time(),"\n\n")

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_DEGs_by_group_comp_cluster_z_score_CON_vs_comp.pdf"), 
    width = 25, height = 60)
{
  for (comp_group in names(comp_groups)){
    
    l1 = bulk_data$DEGs[grepl(comp_group, names(bulk_data$DEGs))]
    pl_genes = unique(unlist(l1))
    
    if(length(pl_genes)>1){
      
      pl_lims_max = 0.5*max(na.omit(bulk_data$gene_Z_scores[[comp_group]]))
      
      p1 = bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores[[comp_group]], 
                            pl_meta = bulk_data$meta,
                            pl_genes = pl_genes,
                            x_col = "cluster_name", 
                            meta_annot_cols = c("cell_type", "cluster_name"),
                            cluster_rows = TRUE, cluster_cols = FALSE,
                            color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                            pl_lims = c(-pl_lims_max, pl_lims_max),  
                            cellwidth = 12, cellheight = 10, fontsize = 10, 
                            title = paste0(comp_group, " - DEGs expression in ",comp_group," (Z-score)"))
      
      pl_genes_ordered = pl_genes[p1$tree_row$order]
      
      bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$PN, 
                       pl_meta = bulk_data$meta,
                       pl_genes = pl_genes_ordered,
                       x_col = "cluster_name", 
                       meta_annot_cols = c("cell_type", "cluster_name"),
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       color = viridis(250),
                       pl_lims = NULL,  cellwidth = 12, cellheight = 10, fontsize = 10, 
                       title = paste0(comp_group, " - DEGs expression in PN (Z-score)"))
      
    }
  }
}
dev.off()

#########################################
#plot genes by GO_term (max top20 GO terms)
#########################################

message("\n\n          *** Plotting DEGs cluster expression Z-score... ", Sys.time(),"\n\n")

pdf(file = paste0(out_dir,script_ind, "Gene_expr_heatmap_GO_DEGs_by_group_comp_cluster_z_score_CON_vs_comp.pdf"), 
    width = 25, height = 60)
{
  for (comp_group in names(comp_groups)){
    
    go_res = bulk_data$GO_results$full[[comp_group]]
    if (nrow(go_res) > 20){go_res = go_res[1:20,]}
    
    go_genes_list = str_split(go_res$geneID, "/")
    names(go_genes_list) = go_res$ID
    
    for (go in names(go_genes_list)){
      
      pl_genes = go_genes_list[[go]]
      
      if(length(pl_genes)>1){
        
        pl_lims_max = 0.5*max(na.omit(bulk_data$gene_Z_scores[[comp_group]]))
        
        p1 = bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores[[comp_group]], 
                              pl_meta = bulk_data$meta,
                              pl_genes = pl_genes,
                              x_col = "cluster_name", 
                              meta_annot_cols = c("cluster_name"),
                              cluster_rows = TRUE, cluster_cols = FALSE,
                              color = colorRampPalette(c("magenta", "black", "yellow"))(250),
                              pl_lims = c(-pl_lims_max, pl_lims_max),  
                              cellwidth = 12, cellheight = 10, fontsize = 10, 
                              title = paste0(go_res$Description[go_res$ID == go]," (", go,
                                             ") \n- DEG expression in ",comp_group," (Z-score)"))
        
        pl_genes_ordered = pl_genes[p1$tree_row$order]
        
        bulkdata_heatmap(pl_mat = bulk_data$gene_Z_scores$PN, 
                         pl_meta = bulk_data$meta,
                         pl_genes = pl_genes_ordered,
                         x_col = "cluster_name", 
                         meta_annot_cols = c("cell_type", "cluster_name"),
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         color = viridis(250),
                         pl_lims = NULL,  cellwidth = 12, cellheight = 10, fontsize = 10, 
                         title = paste0(go_res$Description[go_res$ID == go]," (", go,
                                        ") \n- DEG expression in PN (Z-score)"))
        
      }
    }
  }
}
dev.off()

#get info on version of R, used packages etc
sessionInfo()

message("\n\n##########################################################################\n",
        "# Completed 07 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")


#########################################
######Analysing Treatment Response#######
#########################################

#################################################
# Determine which clusters to analyze
#################################################
# First, extract all available clusters from the deseq_results list
get_available_clusters <- function(bulk_data) {
  # Get all keys from deseq_results
  all_keys <- names(bulk_data$deseq_results)
  
  # Extract cluster names from keys
  cluster_names <- unique(sapply(all_keys, function(key) {
    # Split by underscore and take the first element
    parts <- strsplit(key, "_")[[1]]
    return(parts[1])
  }))
  
  # Remove any non-cluster entries (like 'all')
  return(cluster_names[!cluster_names %in% c("all")])
}

# Get available clusters
comp_clusters <- get_available_clusters(bulk_data)
message("Found ", length(comp_clusters), " clusters to analyze: ", paste(comp_clusters, collapse=", "))


# Fixed function to analyze treatment impact on VCP-affected genes
analyze_treatment_impact <- function(bulk_data, clusters = comp_clusters) {
  # Store results 
  treatment_impact <- list()
  
  for (cl in clusters) {
    message("  * Analyzing treatment impact for cluster: ", cl)
    
    # 1. Get DEGs from VCP vs control comparison (note: direction is flipped in your data)
    ctrl_vs_vcp_key <- paste0(cl, "_VCP.u_vs_ctrl.u")  # This matches your actual keys
    ctrl_vs_vcp_results <- bulk_data$deseq_results[[ctrl_vs_vcp_key]]
    
    if (is.null(ctrl_vs_vcp_results)) {
      message("    - No results found for ", ctrl_vs_vcp_key)
      next
    }
    
    # Get genes up and down in VCP vs control (note the flipped comparison direction)
    vcp_up_genes <- bulk_data$DEGs[[paste0(ctrl_vs_vcp_key, "_down")]]  # These are UP in VCP
    vcp_down_genes <- bulk_data$DEGs[[paste0(ctrl_vs_vcp_key, "_up")]]  # These are DOWN in VCP
    
    message("    - Found ", length(vcp_up_genes), " genes up-regulated in VCP")
    message("    - Found ", length(vcp_down_genes), " genes down-regulated in VCP")
    
    # 2. Get results from VCP+SB vs control comparison 
    ctrl_vs_vcp_sb_key <- paste0(cl, "_VCP.SB_vs_ctrl.u")  # This matches your actual keys
    ctrl_vs_vcp_sb_results <- bulk_data$deseq_results[[ctrl_vs_vcp_sb_key]]
    
    if (is.null(ctrl_vs_vcp_sb_results)) {
      message("    - No results found for ", ctrl_vs_vcp_sb_key)
      next
    }
    
    # 3. Analyze VCP upregulated genes (how do they look in VCP+SB vs control)
    if (length(vcp_up_genes) > 0) {
      # Extract treatment effect on VCP upregulated genes
      vcp_up_treatment <- ctrl_vs_vcp_sb_results[vcp_up_genes, ]
      
      # Classify response by comparing log2FC magnitude and direction
      normalized_genes <- c()
      exacerbated_genes <- c()
      unaffected_genes <- c()
      
      for (gene in vcp_up_genes) {
        if (gene %in% rownames(ctrl_vs_vcp_sb_results)) {
          # Get log2FC values for both comparisons
          l2fc_vcp <- ctrl_vs_vcp_results[gene, "log2FoldChange"]
          l2fc_sb <- ctrl_vs_vcp_sb_results[gene, "log2FoldChange"]
          padj_sb <- ctrl_vs_vcp_sb_results[gene, "padj"]
          
          # Gene is significant in VCP+SB vs control
          if (!is.na(padj_sb) && padj_sb <= 0.05) {
            if (abs(l2fc_sb) < abs(l2fc_vcp) * 0.5) {
              # Effect size reduced by at least 50%
              normalized_genes <- c(normalized_genes, gene)
            } else if (abs(l2fc_sb) > abs(l2fc_vcp) * 1.5) {
              # Effect size increased by at least 50%
              exacerbated_genes <- c(exacerbated_genes, gene)
            } else {
              # Effect similar
              unaffected_genes <- c(unaffected_genes, gene)
            }
          } else {
            # No longer significant with treatment
            normalized_genes <- c(normalized_genes, gene)
          }
        }
      }
      
      message("    - VCP up-regulated genes normalized by treatment: ", length(normalized_genes))
      message("    - VCP up-regulated genes exacerbated by treatment: ", length(exacerbated_genes))
      message("    - VCP up-regulated genes unaffected by treatment: ", length(unaffected_genes))
      
      treatment_impact[[paste0(cl, "_vcp_up_normalized_by_treatment")]] <- normalized_genes
      treatment_impact[[paste0(cl, "_vcp_up_exacerbated_by_treatment")]] <- exacerbated_genes
      treatment_impact[[paste0(cl, "_vcp_up_unaffected_by_treatment")]] <- unaffected_genes
    }
    
    # 4. Analyze VCP downregulated genes (how do they look in VCP+SB vs control)
    if (length(vcp_down_genes) > 0) {
      # Extract treatment effect on VCP downregulated genes
      vcp_down_treatment <- ctrl_vs_vcp_sb_results[vcp_down_genes, ]
      
      # Classify response
      normalized_genes <- c()
      exacerbated_genes <- c()
      unaffected_genes <- c()
      
      for (gene in vcp_down_genes) {
        if (gene %in% rownames(ctrl_vs_vcp_sb_results)) {
          # Get log2FC values for both comparisons
          l2fc_vcp <- ctrl_vs_vcp_results[gene, "log2FoldChange"]
          l2fc_sb <- ctrl_vs_vcp_sb_results[gene, "log2FoldChange"]
          padj_sb <- ctrl_vs_vcp_sb_results[gene, "padj"]
          
          # Gene is significant in VCP+SB vs control
          if (!is.na(padj_sb) && padj_sb <= 0.05) {
            if (abs(l2fc_sb) < abs(l2fc_vcp) * 0.5) {
              # Effect size reduced by at least 50%
              normalized_genes <- c(normalized_genes, gene)
            } else if (abs(l2fc_sb) > abs(l2fc_vcp) * 1.5) {
              # Effect size increased by at least 50%
              exacerbated_genes <- c(exacerbated_genes, gene)
            } else {
              # Effect similar
              unaffected_genes <- c(unaffected_genes, gene)
            }
          } else {
            # No longer significant with treatment
            normalized_genes <- c(normalized_genes, gene)
          }
        }
      }
      
      message("    - VCP down-regulated genes normalized by treatment: ", length(normalized_genes))
      message("    - VCP down-regulated genes exacerbated by treatment: ", length(exacerbated_genes))
      message("    - VCP down-regulated genes unaffected by treatment: ", length(unaffected_genes))
      
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

# Generate summary statistics
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
  # Use the correct key patterns that match your data structure
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

# Create comparison plots to visualize gene expression changes
# Create scatter plots comparing log2FC in ctrl.u vs VCP.u to ctrl.u vs VCP.SB
library(ggplot2)
library(gridExtra)

create_comparison_scatter_plot <- function(bulk_data, cluster) {
  # Get DESeq2 results using correct key patterns
  ctrl_vs_vcp_key <- paste0(cluster, "_VCP.u_vs_ctrl.u")
  ctrl_vs_vcp_sb_key <- paste0(cluster, "_VCP.SB_vs_ctrl.u")
  
  # Extract results
  ctrl_vs_vcp_results <- bulk_data$deseq_results[[ctrl_vs_vcp_key]]
  ctrl_vs_vcp_sb_results <- bulk_data$deseq_results[[ctrl_vs_vcp_sb_key]]
  
  # Check if results exist
  if (is.null(ctrl_vs_vcp_results) || is.null(ctrl_vs_vcp_sb_results)) {
    message("Missing DESeq2 results for ", cluster)
    return(NULL)
  }
  
  # CORRECTED: Fixed the mapping between DEG lists and up/down regulation
  # The comparison is VCP.u vs ctrl.u, so:
  # - Positive log2FC means higher in VCP than control (up in VCP)
  # - Negative log2FC means lower in VCP than control (down in VCP)
  # "_up" in DEGs means genes with positive log2FC (up in the first condition vs second)
  # "_down" in DEGs means genes with negative log2FC (down in the first condition vs second)
  vcp_up_genes <- bulk_data$DEGs[[paste0(ctrl_vs_vcp_key, "_up")]]    # CORRECTED: up in VCP vs ctrl
  vcp_down_genes <- bulk_data$DEGs[[paste0(ctrl_vs_vcp_key, "_down")]] # CORRECTED: down in VCP vs ctrl
  
  # Create data frame for plotting
  plot_data <- data.frame(
    gene = rownames(ctrl_vs_vcp_results),
    l2fc_vcp = ctrl_vs_vcp_results$log2FoldChange,
    l2fc_sb = NA,
    gene_type = "Not Significant"
  )
  
  # Add VCP_SB log2FC values
  for (gene in plot_data$gene) {
    if (gene %in% rownames(ctrl_vs_vcp_sb_results)) {
      plot_data$l2fc_sb[plot_data$gene == gene] <- ctrl_vs_vcp_sb_results[gene, "log2FoldChange"]
    }
  }
  
  # Remove rows with NA values
  plot_data <- plot_data[!is.na(plot_data$l2fc_sb), ]
  
  # Mark gene types
  plot_data$gene_type[plot_data$gene %in% vcp_up_genes] <- "Up in VCP"
  plot_data$gene_type[plot_data$gene %in% vcp_down_genes] <- "Down in VCP"
  plot_data$gene_type <- factor(plot_data$gene_type, 
                                levels = c("Up in VCP", "Down in VCP", "Not Significant"))
  
  # Add normalized/exacerbated status - also update these to match the corrected gene lists
  plot_data$treatment_effect <- "Unaffected"
  plot_data$treatment_effect[plot_data$gene %in% 
                               treatment_impact[[paste0(cluster, "_vcp_up_normalized_by_treatment")]]] <- "Normalized"
  plot_data$treatment_effect[plot_data$gene %in% 
                               treatment_impact[[paste0(cluster, "_vcp_up_exacerbated_by_treatment")]]] <- "Exacerbated"
  plot_data$treatment_effect[plot_data$gene %in% 
                               treatment_impact[[paste0(cluster, "_vcp_down_normalized_by_treatment")]]] <- "Normalized"
  plot_data$treatment_effect[plot_data$gene %in% 
                               treatment_impact[[paste0(cluster, "_vcp_down_exacerbated_by_treatment")]]] <- "Exacerbated"
  
  # Keep the original color scheme since the gene annotations are now corrected
  p <- ggplot(plot_data, aes(x = l2fc_vcp, y = l2fc_sb, color = gene_type)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "darkgrey") +
    labs(
      title = paste0("Gene expression changes in ", cluster),
      x = "log2FC (VCP.u vs ctrl.u)",
      y = "log2FC (VCP.SB vs ctrl.u)",
      color = "Gene type"
    ) +
    scale_color_manual(values = c("Up in VCP" = "firebrick", 
                                  "Down in VCP" = "steelblue", 
                                  "Not Significant" = "grey80")) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
  
  # Extract genes that have changed significantly after treatment
  # Specifically looking at genes that moved closer to control levels
  normalized_genes <- c(
    treatment_impact[[paste0(cluster, "_vcp_up_normalized_by_treatment")]],
    treatment_impact[[paste0(cluster, "_vcp_down_normalized_by_treatment")]]
  )
  
  # Highlight the top 10 most normalized genes
  if (length(normalized_genes) > 0) {
    normalized_data <- plot_data[plot_data$gene %in% normalized_genes, ]
    
    # Calculate normalization magnitude (distance to diagonal line)
    normalized_data$norm_magnitude <- abs(abs(normalized_data$l2fc_sb) - abs(normalized_data$l2fc_vcp))
    normalized_data <- normalized_data[order(normalized_data$norm_magnitude, decreasing = TRUE), ]
    
    # Take top 10 or all if less than 10
    top_n <- min(10, nrow(normalized_data))
    if (top_n > 0) {  # Check if any normalized data exists
      top_normalized <- normalized_data[1:top_n, ]
      
      # Add labels for top normalized genes
      p <- p + geom_text(
        data = top_normalized,
        aes(label = gene),
        size = 3,
        vjust = -1.5,
        check_overlap = TRUE
      )
    }
  }
  
  return(p)
}

# Create plots for each cluster
scatter_plots <- list()
for (cl in comp_clusters) {
  plot <- create_comparison_scatter_plot(bulk_data, cl)
  
  if (!is.null(plot)) {
    scatter_plots[[cl]] <- plot
    
    # Save individual plots
    ggsave(paste0(main_dir, script_ind, "treatment_comparison_", cl, ".pdf"), 
           scatter_plots[[cl]], width = 8, height = 7)
  }
}

# Combine plots into a multi-page PDF if there are multiple clusters
if (length(scatter_plots) > 1) {
  pdf(paste0(main_dir, script_ind, "treatment_comparison_all_clusters.pdf"), 
      width = 10, height = 8)
  for (cl in names(scatter_plots)) {
    print(scatter_plots[[cl]])
  }
  dev.off()
}

# Create a heatmap to visualize top normalized genes across clusters
# Identify top normalized genes across all clusters
get_top_normalized_genes <- function(treatment_impact, bulk_data, n_per_cluster = 5) {
  top_genes <- list()
  
  for (cl in comp_clusters) {
    # Get normalized genes from both up and down regulated categories
    up_norm_key <- paste0(cl, "_vcp_up_normalized_by_treatment")
    down_norm_key <- paste0(cl, "_vcp_down_normalized_by_treatment")
    
    normalized_genes <- c()
    
    if (up_norm_key %in% names(treatment_impact)) {
      normalized_genes <- c(normalized_genes, treatment_impact[[up_norm_key]])
    }
    
    if (down_norm_key %in% names(treatment_impact)) {
      normalized_genes <- c(normalized_genes, treatment_impact[[down_norm_key]])
    }
    
    if (length(normalized_genes) > 0) {
      # Get log2FC values from both comparisons
      ctrl_vs_vcp_key <- paste0(cl, "_VCP.u_vs_ctrl.u")  # Updated key pattern
      ctrl_vs_vcp_sb_key <- paste0(cl, "_VCP.SB_vs_ctrl.u")  # Updated key pattern
      
      vcp_results <- bulk_data$deseq_results[[ctrl_vs_vcp_key]]
      vcp_sb_results <- bulk_data$deseq_results[[ctrl_vs_vcp_sb_key]]
      
      if (!is.null(vcp_results) && !is.null(vcp_sb_results)) {
        # Calculate normalization magnitude
        norm_magnitude <- numeric(length(normalized_genes))
        names(norm_magnitude) <- normalized_genes
        
        for (i in seq_along(normalized_genes)) {
          gene <- normalized_genes[i]
          if (gene %in% rownames(vcp_results) && gene %in% rownames(vcp_sb_results)) {
            l2fc_vcp <- vcp_results[gene, "log2FoldChange"]
            l2fc_sb <- vcp_sb_results[gene, "log2FoldChange"]
            
            # Calculate how much closer to zero (control level) the gene moved
            norm_magnitude[i] <- abs(l2fc_vcp) - abs(l2fc_sb)
          }
        }
        
        # Sort by normalization magnitude
        norm_magnitude <- sort(norm_magnitude, decreasing = TRUE)
        
        # Take top n genes
        top_genes[[cl]] <- names(head(norm_magnitude, n_per_cluster))
      }
    }
  }
  
  # Combine genes from all clusters
  all_top_genes <- unique(unlist(top_genes))
  return(all_top_genes)
}

# Get top normalized genes
top_normalized_genes <- get_top_normalized_genes(treatment_impact, bulk_data)

# Create a heatmap if there are any top genes
if (length(top_normalized_genes) > 0) {
  library(pheatmap)
  
  # Extract data for heatmap
  # For each gene, get log2FC for both VCP.u_vs_ctrl.u and VCP.SB_vs_ctrl.u comparisons
  heatmap_data <- matrix(NA, nrow = length(top_normalized_genes), 
                         ncol = length(comp_clusters) * 2)
  rownames(heatmap_data) <- top_normalized_genes
  
  colnames_vec <- character(length(comp_clusters) * 2)
  col_idx <- 1
  
  for (cl in comp_clusters) {
    # Get VCP.u_vs_ctrl.u comparison
    ctrl_vs_vcp_key <- paste0(cl, "_VCP.u_vs_ctrl.u")
    deseq_res_1 <- bulk_data$deseq_results[[ctrl_vs_vcp_key]]
    
    # Get VCP.SB_vs_ctrl.u comparison
    ctrl_vs_vcp_sb_key <- paste0(cl, "_VCP.SB_vs_ctrl.u")
    deseq_res_2 <- bulk_data$deseq_results[[ctrl_vs_vcp_sb_key]]
    
    if (!is.null(deseq_res_1) && !is.null(deseq_res_2)) {
      # Extract log2FC for top genes
      for (i in seq_along(top_normalized_genes)) {
        gene <- top_normalized_genes[i]
        
        if (gene %in% rownames(deseq_res_1)) {
          heatmap_data[i, col_idx] <- deseq_res_1[gene, "log2FoldChange"]
        }
        
        if (gene %in% rownames(deseq_res_2)) {
          heatmap_data[i, col_idx + 1] <- deseq_res_2[gene, "log2FoldChange"]
        }
      }
      
      colnames_vec[col_idx] <- paste0(cl, "_VCP.u_vs_ctrl.u")
      colnames_vec[col_idx + 1] <- paste0(cl, "_VCP.SB_vs_ctrl.u")
      
      col_idx <- col_idx + 2
    }
  }
  
  # Remove NA columns if any
  valid_cols <- !is.na(colnames_vec)
  heatmap_data <- heatmap_data[, valid_cols]
  colnames_vec <- colnames_vec[valid_cols]
  
  colnames(heatmap_data) <- colnames_vec
  
  # Create annotation for columns to distinguish VCP effect vs treatment effect
  column_anno <- data.frame(
    comparison_type = factor(
      ifelse(grepl("_VCP.u_vs_ctrl.u$", colnames_vec), "VCP effect", "VCP+SB effect"),
      levels = c("VCP effect", "VCP+SB effect")
    )
  )
  rownames(column_anno) <- colnames_vec
  
  # Create heatmap
  pheatmap(
    heatmap_data,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = column_anno,
    annotation_colors = list(
      comparison_type = c("VCP effect" = "dodgerblue", "VCP+SB effect" = "darkorange")
    ),
    main = "Top normalized genes: VCP vs VCP+SB effect",
    filename = paste0(main_dir, script_ind, "top_normalized_genes_heatmap.pdf"),
    width = 12,
    height = 10
  )
}

# Save the updated bulk_data object
save(bulk_data, file = paste0(main_dir, script_ind, "bulk_data_with_treatment_impact_treatment_analysis.rda"))

# Print summary table
print(treatment_summary)

##############################################
#####Assessing Treatment Normalized Genes#####
##############################################

# Assuming treatment_summary is your data frame from the output
# First, let's reshape the data for better visualization
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(patchwork)

treatment_long <- treatment_summary %>%
  select(-vcp_up_total, -vcp_down_total) %>%
  pivot_longer(
    cols = c(vcp_up_normalized, vcp_up_exacerbated, vcp_up_unaffected,
             vcp_down_normalized, vcp_down_exacerbated, vcp_down_unaffected),
    names_to = "category",
    values_to = "count"
  ) %>%
  # Create separate columns for direction and effect
  mutate(
    direction = ifelse(grepl("up", category), "Upregulated in VCP", "Downregulated in VCP"),
    effect = case_when(
      grepl("normalized", category) ~ "Normalized",
      grepl("exacerbated", category) ~ "Exacerbated",
      grepl("unaffected", category) ~ "Unaffected"
    )
  )

# Calculate percentages within each direction for each cluster
treatment_pct <- treatment_long %>%
  group_by(cluster, direction) %>%
  mutate(total = sum(count),
         percentage = ifelse(total > 0, count / total * 100, 0)) %>%
  ungroup()

# Create a stacked bar chart showing treatment effect by cluster
p1 <- ggplot(treatment_pct, aes(x = cluster, y = count, fill = effect)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ direction, scales = "free_y") +
  scale_fill_manual(values = c("Normalized" = "#1B9E77", "Exacerbated" = "#D95F02", "Unaffected" = "#7570B3")) +
  labs(
    title = "Treatment Impact on VCP-affected Genes",
    x = "Cluster",
    y = "Gene Count",
    fill = "Treatment Effect"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

# Create a percentage chart
p2 <- ggplot(treatment_pct, aes(x = cluster, y = percentage, fill = effect)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ direction) +
  scale_fill_manual(values = c("Normalized" = "#1B9E77", "Exacerbated" = "#D95F02", "Unaffected" = "#7570B3")) +
  labs(
    title = "Percentage of VCP-affected Genes by Treatment Effect",
    x = "Cluster",
    y = "Percentage of Genes (%)",
    fill = "Treatment Effect"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

# Create a heatmap showing the percentage of normalized genes across clusters
heatmap_data <- treatment_pct %>%
  filter(effect == "Normalized") %>%
  select(cluster, direction, percentage)

p3 <- ggplot(heatmap_data, aes(x = cluster, y = direction, fill = percentage)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#1B9E77", name = "% Normalized") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), color = "black", size = 3) +
  labs(
    title = "Percentage of Genes Normalized by Treatment Across Clusters",
    x = "Cluster",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# Create a dot plot showing the ratio of normalized to exacerbated genes
ratio_data <- treatment_long %>%
  filter(effect %in% c("Normalized", "Exacerbated")) %>%
  group_by(cluster, direction, effect) %>%
  summarize(count = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = effect, values_from = count) %>%
  mutate(
    Normalized = ifelse(is.na(Normalized), 0, Normalized),
    Exacerbated = ifelse(is.na(Exacerbated), 0, Exacerbated),
    ratio = ifelse(Exacerbated == 0, 
                   ifelse(Normalized > 0, Inf, 0), 
                   Normalized / Exacerbated),
    label = ifelse(is.finite(ratio), 
                   ifelse(ratio > 100, ">100", sprintf("%.1f", ratio)), 
                   ifelse(Normalized > 0, "∞", "NA"))
  )

p4 <- ggplot(ratio_data, aes(x = cluster, y = ratio, color = direction)) +
  geom_point(size = 3) +
  geom_text(aes(label = label), vjust = -1, size = 3) +
  scale_y_log10(limits = c(0.1, NA), 
                breaks = c(0.1, 1, 10, 100, 1000),
                labels = c("0.1", "1", "10", "100", "1000"),
                oob = scales::squish_infinite) +
  labs(
    title = "Ratio of Normalized to Exacerbated Genes",
    x = "Cluster",
    y = "Normalized:Exacerbated Ratio (log scale)",
    color = "Gene Direction"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

# Print the plots
print(p1)
print(p2)
print(p3)
print(p4)

# To save the plots
ggsave(paste0(main_dir, script_ind, "treatment_impact_counts.pdf"), p1, width = 10, height = 6)
ggsave(paste0(main_dir, script_ind, "treatment_impact_percentages.pdf"), p2, width = 10, height = 6)
ggsave(paste0(main_dir, script_ind, "normalization_heatmap.pdf"), p3, width = 10, height = 4)
ggsave(paste0(main_dir, script_ind, "normalization_ratio.pdf"), p4, width = 10, height = 6)

# Combine multiple plots in one figure with patchwork
combined_plot <- (p1 / p2) | (p3 / p4)
combined_plot + plot_annotation(
  title = "Comprehensive Analysis of Treatment Impact on VCP-affected Genes",
  theme = theme(plot.title = element_text(size = 14, hjust = 0.5))
)

# Save the combined plot
ggsave(paste0(main_dir, script_ind, "treatment_impact_summary.pdf"), combined_plot, width = 16, height = 12)


###New Visualisation####

# This script provides a complete analysis of how SB treatment affects the VCP mutant phenotype
# Just run this after you've loaded your bulk_data object with DESeq2 results

#################################################
# Determine which clusters to analyze
#################################################
# First, extract all available clusters from the deseq_results list
get_available_clusters <- function(bulk_data) {
  # Get all keys from deseq_results
  all_keys <- names(bulk_data$deseq_results)
  
  # Extract cluster names from keys
  cluster_names <- unique(sapply(all_keys, function(key) {
    # Split by underscore and take the first element
    parts <- strsplit(key, "_")[[1]]
    return(parts[1])
  }))
  
  # Remove any non-cluster entries (like 'all')
  return(cluster_names[!cluster_names %in% c("all")])
}

# Get available clusters
comp_clusters <- get_available_clusters(bulk_data)
message("Found ", length(comp_clusters), " clusters to analyze: ", paste(comp_clusters, collapse=", "))

#################################################
# Analyze treatment impact on VCP-affected genes
#################################################

#################################################
# Analyze treatment impact on VCP-affected genes
#################################################

analyze_treatment_impact <- function(bulk_data, clusters) {
  # Store results 
  treatment_impact <- list()
  
  for (cl in clusters) {
    message("  * Analyzing treatment impact for cluster: ", cl)
    
    # Check if required comparisons exist - using the correct naming convention
    ctrl_vs_vcp_key <- paste0(cl, "_VCP.u_vs_ctrl.u")
    ctrl_vs_vcp_sb_key <- paste0(cl, "_VCP.SB_vs_ctrl.u")
    
    if (!(ctrl_vs_vcp_key %in% names(bulk_data$deseq_results)) || 
        !(ctrl_vs_vcp_sb_key %in% names(bulk_data$deseq_results))) {
      message("    - Skipping cluster ", cl, " - missing required comparisons")
      next
    }
    
    # 1. Get DEGs from VCP vs control comparison
    ctrl_vs_vcp_results <- bulk_data$deseq_results[[ctrl_vs_vcp_key]]
    
    # Check if DEG lists exist
    up_key <- paste0(ctrl_vs_vcp_key, "_up")
    down_key <- paste0(ctrl_vs_vcp_key, "_down")
    
    if (!(up_key %in% names(bulk_data$DEGs)) || 
        !(down_key %in% names(bulk_data$DEGs))) {
      message("    - Skipping cluster ", cl, " - missing DEG lists")
      next
    }
    
    # Get genes up and down in VCP vs control (note direction is flipped in your data)
    vcp_up_genes <- bulk_data$DEGs[[down_key]]  # down in VCP.u_vs_ctrl.u means up in VCP
    vcp_down_genes <- bulk_data$DEGs[[up_key]]  # up in VCP.u_vs_ctrl.u means down in VCP
    
    message("    - Found ", length(vcp_up_genes), " genes up-regulated in VCP")
    message("    - Found ", length(vcp_down_genes), " genes down-regulated in VCP")
    
    # 2. Get results from VCP+SB vs control comparison
    ctrl_vs_vcp_sb_results <- bulk_data$deseq_results[[ctrl_vs_vcp_sb_key]]
    
    # 3. Analyze VCP upregulated genes (how do they look in VCP+SB vs control)
    if (length(vcp_up_genes) > 0) {
      # Classify response by comparing log2FC magnitude and direction
      normalized_genes <- c()
      exacerbated_genes <- c()
      unaffected_genes <- c()
      
      for (gene in vcp_up_genes) {
        if (gene %in% rownames(ctrl_vs_vcp_sb_results)) {
          # Get log2FC values for both comparisons
          l2fc_vcp <- ctrl_vs_vcp_results[gene, "log2FoldChange"]
          l2fc_sb <- ctrl_vs_vcp_sb_results[gene, "log2FoldChange"]
          padj_sb <- ctrl_vs_vcp_sb_results[gene, "padj"]
          
          # Gene is significant in VCP+SB vs control
          if (!is.na(padj_sb) && padj_sb <= 0.05) {
            if (abs(l2fc_sb) < abs(l2fc_vcp) * 0.5) {
              # Effect size reduced by at least 50%
              normalized_genes <- c(normalized_genes, gene)
            } else if (abs(l2fc_sb) > abs(l2fc_vcp) * 1.5) {
              # Effect size increased by at least 50%
              exacerbated_genes <- c(exacerbated_genes, gene)
            } else {
              # Effect similar
              unaffected_genes <- c(unaffected_genes, gene)
            }
          } else {
            # No longer significant with treatment
            normalized_genes <- c(normalized_genes, gene)
          }
        }
      }
      
      message("    - VCP up-regulated genes normalized by treatment: ", length(normalized_genes))
      message("    - VCP up-regulated genes exacerbated by treatment: ", length(exacerbated_genes))
      message("    - VCP up-regulated genes unaffected by treatment: ", length(unaffected_genes))
      
      treatment_impact[[paste0(cl, "_vcp_up_normalized_by_treatment")]] <- normalized_genes
      treatment_impact[[paste0(cl, "_vcp_up_exacerbated_by_treatment")]] <- exacerbated_genes
      treatment_impact[[paste0(cl, "_vcp_up_unaffected_by_treatment")]] <- unaffected_genes
    }
    
    # 4. Analyze VCP downregulated genes (how do they look in VCP+SB vs control)
    if (length(vcp_down_genes) > 0) {
      # Classify response
      normalized_genes <- c()
      exacerbated_genes <- c()
      unaffected_genes <- c()
      
      for (gene in vcp_down_genes) {
        if (gene %in% rownames(ctrl_vs_vcp_sb_results)) {
          # Get log2FC values for both comparisons
          l2fc_vcp <- ctrl_vs_vcp_results[gene, "log2FoldChange"]
          l2fc_sb <- ctrl_vs_vcp_sb_results[gene, "log2FoldChange"]
          padj_sb <- ctrl_vs_vcp_sb_results[gene, "padj"]
          
          # Gene is significant in VCP+SB vs control
          if (!is.na(padj_sb) && padj_sb <= 0.05) {
            if (abs(l2fc_sb) < abs(l2fc_vcp) * 0.5) {
              # Effect size reduced by at least 50%
              normalized_genes <- c(normalized_genes, gene)
            } else if (abs(l2fc_sb) > abs(l2fc_vcp) * 1.5) {
              # Effect size increased by at least 50%
              exacerbated_genes <- c(exacerbated_genes, gene)
            } else {
              # Effect similar
              unaffected_genes <- c(unaffected_genes, gene)
            }
          } else {
            # No longer significant with treatment
            normalized_genes <- c(normalized_genes, gene)
          }
        }
      }
      
      message("    - VCP down-regulated genes normalized by treatment: ", length(normalized_genes))
      message("    - VCP down-regulated genes exacerbated by treatment: ", length(exacerbated_genes))
      message("    - VCP down-regulated genes unaffected by treatment: ", length(unaffected_genes))
      
      treatment_impact[[paste0(cl, "_vcp_down_normalized_by_treatment")]] <- normalized_genes
      treatment_impact[[paste0(cl, "_vcp_down_exacerbated_by_treatment")]] <- exacerbated_genes
      treatment_impact[[paste0(cl, "_vcp_down_unaffected_by_treatment")]] <- unaffected_genes
    }
  }
  
  # Print summary of findings across all clusters
  message("\n=== Treatment Impact Summary ===")
  message("Clusters analyzed: ", length(unique(gsub("_.*", "", names(treatment_impact)))))
  
  up_norm_count <- sum(sapply(names(treatment_impact)[grepl("_vcp_up_normalized_by_treatment", names(treatment_impact))], 
                              function(x) length(treatment_impact[[x]])))
  down_norm_count <- sum(sapply(names(treatment_impact)[grepl("_vcp_down_normalized_by_treatment", names(treatment_impact))], 
                                function(x) length(treatment_impact[[x]])))
  
  message("Total up-regulated genes normalized: ", up_norm_count)
  message("Total down-regulated genes normalized: ", down_norm_count)
  message("Total genes normalized: ", up_norm_count + down_norm_count)
  
  return(treatment_impact)
}

# Run the treatment impact analysis
message("\n\n          *** Analyzing treatment impact on VCP-affected genes... ", Sys.time(),"\n\n")
treatment_impact <- analyze_treatment_impact(bulk_data, comp_clusters)
bulk_data$treatment_impact <- treatment_impact

# Generate summary statistics
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
  # Check if this cluster was analyzed
  up_norm_key <- paste0(cl, "_vcp_up_normalized_by_treatment")
  if (!(up_norm_key %in% names(treatment_impact))) {
    next
  }
  
  # Get DEG lists using correct naming convention
  ctrl_vs_vcp_key <- paste0(cl, "_VCP.u_vs_ctrl.u")
  up_key <- paste0(ctrl_vs_vcp_key, "_up")
  down_key <- paste0(ctrl_vs_vcp_key, "_down")
  
  # Count the genes in each category
  vcp_up_total <- length(bulk_data$DEGs[[down_key]])  # down in VCP.u_vs_ctrl.u means up in VCP
  vcp_down_total <- length(bulk_data$DEGs[[up_key]])  # up in VCP.u_vs_ctrl.u means down in VCP
  
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
print(treatment_summary)

#################################################
# Create visualization functions
#################################################
create_normalization_plot <- function(bulk_data, treatment_impact, cluster) {
  # Get DESeq2 results - using correct key patterns
  ctrl_vs_vcp_key <- paste0(cluster, "_VCP.u_vs_ctrl.u")
  ctrl_vs_vcp_sb_key <- paste0(cluster, "_VCP.SB_vs_ctrl.u")
  
  # Check if both comparisons exist in the results
  if (!(ctrl_vs_vcp_key %in% names(bulk_data$deseq_results)) || 
      !(ctrl_vs_vcp_sb_key %in% names(bulk_data$deseq_results))) {
    message("  * Skipping plot for cluster ", cluster, " - missing required comparisons")
    return(NULL)
  }
  
  # Extract results
  ctrl_vs_vcp_results <- bulk_data$deseq_results[[ctrl_vs_vcp_key]]
  ctrl_vs_vcp_sb_results <- bulk_data$deseq_results[[ctrl_vs_vcp_sb_key]]
  
  # Get DEGs in VCP vs control
  up_key <- paste0(ctrl_vs_vcp_key, "_up")
  down_key <- paste0(ctrl_vs_vcp_key, "_down")
  
  # Check if DEGs exist in the results
  if (!(up_key %in% names(bulk_data$DEGs)) || 
      !(down_key %in% names(bulk_data$DEGs))) {
    message("  * Skipping plot for cluster ", cluster, " - missing DEG lists")
    return(NULL)
  }
  
  vcp_degs <- c(
    bulk_data$DEGs[[up_key]], # Up in VCP.u_vs_ctrl.u (down in VCP vs control)
    bulk_data$DEGs[[down_key]] # Down in VCP.u_vs_ctrl.u (up in VCP vs control)
  )
  
  if (length(vcp_degs) == 0) {
    message("  * Skipping plot for cluster ", cluster, " - no DEGs found")
    return(NULL)
  }
  
  # Create a data frame for DEGs only
  plot_data <- data.frame(
    gene = vcp_degs,
    l2fc_vcp = NA,  # Will store log2FC for VCP.u_vs_ctrl.u
    l2fc_sb = NA,   # Will store log2FC for VCP.SB_vs_ctrl.u
    gene_type = NA,  # Will store whether gene is up or down in VCP
    stringsAsFactors = FALSE
  )
  
  # Fill in log2FC values and gene types
  for (i in 1:nrow(plot_data)) {
    gene <- plot_data$gene[i]
    
    # Get log2FC values
    if (gene %in% rownames(ctrl_vs_vcp_results)) {
      plot_data$l2fc_vcp[i] <- ctrl_vs_vcp_results[gene, "log2FoldChange"]
      
      # Determine gene type based on log2FC direction
      if (plot_data$l2fc_vcp[i] > 0) {
        plot_data$gene_type[i] <- "Down in VCP" # Positive log2FC in VCP.u_vs_ctrl.u = lower in VCP
      } else {
        plot_data$gene_type[i] <- "Up in VCP"   # Negative log2FC in VCP.u_vs_ctrl.u = higher in VCP
      }
    }
    
    if (gene %in% rownames(ctrl_vs_vcp_sb_results)) {
      plot_data$l2fc_sb[i] <- ctrl_vs_vcp_sb_results[gene, "log2FoldChange"]
    }
  }
  
  # Remove rows with missing values
  plot_data <- plot_data[!is.na(plot_data$l2fc_vcp) & !is.na(plot_data$l2fc_sb), ]
  
  if (nrow(plot_data) == 0) {
    message("  * Skipping plot for cluster ", cluster, " - no genes with data in both comparisons")
    return(NULL)
  }
  
  # Add a normalization score
  # This measures how much the absolute fold change is reduced by treatment
  # 1 = perfect normalization (gene returned to control level)
  # 0 = no effect (gene still at VCP level)
  # <0 = exacerbation (gene moved further from control)
  plot_data$norm_score <- 1 - (abs(plot_data$l2fc_sb) / abs(plot_data$l2fc_vcp))
  
  # Create a column for treatment effect category
  plot_data$effect <- "No Change"
  plot_data$effect[plot_data$norm_score > 0.25] <- "Normalized"
  plot_data$effect[plot_data$norm_score < -0.25] <- "Exacerbated"
  plot_data$effect <- factor(plot_data$effect, levels = c("Normalized", "No Change", "Exacerbated"))
  
  # Count genes in each category
  n_normalized <- sum(plot_data$effect == "Normalized")
  n_no_change <- sum(plot_data$effect == "No Change")
  n_exacerbated <- sum(plot_data$effect == "Exacerbated")
  
  # Order genes by normalization score for better visualization
  plot_data <- plot_data[order(plot_data$norm_score), ]
  
  # Add row number for plotting
  plot_data$index <- 1:nrow(plot_data)
  
  # Create first plot: Normalization scores
  library(ggplot2)
  p1 <- ggplot(plot_data, aes(x = index, y = norm_score, color = effect)) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Normalized" = "forestgreen", 
                                  "No Change" = "gray60", 
                                  "Exacerbated" = "firebrick")) +
    labs(
      title = paste0("Treatment Effect in Cluster ", cluster),
      subtitle = paste0(n_normalized, " genes normalized (", round(n_normalized/nrow(plot_data)*100), "%), ", 
                        n_exacerbated, " genes exacerbated (", round(n_exacerbated/nrow(plot_data)*100), "%)"),
      x = "Genes (ordered by normalization score)",
      y = "Normalization Score",
      color = "Effect"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    annotate("text", x = nrow(plot_data) * 0.05, y = 0.9, 
             label = "Perfect normalization = 1", hjust = 0, size = 3.5) +
    annotate("text", x = nrow(plot_data) * 0.05, y = 0.1, 
             label = "No effect = 0", hjust = 0, size = 3.5) +
    annotate("text", x = nrow(plot_data) * 0.05, y = -0.35, 
             label = "Exacerbation < 0", hjust = 0, size = 3.5)
  
  # Create second plot: Before-after plot for top normalized genes
  # Select top normalized genes (highest positive norm_score)
  top_n <- min(20, sum(plot_data$effect == "Normalized"))
  
  if (top_n > 0) {
    top_genes <- plot_data[plot_data$effect == "Normalized", ]
    top_genes <- top_genes[order(top_genes$norm_score, decreasing = TRUE), ][1:top_n, ]
    
    # Reshape data for paired plot
    library(reshape2)
    top_genes_long <- melt(top_genes[, c("gene", "l2fc_vcp", "l2fc_sb")],
                           id.vars = "gene",
                           variable.name = "condition",
                           value.name = "log2fc")
    
    # Rename conditions for clarity
    top_genes_long$condition <- as.character(top_genes_long$condition)
    top_genes_long$condition[top_genes_long$condition == "l2fc_vcp"] <- "VCP Mutant"
    top_genes_long$condition[top_genes_long$condition == "l2fc_sb"] <- "VCP+SB Treatment"
    top_genes_long$condition <- factor(top_genes_long$condition, 
                                       levels = c("VCP Mutant", "VCP+SB Treatment"))
    
    # Create "distance from control" plot
    p2 <- ggplot(top_genes_long, aes(x = reorder(gene, -abs(log2fc)), y = abs(log2fc), 
                                     fill = condition, group = gene)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      scale_fill_manual(values = c("VCP Mutant" = "firebrick", "VCP+SB Treatment" = "forestgreen")) +
      labs(
        title = paste0("Top Normalized Genes in Cluster ", cluster),
        subtitle = "Showing reduction in differential expression with treatment",
        x = "Gene",
        y = "Distance from Control (|log2FC|)",
        fill = "Condition"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  } else {
    # If no normalized genes, create empty plot with message
    p2 <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = "No normalized genes found") +
      theme_void()
  }
  
  # Create third plot: Direction of change with SB treatment
  # Prepare data frame for direction plot
  arrows_data <- plot_data
  arrows_data$direction <- "Towards Control"
  arrows_data$direction[sign(arrows_data$l2fc_vcp) == sign(arrows_data$l2fc_sb) & 
                          abs(arrows_data$l2fc_sb) > abs(arrows_data$l2fc_vcp)] <- "Away from Control"
  arrows_data$direction[sign(arrows_data$l2fc_vcp) != sign(arrows_data$l2fc_sb)] <- "Past Control"
  
  # Count genes in each direction
  n_towards <- sum(arrows_data$direction == "Towards Control")
  n_away <- sum(arrows_data$direction == "Away from Control")
  n_past <- sum(arrows_data$direction == "Past Control")
  
  # Create pie chart
  direction_counts <- data.frame(
    direction = c("Towards Control", "Away from Control", "Past Control"),
    count = c(n_towards, n_away, n_past)
  )
  direction_counts$percentage <- direction_counts$count / sum(direction_counts$count) * 100
  direction_counts$label <- paste0(direction_counts$direction, "\n", 
                                   round(direction_counts$percentage), "%")
  
  p3 <- ggplot(direction_counts, aes(x = "", y = count, fill = direction)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c("Towards Control" = "forestgreen", 
                                 "Away from Control" = "firebrick", 
                                 "Past Control" = "goldenrod")) +
    labs(
      title = paste0("Direction of Gene Expression Changes\nwith SB Treatment in Cluster ", cluster),
      fill = "Direction"
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    ) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5))
  
  # Create main scatter plot comparing log2FC values
  p4 <- ggplot(plot_data, aes(x = l2fc_vcp, y = l2fc_sb, color = effect)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "darkgrey") +
    scale_color_manual(values = c("Normalized" = "forestgreen", 
                                  "No Change" = "gray60", 
                                  "Exacerbated" = "firebrick")) +
    labs(
      title = paste0("Gene Expression Changes in Cluster ", cluster),
      subtitle = "Points below diagonal line show movement toward control expression",
      x = "log2FC (VCP.u vs ctrl.u)",
      y = "log2FC (VCP.SB vs ctrl.u)",
      color = "Treatment Effect"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  # Check if we should also add gene labels on the scatter plot
  # Add labels for the top 10 most normalized genes
  if (n_normalized >= 5) {
    # Get top normalized genes (highest norm_score)
    label_genes <- plot_data[plot_data$effect == "Normalized", ]
    label_genes <- label_genes[order(label_genes$norm_score, decreasing = TRUE), ]
    label_genes <- head(label_genes, 10)
    
    # Add labels to scatter plot
    p4 <- p4 + 
      geom_text_repel(
        data = label_genes,
        aes(label = gene),
        size = 3,
        box.padding = 0.5,
        point.padding = 0.2,
        force = 2,
        segment.color = "gray50",
        max.overlaps = 15
      )
  }
  
  # Return all plots
  return(list(
    normalization_plot = p1,
    top_genes_plot = p2,
    direction_plot = p3,
    scatter_plot = p4
  ))
}

#################################################
# Create plots for each cluster
#################################################
message("\n\n          *** Creating normalization plots... ", Sys.time(),"\n\n")

# Create output directory if it doesn't exist
if (!dir.exists(paste0(main_dir, script_ind))) {
  dir.create(paste0(main_dir, script_ind), recursive = TRUE)
}

# Process each cluster
library(gridExtra)
for (cl in comp_clusters) {
  message("  * Creating plots for cluster: ", cl)
  
  # Skip if this cluster wasn't analyzed in treatment_impact
  up_norm_key <- paste0(cl, "_vcp_up_normalized_by_treatment")
  if (!(up_norm_key %in% names(treatment_impact))) {
    message("    - Skipping plot for cluster ", cl, " - not analyzed in treatment_impact")
    next
  }
  
  plots <- create_normalization_plot(bulk_data, treatment_impact, cl)
  
  # Skip if no plots were generated
  if (is.null(plots)) {
    message("    - No plots generated for cluster ", cl)
    next
  }
  
  # Combine plots into a single figure
  combined_plot <- grid.arrange(
    plots$scatter_plot,
    plots$direction_plot,
    plots$normalization_plot, 
    plots$top_genes_plot,
    layout_matrix = rbind(c(1,2), c(3,4)),
    heights = c(1, 1)
  )
  
  # Save the combined plot
  output_file <- paste0(main_dir, script_ind, "treatment_effect_", cl, ".pdf")
  ggsave(output_file, combined_plot, width = 12, height = 12)
  message("    - Saved plot to: ", output_file)
}

#################################################
# Create a summary plot across all clusters
#################################################
# Create a stacked bar plot for all clusters
summary_data <- data.frame(
  cluster = character(),
  normalized_pct = numeric(),
  no_change_pct = numeric(),
  exacerbated_pct = numeric(),
  stringsAsFactors = FALSE
)

for (cl in row.names(treatment_summary)) {
  # Calculate percentages
  total_deg <- treatment_summary[cl, "vcp_up_total"] + treatment_summary[cl, "vcp_down_total"]
  normalized <- treatment_summary[cl, "vcp_up_normalized"] + treatment_summary[cl, "vcp_down_normalized"]
  exacerbated <- treatment_summary[cl, "vcp_up_exacerbated"] + treatment_summary[cl, "vcp_down_exacerbated"]
  unaffected <- treatment_summary[cl, "vcp_up_unaffected"] + treatment_summary[cl, "vcp_down_unaffected"]
  
  if (total_deg > 0) {  # Avoid division by zero
    # Add to summary data
    summary_data <- rbind(summary_data, data.frame(
      cluster = cl,
      normalized_pct = normalized / total_deg * 100,
      no_change_pct = unaffected / total_deg * 100,
      exacerbated_pct = exacerbated / total_deg * 100
    ))
  }
}

# Create the summary plot
if (nrow(summary_data) > 0) {
  library(reshape2)
  summary_long <- melt(summary_data, id.vars = "cluster", 
                       variable.name = "effect", value.name = "percentage")
  summary_long$effect <- gsub("_pct", "", summary_long$effect)
  summary_long$effect <- factor(summary_long$effect, 
                                levels = c("normalized", "no_change", "exacerbated"),
                                labels = c("Normalized", "No Change", "Exacerbated"))
  
  summary_plot <- ggplot(summary_long, aes(x = reorder(cluster, percentage, 
                                                       function(x) sum(x[summary_long$effect[summary_long$cluster == 
                                                                                               unique(summary_long$cluster)[which(unique(summary_long$cluster) == 
                                                                                                                                    as.character(unique(cluster)))]] == "Normalized"])), 
                                           y = percentage, fill = effect)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Normalized" = "forestgreen", 
                                 "No Change" = "gray60", 
                                 "Exacerbated" = "firebrick")) +
    labs(
      title = "SB Treatment Effect Across All Clusters",
      subtitle = "Percentage of VCP-dysregulated genes that are normalized by SB treatment",
      x = "Cell Cluster",
      y = "Percentage of Dysregulated Genes",
      fill = "Effect of SB Treatment"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Save the summary plot
  ggsave(paste0(main_dir, script_ind, "treatment_effect_summary.pdf"), 
         summary_plot, width = 10, height = 7)
  
  # Print the plot
  print(summary_plot)
}

message("\n\n          *** All analyses complete! ", Sys.time(),"\n\n")

# Save the updated bulk_data object
save(bulk_data, file = paste0(main_dir, script_ind, "bulk_data_with_treatment_analysis_treatment_analysis.rda"))







