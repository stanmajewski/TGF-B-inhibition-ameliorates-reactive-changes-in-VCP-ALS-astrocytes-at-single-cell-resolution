message("\n\n##########################################################################\n",
        "# Start 06 Pseudobulk DESeq2 analysis: ", Sys.time(),
        "\n##########################################################################\n",
        "\n   contrasts to extract DEGs by cluster",
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
library(DESeq2)

#specify script/output index as prefix for file names
script_ind = "06_"

#specify output directory
out_dir = paste0(main_dir, "/results/")

#load group and file info
gr_tab = read_csv("results/02_gr_tab_filtered.csv")

#load dataset
load(file = paste0(main_dir,"results/05_bulk_data.rda")) 

#correct cell names in count matrix
#t1 = colnames(bulk_data$counts)
#t2 = str_replace_all(t1, "snRNA.", "snRNA-")
#colnames(bulk_data$counts) = t2

#get marker gene panels
GOI = list()
t1 = read_csv(paste0(main_dir,"reference_data/transcription_factors.csv"))
GOI$TF = t1$Symbol


#define group comparisons for each cluster
comp_groups = list(VCP.u_vs_ctrl.u = c("VCP_u","ctrl_u"),
                   ctrl.SB_vs_ctrl.u = c("ctrl_SB","ctrl_u"),
                   VCP.SB_vs_VCP.u = c("VCP_SB", "VCP_u"),
                   VCP.SB_vs_ctrl.u = c("VCP_SB", "ctrl_u"))


###########################################################
# functions
###########################################################


###########################################################
# remove clusters with each <2 pseudobulks per cluster per group
###########################################################

message("\n\n          *** Preparing dataset for DESeq2 analysis... ", Sys.time(),"\n\n")

set.seed(12345)

t1 = bulk_data$meta
t2 = t1 %>% group_by(cluster_name, group) %>% summarise(N_pseudobulks = n())

comp_clusters = NULL

for (cl in unique(t1$cluster_name)){
  t3 = t2[t2$cluster_name == cl,]
  if (nrow(t3)==length(unique(t1$group))){
    if(min(t3$N_pseudobulks >=2)){
      comp_clusters = c(comp_clusters, cl)
    }
  }
}

comp_meta = t1[t1$cluster_name %in% comp_clusters,]

# define grouping variable to extract comparisons by group for each cluster
comp_meta$cluster_group = paste0(comp_meta$cluster_name, "_", comp_meta$group)

rownames(comp_meta) = comp_meta$pseudobulk

bulk_data$meta = comp_meta


#remove genes from countmatrix with <0.1 counts/cell in all pseudobulks  

pseudobulks = unique(comp_meta$cluster_sample)

# Add X prefix to match count matrix column names
pseudobulks <- paste0("X", unique(comp_meta$cluster_sample))

#m1 = bulk_data$counts[,pseudobulks]

N_pb_cells = lengths(bulk_data$cells)

#for (pb in pseudobulks){
#  m1[,pb] = m1[,pb]/N_pb_cells[pb]
#}

#keep_genes = rownames(m1)[apply(m1, 1, max)>0.1]

#comp_counts = bulk_data$counts[keep_genes ,comp_meta$pseudobulk]

# Start fresh with m1
m1 <- bulk_data$counts[, pseudobulks]

# Fix the normalization by using pseudobulks without X prefix for N_pb_cells lookup
for (pb in pseudobulks){
  pb_no_prefix <- gsub("^X", "", pb)  # Remove X prefix to match N_pb_cells names
  m1[,pb] <- m1[,pb] / N_pb_cells[pb_no_prefix]
}

# Now filter genes properly
keep_genes <- rownames(m1)[apply(m1, 1, max) > 0.1]

print("Number of genes passing filter:")
print(length(keep_genes))

# Create final count matrix
comp_counts <- bulk_data$counts[keep_genes, pseudobulks]

print("Final comp_counts dimensions:")
print(dim(comp_counts))

#Set rownames properly using the pseudobulk/cluster_sample column with X prefix
comp_meta <- as.data.frame(comp_meta)
rownames(comp_meta) <- paste0("X", comp_meta$cluster_sample)


#################################################
# DEG analysis for count data with DESeq2
#################################################

#run DESeq2 analysis (LRT test to assess full model vs reduced model to identify sign changes)

message("\n\n          *** Running DESeq2 analysis... ", Sys.time(),"\n\n")


dds = DESeqDataSetFromMatrix(comp_counts, colData = comp_meta, 
                             design = ~cluster_group)
dds = DESeq(dds) #standard Wald test

bulk_data$deseq_dataset = dds

bulk_data$vst_mat = assay(vst(dds))


#save bulk_dataset with DESeq results

save(bulk_data, file = paste0(main_dir,script_ind, "bulk_data_with_DESeq_results.rda")) 


#################################################
# Extract DESeq2 results
#################################################

message("\n\n          *** Extracting DESeq2 results... ", Sys.time(),"\n\n")

#extract stats for each comparison (genes signif dereg and with tendency of dereg)
deseq_res_list = NULL
genes_reg_list = NULL

for (cl in comp_clusters){
  message("  * extracting DESeq results for ", cl)
  
  for (comp_group in names(comp_groups)){
    comp1 = paste0(cl, "_", comp_groups[[comp_group]][1])
    comp2 = paste0(cl, "_", comp_groups[[comp_group]][2])
    t1 = as.data.frame(results(dds, contrast=c("cluster_group", comp1, comp2)))
    
    deseq_res_list[[paste0(cl, "_", comp_group)]] = t1
    
    v3 = rownames(t1)[!is.na(t1$padj) & 
                        t1$padj <= 0.05 & t1$log2FoldChange >= +log2(1.2)]
    genes_reg_list[[paste0(cl, "_", comp_group,"_up")]] = v3
    v3 = rownames(t1)[!is.na(t1$padj) & 
                        t1$padj <= 0.05 & t1$log2FoldChange <= -log2(1.2)]
    genes_reg_list[[paste0(cl, "_", comp_group,"_down")]] = v3
  }
  
}

bulk_data$deseq_results = deseq_res_list
bulk_data$DEGs = genes_reg_list


#save bulk_dataset with DESeq results

save(bulk_data, file = paste0(main_dir, "06_bulk_data_with_DESeq_results.rda")) 



# save table with DEGs by comparison including TFs and HSA21 genes

l1 = genes_reg_list
l2 = lapply(l1, function(x){x = x[x%in% GOI$TF]})
names(l2) = paste0(names(l1), "_TF")

l3 = c(l1, l2)

m1 = matrix(nrow = max(lengths(l3)), ncol = length(l3))
colnames(m1) = names(l3)

for (i in names(l3)){
  v1 = l3[[i]]
  if(length(v1)>0){
    m1[1:length(v1),i] = v1
  }
}
m1[is.na(m1)] = ""

write_csv(as_tibble(m1), file = paste0(main_dir,script_ind, "06_DESeq2_up_down_by_cluster_genes.csv"))

t1 = tibble(gene_set = names(l3), N_genes = lengths(l3))

write_csv(t1, file = paste0(main_dir,script_ind, "06_DESeq2_up_down_by_cluster_N_genes.csv"))


message("\n\n##########################################################################\n",
        "# Completed 06 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")






