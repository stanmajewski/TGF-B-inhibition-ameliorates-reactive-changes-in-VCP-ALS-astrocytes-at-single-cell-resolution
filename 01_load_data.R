message("\n\n##########################################################################\n",
        "# Start 01: Load dataset into Seurat ", Sys.time(),
        "\n##########################################################################\n",
        "\n   Plot cellranger QC, load count data, merge datasets, plot Seurat QC, \n",
        "\n##########################################################################\n\n")

#set environment and main directory

main_dir = getwd()
setwd(main_dir)

# Open packages necessary for analysis.
library(tidyverse)
library(Seurat)
library(patchwork)


#specify script/output index as prefix for file names
script_ind = "01_"

#specify output directory
out_dir = paste0(main_dir, "/results/")

# Create output directory if it doesn't exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


#load group and file info for analysis dataset
gr_tab_data = read_csv("data/metadata/sample_metadata.csv")

#######################################
#load cellranger QC data
#######################################

message("\n\n     *** Loading/plotting Cellranger QC data (", Sys.time(), ") \n\n")

set.seed(1234)

#load cellranger QC data

cellranger_QC = NULL
i = 1
for (i in 1:nrow(gr_tab_data)){
  
  t1 = read_csv(paste0(gr_tab_data$dataset_folder[i],
                       gr_tab_data$geo_accession[i],"/outs/metrics_summary.csv"))
  t2 = cbind(sample = gr_tab_data$sample[i], 
             geo_accession = gr_tab_data$geo_accession[i], 
             t1)
  cellranger_QC = rbind(cellranger_QC, t2)
  
}

write_csv(cellranger_QC, paste0(out_dir,script_ind,"cellranger_QC.csv"))


### plot cellranger QC

t1 = cellranger_QC
names(t1)

plot_cols = c("Estimated Number of Cells", "Mean Reads per Cell", 
              "Median UMI Counts per Cell","Median Genes per Cell",
              "Fraction Reads in Cells",
              "Reads Mapped Confidently to Genome", 
              "Reads Mapped Confidently to Transcriptome")


pl =  lapply(plot_cols, function(plot_col){
  t2 = t1[,c('sample', plot_col)]
  names(t2) = c("sample", "pl_col")
  t2$pl_col = as.numeric(str_remove(t2$pl_col, "%"))
  
  g1 = ggplot(t2, aes(x = sample, y = as.numeric(pl_col)))+geom_col()+
    theme_classic()+ labs(title = plot_col, y = plot_col)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
})           


pdf(file = paste0(out_dir,script_ind,"cellranger_QC_barplots.pdf"), 
    width = 10, height = 5)
lapply(pl, function(x){x})
dev.off()




#######################################
#load RNA, ATAC data, create Seurat objects and load cellranger QC data
#######################################

#load RNA and ATAC data
count_list = list()
seur_list = list()

for (i in 1:nrow(gr_tab_data)){
  
  message("\n\n     ***Loading count data for ", gr_tab_data$sample[i], 
          " (sample ",i, " of ", nrow(gr_tab_data)," )\n\n")
  
  counts = 
    Read10X(paste0(gr_tab_data$dataset_folder[i]))
  
  # create a Seurat object containing the RNA data
  seur_temp <- CreateSeuratObject(
    counts = CreateAssayObject(counts = counts),
    assay = "RNA"
  )
  
  #add sample metadata
  for (meta1 in colnames(gr_tab_data)){
    seur_temp[[meta1]] = gr_tab_data[[meta1]][i]
  }
  
  #add percent.mt QC
  v1 = seur_temp@assays$RNA@counts@Dimnames[[1]]
  v2 = v1[grepl(v1, pattern = "^MT-")]
  seur_temp[["percent.mt"]] <- PercentageFeatureSet(seur_temp, features = v2, assay = "RNA")
  
  #rename cells to make names unique to allow merging
  seur_temp <- RenameCells(object = seur_temp, add.cell.id = gr_tab_data$sample[i])
  
  seur_list[[ gr_tab_data$sample[i] ]] = seur_temp

}

#checkpoint
gc(verbose = TRUE, reset = FALSE, full = TRUE)
message("\n\n     *** Count data loaded. Saving Seurat objects... (", Sys.time(), ") \n\n")

save(seur_list, 
     file = paste0(out_dir, "seur_list_patient.rda")) 

gc(verbose = TRUE, reset = FALSE, full = TRUE)



#########################################
# Merge datasets 
#########################################

message("\n\n     *** Seurat objects saved. Merging Seurat objects... (", Sys.time(), ") \n\n")

#merge datasets preliminarily for QC plots

seur <- merge(seur_list[[1]], y = seur_list[-1], merge.data = TRUE)

#rm(seur_list) #free up space
gc()

message("\n\n     *** Seurat objects merged. Saving combined dataset... (", Sys.time(), ") \n\n")

save(seur, file = paste0(out_dir,script_ind,"seur_merge_patient.rda")) 



#########################################
# QC plots for merged dataset 
#########################################

message("\n\n     *** Combined dataset saved. Plotting Seurat QC... (", Sys.time(), ") \n\n")


###QC plots

#for large datasets, use subsample for plotting (else plots get  too large)
if (nrow(seur@meta.data)>50000){
  seur = seur[, sample(colnames(seur), size =50000, replace=F)]
} 

meta = seur@meta.data

pl = list()

pl[["nCount_RNA_max50k"]] = VlnPlot(seur, features = "nCount_RNA", group.by = "sample", pt.size = 0.01) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  NoLegend()+ ylim(0, 50000)+ labs(title = "nCount_RNA_max50k") 
pl[["nCount_RNA_max3000"]] = VlnPlot(seur, features = "nCount_RNA", group.by = "sample", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  +
  NoLegend()+ ylim(0, 3000)+labs(title = "nCount_RNA_max3k")
pl[["nFeature_RNA"]] = VlnPlot(seur, features = "nFeature_RNA", group.by = "sample", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  + 
  NoLegend()+ labs(title = "nFeature_RNA")
pl[["percent.mt"]] = VlnPlot(seur, features = "percent.mt", group.by = "sample", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  + 
  NoLegend()+ labs(title = "percent.mt")
pl[["percent.mt_max5"]] = VlnPlot(seur, features = "percent.mt", group.by = "sample", pt.size = 0.01)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  NoLegend()+ylim(0,5)+ labs(title = "percent.mt_max5%")


pdf(file = paste0(out_dir,script_ind,"QC_plots.pdf"), width = 10, height = 6)
lapply(pl, function(x){x} )
dev.off()


message("\n\n##########################################################################\n",
        "# Completed 01 ", Sys.time(),
        "\n##########################################################################\n",
        "\n##########################################################################\n\n\n")

