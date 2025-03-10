library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

TrioBrain.integrated_slim <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim.rds")
TrioBrain.integrated_slim_ISub_DF_labeled <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_ISub_DF_labeled.rds")

nodoublet_id <- WhichCells(TrioBrain.integrated_slim_ISub_DF_labeled, idents='Doublets', invert=T)

load("Processed_Data/cell_type_lists.RData")

TrioBrain.integrated_slim_labeled_final <- TrioBrain.integrated_slim

label_transfer_integrated_seurat <- function(subcluster_seurat=list_subcluster_celltypes_glia, integrated_seurat=TrioBrain.integrated_slim_labeled_final) {
  
  subcluster_idents <- as.character(unique(Idents(subcluster_seurat[[5]])))
  
  for (i in 1:length(subcluster_idents)) {
    
    sub_id <- subcluster_idents[i]
    cell_id <- WhichCells(subcluster_seurat[[5]], idents=sub_id)
    Idents(integrated_seurat, cells = cell_id ) <- sub_id
    
  }
  
  return(integrated_seurat)
  
}

subcluster_list <- list(list_subcluster_celltypes_glia, list_subcluster_celltypes_KC, list_subcluster_celltypes_MA, 
  list_subcluster_celltypes_Clock, list_subcluster_celltypes_Poxn, 
  list_subcluster_celltypes_OPN,
  list_subcluster_celltypes_fru, list_subcluster_celltypes_NP, list_subcluster_celltypes_Ach, 
  list_subcluster_celltypes_Glu, list_subcluster_celltypes_GABA)

for (nl in 1:length(subcluster_list)) {
  
  TrioBrain.integrated_slim_labeled_final <- label_transfer_integrated_seurat(subcluster_seurat=subcluster_list[[nl]], integrated_seurat=TrioBrain.integrated_slim_labeled_final)
  
}

DimPlot(TrioBrain.integrated_slim_labeled_final, label=T, reduction='tsne') + NoAxes() + NoLegend()

saveRDS(TrioBrain.integrated_slim_labeled_final, file = "Processed_Data/TrioBrain.integrated_slim_labeled_final.rds")

Final_idents <- as.character(unique(Idents(TrioBrain.integrated_slim_labeled_final)))
exclude <- c(Final_idents[grep('Unanno', Final_idents)], 0:100)

Anno_idents <- setdiff(Final_idents,exclude)

TrioBrain.integrated_slim_labeled_final_anno_only <- subset(TrioBrain.integrated_slim_labeled_final, idents = Anno_idents)

DimPlot(TrioBrain.integrated_slim_labeled_final_anno_only, label=T, reduction='tsne') + NoAxes() + NoLegend()
saveRDS(TrioBrain.integrated_slim_labeled_final_anno_only, file = "Processed_Data/TrioBrain.integrated_slim_labeled_final_anno_only.rds")

## DimPlots

p1 <- DimPlot(TrioBrain.integrated_slim_labeled_final, reduction = "tsne", label = TRUE, raster=FALSE, pt.size = 0.01) + NoLegend() + NoAxes()
p2 <- DimPlot(TrioBrain.integrated_slim_labeled_final_anno_only, reduction = "tsne", label = TRUE, raster=FALSE, pt.size = 0.01) + NoLegend() + NoAxes()
plot_grid(p1, p2)

ggsave("Plots/TrioBrain.integrated_slim_labeled_final_tsne.png", width=17, height = 8)

