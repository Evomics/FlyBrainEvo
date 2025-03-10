library(Seurat)
library(tidyverse)
library(ggplot2)

TrioBrain.integrated_slim_labeled_final <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_labeled_final.rds")

TrioBrain.integrated_slim_labeled_final_species <- TrioBrain.integrated_slim_labeled_final
TrioBrain.integrated_slim_labeled_final_species$orig.ident <- gsub("_rep[1-6]", "", TrioBrain.integrated_slim_labeled_final_species$orig.ident)
TrioBrain.integrated_slim_labeled_final_species@meta.data$orig.ident <- factor(TrioBrain.integrated_slim_labeled_final_species@meta.data$orig.ident, 
                                                                               levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))

df <- subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident != "DsecNoni")




DimPlot(df, reduction = "tsne", shuffle = T, raster=F, pt.size = 0.01, label=T) + NoLegend() + labs(title=NULL) + NoAxes() + 
  coord_cartesian(xlim=c(-15,-5),ylim=c(-7,3))

DimPlot(df, 
        reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F,
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8")) + 
  coord_cartesian(xlim=c(-15,-5),ylim=c(-7,5))

DimPlot(subset(df, idents=0), reduction = "tsne") + 
  coord_cartesian(xlim=c(-15,-5),ylim=c(-7,5))

DimPlot(subset(df, idents=0), reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F,
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8"))+ labs(title=NULL)

ggsave("Plots/Manuscript/Cluster0_tsne.png", width=7, height = 6)

DimPlot(subset(df, idents=0), 
        reduction = "umap", group.by = "orig.ident", shuffle = T, raster=F,
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8"))+ labs(title=NULL)

ggsave("Plots/Manuscript/Cluster0_umap.png", width=7, height = 6)


DimPlot(df, 
        reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F,
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8"))

