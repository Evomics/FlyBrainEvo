library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(sctransform)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/shared_genes.RData")
TrioBrain.integrated_slim<- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim.rds")

ISub <- function(Seurat_object, cluster, method='rpca', resolution = 0.1, npc=50, genelist=shared_genes, kweight=100) {
  
  Seurat <- subset(Seurat_object, idents = cluster)
  
  DefaultAssay(Seurat) <- "RNA"
  
  Dmel <- subset(Seurat, subset = orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"))
  Dsim <- subset(Seurat, subset = orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"))
  Dsec <- subset(Seurat, subset = orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"))
  DsecNoni <- subset(Seurat, subset = orig.ident %in% c("DsecNoni_rep1", "DsecNoni_rep2", "DsecNoni_rep3", "DsecNoni_rep4", "DsecNoni_rep5", "DsecNoni_rep6"))  
  
  Seurat_list <- list(Dmel, Dsim, Dsec, DsecNoni) 
  
  Seurat_list <- lapply(X = Seurat_list, FUN = SCTransform, method = "glmGamPoi", verbose = FALSE)
  
  Int_features <- SelectIntegrationFeatures(object.list = Seurat_list, nfeatures = 3000, verbose = FALSE)
  
  Seurat_list <- PrepSCTIntegration(object.list = Seurat_list, assay = "SCT", 
                                    anchor.features = Int_features,
                                    verbose = F)
  
  Seurat_list <- lapply(X = Seurat_list, FUN = RunPCA, features = Int_features, verbose = FALSE)
  
  anchors<- FindIntegrationAnchors(object.list = Seurat_list, normalization.method = "SCT",
                                   anchor.features = Int_features, verbose = F, reduction = method, reference=1)
  
  Seurat_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = F,
                                     features.to.integrate = genelist,
                                     k.weight = kweight)
  
  Seurat_integrated <- RunPCA(Seurat_integrated, npcs = npc, verbose = F)
  Seurat_integrated <- FindNeighbors(Seurat_integrated, reduction = "pca", dims = 1:npc, verbose = F)
  Seurat_integrated <- FindClusters(Seurat_integrated, resolution = resolution, verbose = F)
  
  DefaultAssay(Seurat_integrated) <- 'RNA'
  Seurat_integrated <- DietSeurat(Seurat_integrated, assays = c("RNA"))
  
  return(Seurat_integrated)
  
}

Cluster14 <- ISub(TrioBrain.integrated_slim, 14)

Cluster14 <- SCTransform(Cluster14, method = "glmGamPoi", verbose = FALSE)
Cluster14 <- RunPCA(Cluster14, npcs = 50, verbose = F)
Cluster14 <- RunUMAP(Cluster14, dims = 1:50, verbose = F)
  
FeaturePlot(Cluster14, feature = 'Ilp6', order=T)
DotPlot(Cluster14, feature = 'Ilp6', group.by='orig.ident')

Dmel_metadata <- Trio_ISub_DF_list_labeled[[1]]@meta.data

Ilp6_cells <- WhichCells(Trio_ISub_DF_list_labeled[[1]], idents = 'Ilp6')
DimPlot(Trio_ISub_DF_list_labeled[[1]], cells.highlight=Ilp6_cells, reduction='tsne') + NoLegend()

Dmel_markers <- Trio_ISub_DF_marker_list[[1]]

### Ilp6 population is 'Fat' cell type, only identified in Dmel ###

