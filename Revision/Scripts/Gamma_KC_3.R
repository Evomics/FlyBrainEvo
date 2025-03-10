library(ggplot2)
library(dplyr)
library(ggsankey)
library(Seurat)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/shared_genes.RData")

TrioBrain.integrated_slim_ISub_DF_labeled <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_ISub_DF_labeled.rds")
Trio_ISub_DF_list_labeled <- readRDS(file = "Processed_Data/Trio_ISub_DF_list_labeled.rds")

Dmel_gamma_KC <- subset(Trio_ISub_DF_list_labeled[[1]], idents = c('γ-KC_1','γ-KC_2'))

DimPlot(Dmel_gamma_KC)+ NoAxes()
DimPlot(Dmel_gamma_KC, group.by='orig.ident')+ NoAxes()

DsecNoni_gamma_KC <- subset(Trio_ISub_DF_list_labeled[[4]], idents = c('γ-KC'))

DimPlot(DsecNoni_gamma_KC)+ NoAxes()
DimPlot(DsecNoni_gamma_KC, group.by='orig.ident')+ NoAxes()

Dsec_gamma_KC <- subset(Trio_ISub_DF_list_labeled[[3]], idents = c('γ-KC_1','γ-KC_2'))

DimPlot(Dsec_gamma_KC)+ NoAxes()
DimPlot(Dsec_gamma_KC, group.by='orig.ident')+ NoAxes()




Dsec_gamma_KC_markers <- FindMarkers(Dsec_gamma_KC, only.pos = F, min.pct = 0.15, 
                                                              logfc.threshold = 0.25)

Dsec_gamma_KC_markers_top10 <- Dsec_gamma_KC_markers %>%
  dplyr::filter(cluster == "γ-KC_2") %>%
  top_n(n = 10, wt = -p_val_adj) %>%
  top_n(n = 10, wt = avg_log2FC)

FeaturePlot(Dsec_gamma_KC, feature=Dsec_gamma_KC_markers_top10$gene) * NoAxes()

DotPlot(Dsec_gamma_KC, feature=Dsec_gamma_KC_markers_top10$gene)



DotPlot(TrioBrain.integrated_slim_ISub_DF_labeled, feature=c('14-3-3zeta', 'slo', 'amon', 'CG42268', 'mmd', 'mamo', 'ctp', 'Tomosyn'))

DotPlot(subset(TrioBrain.integrated_slim_ISub_DF_labeled, idents = 'AST'), feature=c('14-3-3zeta', 'slo', 'amon', 'CG42268', 'mmd', 'mamo', 'ctp', 'Tomosyn'))



TrioBrain_ISub_DF.integrated.markers.rds <- readRDS(file="Processed_Data/TrioBrain_ISub_DF.integrated.markers.rds")

DimPlot(TrioBrain.integrated_slim_ISub_DF_labeled, reduction='tsne', label=T) + NoLegend() + NoAxes()

metadata_integrated <- TrioBrain.integrated_slim_ISub_DF_labeled@meta.data

gamma_KC_subset <- subset(TrioBrain.integrated_slim_ISub_DF_labeled, idents = c('γ-KC_1','γ-KC_2','γ-KC_3','γ-KC_4'))
gamma_KC_subset@meta.data$orig.ident <- gsub("_rep[1-6]", "", gamma_KC_subset@meta.data$orig.ident)

gamma_KC_subset@meta.data$orig.ident <- factor(gamma_KC_subset@meta.data$orig.ident, 
                                                                         levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))

DimPlot(gamma_KC_subset, reduction='tsne', label=T, group.by = 'orig.ident') + NoAxes()

TrioBrain.integrated_slim <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim.rds")
DefaultAssay(TrioBrain.integrated_slim) <- 'RNA'
TrioBrain.integrated_slimer <- DietSeurat(TrioBrain.integrated_slim, assays = c("RNA"))

Get_subcluster <- function(Seurat_object, cluster, method='rpca', resolution = 0.1, npc=50, genelist=shared_genes, kweight=100) {
  
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
  Seurat_integrated <- RunUMAP(Seurat_integrated, dims = 1:npc, verbose = F)
  Seurat_integrated <- RunTSNE(Seurat_integrated, dims = 1:npc, verbose = F)

  return(Seurat_integrated)
  
}

gamma_KC <- Get_subcluster(TrioBrain.integrated_slimer, cluster=6)

DimPlot(gamma_KC, reduction='tsne', label=T) + NoLegend() + NoAxes()

gamma_KC_0 <- Get_subcluster(gamma_KC, cluster=0, kweight=60)

DimPlot(gamma_KC_0, reduction='tsne', label=T) + NoLegend() + NoAxes()

gamma_KC_0_1 <- Get_subcluster(gamma_KC_0, cluster=1, kweight=100)

DimPlot(gamma_KC_0_1, reduction='tsne', label=T) + NoLegend() + NoAxes()

gamma_KC_0_1@meta.data$orig.ident <- gsub("_rep[1-6]", "", gamma_KC_0_1@meta.data$orig.ident)

DimPlot(gamma_KC_0_1, reduction='tsne', label=T, group.by = 'orig.ident') + NoLegend() + NoAxes()
