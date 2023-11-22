library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(sctransform)
library(cowplot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/shared_genes.RData")

#TrioBrain.integrated_slim <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim.rds")
TrioBrain.integrated_slim_labeled<- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled.rds")
TrioBrain.integrated.markers <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated.markers.rds")

Subclustering <- function(Seurat_object, cluster, method='rpca', resolution = 0.8, npc=50, genelist=shared_genes, ptsize=0.1) {
  
  Seurat <- subset(Seurat_object, idents = cluster)

  DefaultAssay(Seurat) <- "RNA"
  
  Dmel <- subset(Seurat, subset = orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"))
  Dsim <- subset(Seurat, subset = orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"))
  Dsec <- subset(Seurat, subset = orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"))
  DsecNoni <- subset(Seurat, subset = orig.ident %in% c("DsecNoni_rep1", "DsecNoni_rep2", "DsecNoni_rep3", "DsecNoni_rep4", "DsecNoni_rep5", "DsecNoni_rep6"))  
  
  Seurat_list <- list(Dmel, Dsim, Dsec, DsecNoni) 
  
  Seurat_list <- lapply(X = Seurat_list, FUN = SCTransform, method = "glmGamPoi", verbose = FALSE)
   
  Int_features <- SelectIntegrationFeatures(object.list = Seurat_list, nfeatures = 3000)
  
  Seurat_list <- PrepSCTIntegration(object.list = Seurat_list, assay = "SCT", 
                                                       anchor.features = Int_features,
                                                       verbose = TRUE)
  
  Seurat_list <- lapply(X = Seurat_list, FUN = RunPCA, features = Int_features)
  
  anchors<- FindIntegrationAnchors(object.list = Seurat_list, normalization.method = "SCT",
                                                             anchor.features = Int_features, verbose = T, reduction = method, reference=1)
  
  Seurat_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = TRUE,
                                                features.to.integrate = genelist)
  
  Seurat_integrated <- RunPCA(Seurat_integrated, npcs = npc, verbose = T)
  
  Seurat_integrated <- RunUMAP(Seurat_integrated, dims = 1:npc)
  Seurat_integrated <- RunTSNE(Seurat_integrated, dims = 1:npc)
  Seurat_integrated <- FindNeighbors(Seurat_integrated, reduction = "pca", dims = 1:npc)
  Seurat_integrated <- FindClusters(Seurat_integrated, resolution = resolution)
  
  Seurat_integrated_species <- Seurat_integrated
  
  Seurat_integrated_species@meta.data$orig.ident <- gsub("_rep1", "", Seurat_integrated_species@meta.data$orig.ident)
  Seurat_integrated_species@meta.data$orig.ident <- gsub("_rep2", "", Seurat_integrated_species@meta.data$orig.ident)
  Seurat_integrated_species@meta.data$orig.ident <- gsub("_rep3", "", Seurat_integrated_species@meta.data$orig.ident)
  Seurat_integrated_species@meta.data$orig.ident <- gsub("_rep4", "", Seurat_integrated_species@meta.data$orig.ident)
  Seurat_integrated_species@meta.data$orig.ident <- gsub("_rep5", "", Seurat_integrated_species@meta.data$orig.ident)
  Seurat_integrated_species@meta.data$orig.ident <- gsub("_rep6", "", Seurat_integrated_species@meta.data$orig.ident)
  
  Seurat_integrated_species@meta.data$orig.ident <- factor(Seurat_integrated_species@meta.data$orig.ident, 
                                                                   levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))
  
  plot1 <- DimPlot(Seurat_integrated_species, reduction = "umap", 
                     group.by = "orig.ident", shuffle = T, pt.size = ptsize, raster=FALSE) + labs(title=NULL)
  plot2 <- DimPlot(Seurat_integrated_species, reduction = "umap", label = TRUE, raster=FALSE, pt.size = ptsize)
  
  plot <- plot_grid(plot1, plot2)
  
  print(plot)
  
    df_Seurat_integrated_metadata <- Seurat_integrated@meta.data %>%
    dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                                 ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                        ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                               "DsecNoni")))) %>%
    dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))
  
  
  df_Seurat_integrated_metadata_summary_rep <- df_Seurat_integrated_metadata %>%
    dplyr::group_by(species, seurat_clusters, orig.ident) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(species, orig.ident) %>%
    dplyr::mutate(total=sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(percent = n/total*100)
  
  df_Seurat_integrated_metadata_summary <- df_Seurat_integrated_metadata_summary_rep %>%
    dplyr::group_by(species, seurat_clusters) %>%
    dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
    dplyr::ungroup()
  
  plot_quant <- df_Seurat_integrated_metadata_summary %>%
    ggplot(.) +
    geom_col(aes(y=seurat_clusters, x=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.9)) + 
    geom_errorbar(aes(y=seurat_clusters, x=percent_combined,
                      xmin=ifelse(percent_combined-sem<0,0,percent_combined-sem), xmax=percent_combined+sem, group=species), width=0.9, position=position_dodge(width=0.9)) +
    geom_point(data=df_Seurat_integrated_metadata_summary_rep, size = 1, alpha =0.8,  
               aes(y=seurat_clusters, x=percent, group=species), fill='red', shape=21, position=position_dodge(width=0.9)) +
    #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=seurat_clusters, y=n/total*100)) +
    theme_bw() +
    theme(axis.title.x = element_text(size=11, color='black'),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=10, color='black'),
          axis.text.y=element_text(size=10, color='black'),
          panel.grid = element_blank()) +
    labs(y="Cluster", x="Percent of cell type (%)", fill="") +
    scale_fill_brewer(palette = 'Set1')
  
  print(plot_quant)
  
  return_list <- list(Seurat_integrated, plot, plot_quant, df_Seurat_integrated_metadata_summary, df_Seurat_integrated_metadata_summary_rep)
  
  return(return_list)
  
}

## IPC+Dh44 ##

Cluster17_res1 <- Subclustering(Seurat_object=Cluster17_seurat, cluster = 0:11, resolution = 1, method='rpca', npc=10)

DefaultAssay(Cluster17_res1[[1]]) <- "SCT"

FeaturePlot(Cluster17_res1[[1]], feature=c("Dh44", "Ilp2"), pt.size=0.3, min.cutoff=3)

Cluster17_res1[[2]]
Cluster17_res1[[3]]
View(Cluster17_res1[[4]])

Cluster17_res1_seurat <- DietSeurat(Cluster17_res1[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))

Cluster17_res1_seurat <- PrepSCTFindMarkers(Cluster17_res1_seurat, assay = "SCT", verbose = TRUE)

Cluster17_res1_markers <- FindAllMarkers(Cluster17_res1_seurat, only.pos = F, min.pct = 0.15, 
                               logfc.threshold = 0.25, test.use = "MAST")


top24_Cluster17_res1 <- Cluster17_res1_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster17_res1$cluster)) {
  
  top24_Cluster17_res1_sub <- top24_Cluster17_res1 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster17_res1_seurat, features = top24_Cluster17_res1_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster17_res1_subclusters{i}.png"), width=24, height = 24)
  
}

## Dh44 ##

Dh44 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'Dh44', resolution = 0.1, method='cca')

FeaturePlot(Dh44[[1]], feature='Dh44')

## subset large cluster ##

Cluster0 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'Cluster0', resolution = 0.1)

DefaultAssay(Cluster0[[1]]) <- "SCT"

Cluster0_seurat <- DietSeurat(Cluster0[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))

DimPlot(Cluster0_seurat, label=T, repel=T)

Cluster0_seurat <- PrepSCTFindMarkers(Cluster0_seurat, assay = "SCT", verbose = TRUE)

Cluster0_markers <- FindAllMarkers(Cluster0_seurat, only.pos = F, min.pct = 0.15, 
                               logfc.threshold = 0.25, test.use = "MAST")

FeaturePlot(Cluster0_seurat, features=c("FMRFa"))
VlnPlot(Cluster0_seurat, features=c("FMRFa"))
FeaturePlot(TrioBrain.integrated_slim_labeled, features=c("FMRFa"))

#save(Cluster0_markers, Cluster0_seurat,file="Processed_Data/Cluster0.RData")
load("Processed_Data/Cluster0.RData")

## Cluster0- subclustesr markers ##

top24_Cluster0 <- Cluster0_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster0$cluster)) {
  
  top24_Cluster0_sub <- top24_Cluster0 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster0_seurat, features = top24_Cluster0_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster0_subclusters{i}.png"), width=24, height = 24)
  
}

## KD intersection for annotation

## KD (2018) annotation

Cluster0_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=celltype) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~cluster, scale='free', ncol=6)

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50_Cluster0.png"), width=20, height =20)

## frequency ##

Seurat_integrated <- Cluster0_seurat

df_Seurat_integrated_metadata <- Seurat_integrated@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))


df_Seurat_integrated_metadata_summary_rep <- df_Seurat_integrated_metadata %>%
  dplyr::group_by(species, seurat_clusters, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

df_Seurat_integrated_metadata_summary <- df_Seurat_integrated_metadata_summary_rep %>%
  dplyr::group_by(species, seurat_clusters) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

plot_quant <- df_Seurat_integrated_metadata_summary %>%
  ggplot(.) +
  geom_col(aes(y=seurat_clusters, x=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.9)) + 
  geom_errorbar(aes(y=seurat_clusters, x=percent_combined,
                    xmin=ifelse(percent_combined-sem<0,0,percent_combined-sem), xmax=percent_combined+sem, group=species), width=0.9, position=position_dodge(width=0.9)) +
  geom_point(data=df_Seurat_integrated_metadata_summary_rep, size = 1, alpha =0.8,  
             aes(y=seurat_clusters, x=percent, group=species), fill='red', shape=21, position=position_dodge(width=0.9)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=seurat_clusters, y=n/total*100)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11, color='black'),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=10, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        panel.grid = element_blank()) +
  labs(y="Cluster", x="Percent of cell type (%)", fill="") +
  scale_fill_brewer(palette = 'Set1')

plot_quant

## Cluster 1 ##

Cluster1 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'Cluster1', resolution = 0.1)

DefaultAssay(Cluster1[[1]]) <- "SCT"

Cluster1_seurat <- DietSeurat(Cluster1[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))

DimPlot(Cluster1_seurat, label=T, repel=T)

Cluster1_seurat <- PrepSCTFindMarkers(Cluster1_seurat, assay = "SCT", verbose = TRUE)

Cluster1_markers <- FindAllMarkers(Cluster1_seurat, only.pos = F, min.pct = 0.15, 
                                   logfc.threshold = 0.25, test.use = "MAST")


save(Cluster1_markers, Cluster1_seurat,file="Processed_Data/Cluster1.RData")

load("Processed_Data/Cluster1.RData")

subset(Cluster1_seurat, ident = 15)
subset(TrioBrain.integrated_slim_labeled, ident = 'ort')

Cluster1[[3]]

FeaturePlot(TrioBrain.integrated_slim_labeled, features = c("Oaz", "Proc"), order=T)

## cluster1 - subclustesr markers ##

top24_Cluster1 <- Cluster1_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster1$cluster)) {
  
  top24_Cluster1_sub <- top24_Cluster1 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster1_seurat, features = top24_Cluster1_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster1_subclusters{i}.png"), width=24, height = 24)
  
}

## KD intersection for annotation

## KD (2018) annotation

Cluster1_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=celltype) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~cluster, scale='free', ncol=6)

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50_Cluster1.png"), width=25, height =20)


## Cluster 2 ##

Cluster2 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'GABA', resolution = 0.1)

DefaultAssay(Cluster2[[1]]) <- "SCT"

Cluster2_seurat <- DietSeurat(Cluster2[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))

DimPlot(Cluster2_seurat, label=T, repel=T)

Cluster2_seurat <- PrepSCTFindMarkers(Cluster2_seurat, assay = "SCT", verbose = TRUE)

Cluster2_markers <- FindAllMarkers(Cluster2_seurat, only.pos = F, min.pct = 0.15, 
                                   logfc.threshold = 0.25, test.use = "MAST")


save(Cluster2_markers, Cluster2_seurat,file="Processed_Data/Cluster2.RData")

Cluster2[[3]]

subset(Cluster2_seurat, ident = 3)
subset(TrioBrain.integrated_slim_labeled, ident = 'Dh31(I)')
subset(TrioBrain.integrated_slim_labeled, ident = 'Poxn')

subset(Cluster2_seurat, ident = 7)
subset(TrioBrain.integrated_slim_labeled, ident = 'fru')

subset(Cluster2_seurat, ident = 8) # Da6(I)
subset(Cluster2_seurat, ident = 11) # sNPF(I)

FeaturePlot(TrioBrain.integrated_slim_labeled, features = c("Oaz", "Proc"), order=T)

## Cluster2 - subclustesr markers ##

top24_Cluster2 <- Cluster2_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster2$cluster)) {
  
  top24_Cluster2_sub <- top24_Cluster2 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster2_seurat, features = top24_Cluster2_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster2_subclusters{i}.png"), width=24, height = 24)
  
}

## KD intersection for annotation

## KD (2018) annotation

Cluster2_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=celltype) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~cluster, scale='free', ncol=6)

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50_Cluster2.png"), width=25, height =20)

### Cluster 6 ###

Cluster6 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'Cluster6', resolution = 0.1, npc=10)

DefaultAssay(Cluster6[[1]]) <- "SCT"

Cluster6_seurat <- DietSeurat(Cluster6[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))

DimPlot(Cluster6_seurat, label=T, repel=T)

Cluster6_seurat <- PrepSCTFindMarkers(Cluster6_seurat, assay = "SCT", verbose = TRUE)

Cluster6_markers <- FindAllMarkers(Cluster6_seurat, only.pos = F, min.pct = 0.15, 
                                   logfc.threshold = 0.25, test.use = "MAST")


save(Cluster6_markers, Cluster6_seurat,file="Processed_Data/Cluster6.RData")
#load("Processed_Data/Cluster6.RData")

DimPlot(Cluster6_seurat)

## Cluster6 - subclustesr markers ##

top24_Cluster6 <- Cluster6_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster6$cluster)) {
  
  top24_Cluster6_sub <- top24_Cluster6 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster6_seurat, features = top24_Cluster6_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster6_subclusters{i}.png"), width=24, height = 24)
  
}

## KD intersection for annotation

Cluster6_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=celltype) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~cluster, scale='free', ncol=6)

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50_Cluster6.png"), width=25, height =20)

### Cluster 7 ###

Cluster7 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'Cluster7', resolution = 0.1)

DefaultAssay(Cluster7[[1]]) <- "SCT"

Cluster7_seurat <- DietSeurat(Cluster7[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))

DimPlot(Cluster7_seurat, label=T, repel=T)

Cluster7_seurat <- PrepSCTFindMarkers(Cluster7_seurat, assay = "SCT", verbose = TRUE)

Cluster7_markers <- FindAllMarkers(Cluster7_seurat, only.pos = F, min.pct = 0.15, 
                                   logfc.threshold = 0.25, test.use = "MAST")


save(Cluster7_markers, Cluster7_seurat,file="Processed_Data/Cluster7.RData")

Cluster7[[2]]

## Cluster7- subclustesr markers ##

top24_Cluster7 <- Cluster7_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster7$cluster)) {
  
  top24_Cluster7_sub <- top24_Cluster7 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster7_seurat, features = top24_Cluster7_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster7_subclusters{i}.png"), width=24, height = 24)
  
}

## KD intersection for annotation

## KD (2018) annotation

Cluster7_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=celltype) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~cluster, scale='free', ncol=4)

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50_Cluster7.png"), width=15, height =10)

FeaturePlot(TrioBrain.integrated_slim_labeled, feature="Mip")


### Cluster 9 (Mon) ###

Cluster9 <- Subclustering(Seurat_object=TrioBrain.integrated_slim, cluster = 9, resolution = 0.1, npc=8)

DefaultAssay(Cluster9[[1]]) <- "SCT"

Cluster9_seurat <- DietSeurat(Cluster9[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))

DimPlot(Cluster9_seurat, label=T, repel=T)

Cluster9_seurat <- PrepSCTFindMarkers(Cluster9_seurat, assay = "SCT", verbose = TRUE)

Cluster9_markers <- FindAllMarkers(Cluster9_seurat, only.pos = F, min.pct = 0.15, 
                                   logfc.threshold = 0.25, test.use = "MAST")


save(Cluster9_markers, Cluster9_seurat,file="Processed_Data/Cluster9.RData")

Cluster9[[2]]

load("Processed_Data/Cluster9.RData")

## Cluster9- subclustesr markers ##

top24_Cluster9 <- Cluster9_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster9$cluster)) {
  
  top24_Cluster9_sub <- top24_Cluster9 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster9_seurat, features = top24_Cluster9_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster9_subclusters{i}.png"), width=24, height = 24)
  
}

## KD intersection for annotation

## KD (2018) annotation

Cluster9_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=celltype) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~cluster, scale='free', ncol=6)

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50_Cluster9.png"), width=20, height =10)

FeaturePlot(Cluster9_seurat, feature = c("Tdc2", "Tbh", "Hdc"), pt.size=0.2, label=T, order=T)

FeaturePlot(Cluster9_seurat, feature = c("ple", "DAT", "Tdc2", "Tbh", "SerT", "Trh", "Hdc"), pt.size=0.2, label=T)

FeaturePlot(subset(TrioBrain.integrated_slim_labeled, idents = "MON"), feature = c("ple", "Tdc2", "Tbh", "SerT"), pt.size=0.01)

### SubCluster 9-2 (SER) ###

Cluster9_2 <- subset(Cluster9_seurat, idents = 2)

Cluster9_2 <- SCTransform(Cluster9_2)
Cluster9_2 <- RunPCA(Cluster9_2, npcs = 50, verbose = T)
Cluster9_2 <- RunUMAP(Cluster9_2, dims = 1:50)
Cluster9_2 <- FindNeighbors(Cluster9_2, reduction = "pca", dims = 1:50)
Cluster9_2 <- FindClusters(Cluster9_2, resolution = 0.8)

FeaturePlot(Cluster9_2, feature=c("Hdc", "SerT"))

DimPlot(Cluster9_2, label=T, repel=T)


### Cluster11 (fru) ###

Cluster11 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'fru', resolution = 0.1, npc=10)

DefaultAssay(Cluster11[[1]]) <- "SCT"

Cluster11_seurat <- DietSeurat(Cluster11[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))

DimPlot(Cluster11_seurat, label=T, repel=T)

Cluster11_seurat <- PrepSCTFindMarkers(Cluster11_seurat, assay = "SCT", verbose = TRUE)

Cluster11_markers <- FindAllMarkers(Cluster11_seurat, only.pos = T, min.pct = 0.15, 
                                    logfc.threshold = 0.25, test.use = "MAST")


save(Cluster11_markers, Cluster11_seurat,file="Processed_Data/Cluster11.RData")

Cluster11[[2]]
FeaturePlot(Cluster11_seurat, feature="Ms")

## Cluster11- subclustesr markers ##

top24_Cluster11 <- Cluster11_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster11$cluster)) {
  
  top24_Cluster11_sub <- top24_Cluster11 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster11_seurat, features = top24_Cluster11_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster11_subclusters{i}.png"), width=24, height = 24)
  
}

## KD intersection for annotation

## KD (2018) annotation

df_KD_cluster_select <- df_KD_cluster %>%
  dplyr::select(gene, celltype=Annotation)

Cluster11_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=celltype) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~cluster, scale='free', ncol=6)

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50_Cluster11.png"), width=20, height =10)




### Cluster12 (sNPF) ###

Cluster12 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'sNPF', resolution = 0.1, npc=10)

DefaultAssay(Cluster12[[1]]) <- "SCT"

Cluster12_seurat <- DietSeurat(Cluster12[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))

DimPlot(Cluster12_seurat, label=T, repel=T)

Cluster12_seurat <- PrepSCTFindMarkers(Cluster12_seurat, assay = "SCT", verbose = TRUE)

Cluster12_markers <- FindAllMarkers(Cluster12_seurat, only.pos = T, min.pct = 0.15, 
                                    logfc.threshold = 0.25, test.use = "MAST")


save(Cluster12_markers, Cluster12_seurat,file="Processed_Data/Cluster12.RData")

load("Processed_Data/Cluster12.RData")

Cluster12[[2]]
FeaturePlot(Cluster12_seurat, feature="CCHa2")
FeaturePlot(Cluster12_seurat, feature="AstC")
FeaturePlot(Cluster12_seurat, feature=c("sNPF", "FMRFa"))

FeaturePlot(subset(TrioBrain.integrated_slim_labeled, ident="sNPF"), feature=c("VGlut"))
FeaturePlot(subset(TrioBrain.integrated_slim_labeled, ident="sNPF"), feature=c("sNPF"))
FeaturePlot(subset(TrioBrain.integrated_slim_labeled, ident="sNPF"), feature=c("CCHa2"))
FeaturePlot(subset(TrioBrain.integrated_slim_labeled, ident="sNPF"), feature=c("Nos"))

subset(Cluster0_seurat, ident = 5)
subset(Cluster12_seurat, ident = 4)

FeaturePlot(TrioBrain.integrated_slim_labeled, feature=c("FMRFa"), reduction = 'tsne')
FeaturePlot(TrioBrain.integrated_slim_labeled, feature=c("sNPF"), reduction = 'tsne')
FeaturePlot(subset(TrioBrain.integrated_slim_labeled, ident="sNPF"), feature=c("VGlut"), reduction = 'tsne', order=T)
FeaturePlot(subset(TrioBrain.integrated_slim_labeled, ident="sNPF"), feature=c("CCHa2"), reduction = 'tsne', order=T)
FeaturePlot(subset(TrioBrain.integrated_slim_labeled, ident="sNPF"), feature=c("Nep1"), reduction = 'tsne', order=T)

VlnPlot(TrioBrain.integrated_slim_labeled, feature=c("sNPF", "FMRFa", "CCHa2"), raster=FALSE)
VlnPlot(TrioBrain.integrated_slim_labeled, feature=c("CCHa2"), raster=FALSE)

subset(Cluster12_seurat, ident = 3)
subset(TrioBrain.integrated_slim_labeled, ident = 'CCHa2')
subset(TrioBrain.integrated_slim_labeled, ident = 'Tk')


## Cluster12- subclustesr markers ##

top24_Cluster12 <- Cluster12_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster12$cluster)) {
  
  top24_Cluster12_sub <- top24_Cluster12 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster12_seurat, features = top24_Cluster12_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster12_subclusters{i}.png"), width=24, height = 24)
  
}

## KD intersection for annotation

## KD (2018) annotation

df_KD_cluster_select <- df_KD_cluster %>%
  dplyr::select(gene, celltype=Annotation)

Cluster12_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=celltype) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~cluster, scale='free', ncol=6)

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50_Cluster12.png"), width=20, height =10)

### Cluster14 (PRN) ###

Cluster14 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'PRN', resolution = 0.1, npc=10)

DefaultAssay(Cluster14[[1]]) <- "SCT"

Cluster14_seurat <- DietSeurat(Cluster14[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))

DimPlot(Cluster14_seurat, label=T, repel=T)

Cluster14_seurat <- PrepSCTFindMarkers(Cluster14_seurat, assay = "SCT", verbose = TRUE)

Cluster14_markers <- FindAllMarkers(Cluster14_seurat, only.pos = T, min.pct = 0.15, 
                                    logfc.threshold = 0.25, test.use = "MAST")


save(Cluster14_markers, Cluster14_seurat,file="Processed_Data/Cluster14.RData")

Cluster14[[2]]

## Cluster14- subclustesr markers ##

top24_Cluster14 <- Cluster14_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster14$cluster)) {
  
  top24_Cluster14_sub <- top24_Cluster14 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster14_seurat, features = top24_Cluster14_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster14_subclusters{i}.png"), width=24, height = 24)
  
}

## KD intersection for annotation

## KD (2018) annotation

df_KD_cluster_select <- df_KD_cluster %>%
  dplyr::select(gene, celltype=Annotation)

Cluster14_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=celltype) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~cluster, scale='free', ncol=6)

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50_Cluster14.png"), width=20, height =10)

FeaturePlot(Cluster14_seurat, feature=c("Tret1-1","Syt1", "Gs2"), pt.size=0.1, ncol=3, label=T)

VlnPlot(Cluster14_seurat, feature="Tret1-1")

### Cluster 15 ###

Cluster15 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'Cluster15', resolution = 0.1, npc=50)

DefaultAssay(Cluster15[[1]]) <- "SCT"

Cluster15_seurat <- DietSeurat(Cluster15[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))

DimPlot(Cluster15_seurat, label=T, repel=T)

Cluster15_seurat <- PrepSCTFindMarkers(Cluster15_seurat, assay = "SCT", verbose = TRUE)

Cluster15_markers <- FindAllMarkers(Cluster15_seurat, only.pos = T, min.pct = 0.15, 
                                   logfc.threshold = 0.25, test.use = "MAST")


save(Cluster15_markers, Cluster15_seurat,file="Processed_Data/Cluster15.RData")

Cluster15[[2]]

## Cluster15- subclustesr markers ##

top24_Cluster15 <- Cluster15_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster15$cluster)) {
  
  top24_Cluster15_sub <- top24_Cluster15 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster15_seurat, features = top24_Cluster15_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster15_subclusters{i}.png"), width=24, height = 24)
  
}

## KD intersection for annotation

## KD (2018) annotation

Cluster15_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=celltype) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~cluster, scale='free', ncol=6)

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50_Cluster15.png"), width=20, height =20)


### Cluster17 ###

Cluster17 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'PEP', resolution = 0.1, npc=10)

DefaultAssay(Cluster17[[1]]) <- "SCT"

Cluster17_seurat <- DietSeurat(Cluster17[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))

DimPlot(Cluster17_seurat, label=T, repel=T)

Cluster17_seurat <- PrepSCTFindMarkers(Cluster17_seurat, assay = "SCT", verbose = TRUE)

Cluster17_markers <- FindAllMarkers(Cluster17_seurat, only.pos = T, min.pct = 0.15, 
                                    logfc.threshold = 0.25, test.use = "MAST")


save(Cluster17_markers, Cluster17_seurat,file="Processed_Data/Cluster17.RData")

Cluster17[[2]]

## Cluster17- subclustesr markers ##

top24_Cluster17 <- Cluster17_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster17$cluster)) {
  
  top24_Cluster17_sub <- top24_Cluster17 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster17_seurat, features = top24_Cluster17_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster17_subclusters{i}.png"), width=24, height = 24)
  
}

## KD intersection for annotation

## KD (2018) annotation

df_KD_cluster_select <- df_KD_cluster %>%
  dplyr::select(gene, celltype=Annotation)

Cluster17_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=celltype) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~cluster, scale='free', ncol=6)

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50_Cluster17.png"), width=20, height =20)


### Cluster19 ###

Cluster19 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'Cluster19', resolution = 0.1, npc=10)

DefaultAssay(Cluster19[[1]]) <- "SCT"

Cluster19_seurat <- DietSeurat(Cluster19[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))

DimPlot(Cluster19_seurat, label=T, repel=T)

Cluster19_seurat <- PrepSCTFindMarkers(Cluster19_seurat, assay = "SCT", verbose = TRUE)

Cluster19_markers <- FindAllMarkers(Cluster19_seurat, only.pos = T, min.pct = 0.15, 
                                    logfc.threshold = 0.25, test.use = "MAST")


save(Cluster19_markers, Cluster19_seurat,file="Processed_Data/Cluster19.RData")

Cluster19[[2]]

load("Processed_Data/Cluster19.RData")

## Cluster19- subclustesr markers ##

top24_Cluster19 <- Cluster19_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster19$cluster)) {
  
  top24_Cluster19_sub <- top24_Cluster19 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster19_seurat, features = top24_Cluster19_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster19_subclusters{i}.png"), width=24, height = 24)
  
}

## KD intersection for annotation

## KD (2018) annotation

df_KD_cluster_select <- df_KD_cluster %>%
  dplyr::select(gene, celltype=Annotation)

Cluster19_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=celltype) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~cluster, scale='free', ncol=6)

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50_Cluster19.png"), width=20, height =20)

### Cluster20 ###

Cluster20 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'Clock', resolution = 0.1, npc=10)

DefaultAssay(Cluster20[[1]]) <- "SCT"

Cluster20_seurat <- DietSeurat(Cluster20[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))



Cluster20_seurat <- PrepSCTFindMarkers(Cluster20_seurat, assay = "SCT", verbose = TRUE)

Cluster20_markers <- FindAllMarkers(Cluster20_seurat, only.pos = T, min.pct = 0.15, 
                                    logfc.threshold = 0.25, test.use = "MAST")

#save(Cluster20_markers, Cluster20_seurat,file="Processed_Data/Cluster20.RData")

Cluster20[[2]]

DimPlot(Cluster20_seurat, label=T, repel=T)
VlnPlot(Cluster20_seurat, feature=c('Clk','tim', 'cry') )
FeaturePlot(Cluster20_seurat, feature=c('Clk'), order=T )
FeaturePlot(Cluster20_seurat, feature=c('tim'), order=T )

FeaturePlot(subset(TrioBrain.integrated_slim_labeled, ident="Clock"), feature=c("Clk", "tim"), order=T)
FeaturePlot(subset(TrioBrain.integrated_slim_labeled, ident="Clock"), feature=c("pros", "Imp"), order=T)

subset(Cluster20_seurat, ident=3)

VlnPlot(TrioBrain.integrated_slim_labeled, feature='Dh44')

#load("Processed_Data/Cluster20.RData")

## Cluster20- subclustesr markers ##

top24_Cluster20 <- Cluster20_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster20$cluster)) {
  
  top24_Cluster20_sub <- top24_Cluster20 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster20_seurat, features = top24_Cluster20_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster20_subclusters{i}.png"), width=24, height = 24)
  
}

## KD intersection for annotation

## KD (2018) annotation

df_KD_cluster_select <- df_KD_cluster %>%
  dplyr::select(gene, celltype=Annotation)

Cluster20_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=celltype) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~cluster, scale='free', ncol=6)

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50_Cluster20.png"), width=20, height =20)


### Cluster22 ###

Cluster22 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'Cluster22', resolution = 0.1, npc=10)

DefaultAssay(Cluster22[[1]]) <- "SCT"

Cluster22_seurat <- DietSeurat(Cluster22[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))


Cluster22_seurat <- PrepSCTFindMarkers(Cluster22_seurat, assay = "SCT", verbose = TRUE)

Cluster22_markers <- FindAllMarkers(Cluster22_seurat, only.pos = T, min.pct = 0.15, 
                                    logfc.threshold = 0.25, test.use = "MAST")

save(Cluster22_markers, Cluster22_seurat,file="Processed_Data/Cluster22.RData")

Cluster22[[2]]

#load("Processed_Data/Cluster22.RData")

## Cluster22- subclustesr markers ##

top24_Cluster22 <- Cluster22_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster22$cluster)) {
  
  top24_Cluster22_sub <- top24_Cluster22 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster22_seurat, features = top24_Cluster22_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster22_subclusters{i}.png"), width=24, height = 24)
  
}

### Cluster23 ###

Cluster23 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'Cluster23', resolution = 0.1, npc=10)

DefaultAssay(Cluster23[[1]]) <- "SCT"

Cluster23_seurat <- DietSeurat(Cluster23[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))


Cluster23_seurat <- PrepSCTFindMarkers(Cluster23_seurat, assay = "SCT", verbose = TRUE)

Cluster23_markers <- FindAllMarkers(Cluster23_seurat, only.pos = T, min.pct = 0.15, 
                                    logfc.threshold = 0.25, test.use = "MAST")

save(Cluster23_markers, Cluster23_seurat,file="Processed_Data/Cluster23.RData")

Cluster23[[2]]

DimPlot(Cluster23_seurat, label=T)
subset(Cluster23_seurat, ident=0)
FeaturePlot(Cluster23_seurat, feature=c('Tdc2', "Tbh"))

#load("Processed_Data/Cluster23.RData")

## Cluster23- subclustesr markers ##

top24_Cluster23 <- Cluster23_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster23$cluster)) {
  
  top24_Cluster23_sub <- top24_Cluster23 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster23_seurat, features = top24_Cluster23_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster23_subclusters{i}.png"), width=24, height = 24)
  
}





### Cluster24 (CTX) ###

Cluster24 <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'CTX', resolution = 0.1, npc=5)

DefaultAssay(Cluster24[[1]]) <- "SCT"

Cluster24_seurat <- DietSeurat(Cluster24[[1]], assays = c("RNA","SCT"), scale.data = T, dimreducs = c("umap", "tsne"))

DimPlot(Cluster24_seurat, label=T, repel=T)

Cluster24_seurat <- PrepSCTFindMarkers(Cluster24_seurat, assay = "SCT", verbose = TRUE)

Cluster24_markers <- FindAllMarkers(Cluster24_seurat, only.pos = T, min.pct = 0.15, 
                                    logfc.threshold = 0.25, test.use = "MAST")


save(Cluster24_markers, Cluster24_seurat,file="Processed_Data/Cluster24.RData")

Cluster24[[2]]

## Cluster24- subclustesr markers ##

top24_Cluster24 <- Cluster24_markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

for (i in unique(top24_Cluster24$cluster)) {
  
  top24_Cluster24_sub <- top24_Cluster24 %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(Cluster24_seurat, features = top24_Cluster24_sub$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/Subcluster/TrioBrain_markers_Cluster24_subclusters{i}.png"), width=24, height = 24)
  
}

## KD intersection for annotation

## KD (2018) annotation

df_KD_cluster_select <- df_KD_cluster %>%
  dplyr::select(gene, celltype=Annotation)

Cluster24_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_bar() +
  aes(x=celltype) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~cluster, scale='free', ncol=6)

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50_Cluster24.png"), width=20, height =20)

## check large clusters ##

### Cluster11 (fru) ###

FeaturePlot(subset(TrioBrain.integrated_slim_labeled, idents='fru'), feature = "fru")
FeaturePlot(subset(TrioBrain.integrated_slim_labeled, idents='sNPF'), feature = "sNPF")
FeaturePlot(subset(TrioBrain.integrated_slim_labeled, idents='Dh31(E)'), feature = "Dh31")
FeaturePlot(subset(TrioBrain.integrated_slim_labeled, idents='Clock'), feature = "Clk", order=T)
FeaturePlot(subset(TrioBrain.integrated_slim_labeled, idents='Clock'), feature = "tim", order=T)



### OPNs ###

OPN_markers_KD <- df_KD_cluster_markers %>%
  dplyr::left_join(., df_KD_cluster_annotation, by='cluster') %>%
  na.omit() %>%
  dplyr::filter(grepl("OPN", Annotation)) %>%
  top_n(n = 50, wt = -p.value) %>%
  top_n(n = 50, wt = mean_diff) %>%
  dplyr::select(gene, celltype=Annotation)

OPN_markers <- OPN_markers_KD$gene

Cluster9[[1]] <- AddModuleScore(object = Cluster9[[1]], features = list(OPN_markers), name = 'OPN')

VlnPlot(Cluster9[[1]], features=c("OPN1"), pt.size=0) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank()) +
  NoLegend() 

VlnPlot(Cluster9[[1]], features=c("acj6", "Oaz", "vvl", "C15", "unpg","kn", "mirr", "fru"), pt.size=0) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank()) +
  NoLegend() 

Cluster9[[2]]

DefaultAssay(Cluster9[[1]]) <- "integrated"

OPN <- Subclustering(Seurat_object=Cluster9[[1]], cluster = 2:7, resolution = 0.1)

OPN[[1]]

OPN[[1]]  <- SCTransform(OPN[[1]], method = "glmGamPoi", verbose = FALSE)

OPN_submarkers <- FindAllMarkers(OPN[[1]], only.pos = F, min.pct = 0.15, 
                                 logfc.threshold = 0.25, test.use = "MAST")

DotPlot(OPN[[1]], features=c("acj6", "Oaz", "vvl", "C15", "unpg","kn", "mirr", "fru")) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank()) +
  NoLegend() 

OPN[[2]]

## OPN subcluster DEGs ###

DotPlot(Cluster20[[1]], feature = c("tim", "brk"))

test <- Cluster20[[1]]

test@meta.data$orig.ident <- gsub("_rep1", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep2", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep3", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep4", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep5", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep6", "", test@meta.data$orig.ident)

test@meta.data$orig.ident <- factor(test@meta.data$orig.ident, 
                                                                 levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))

Idents(test) <- test$orig.ident

test_degs <- FindMarkers(test, ident.1 = "Dmel", ident.2 = "Dsec", verbose = T, test.use = "MAST") %>%
  tibble::rownames_to_column('gene') %>%
  dplyr::filter(p_val_adj < 0.05)



## KCab ##

KCab <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 3, resolution = 0.1)

KCab[[1]]  <- SCTransform(KCab[[1]], method = "glmGamPoi", verbose = FALSE)

KCab_markers <- FindAllMarkers(KCab[[1]], only.pos = F, min.pct = 0.15, 
                                                    logfc.threshold = 0.25, test.use = "MAST")

FeaturePlot(KCab[[1]], feature=c("rn", "Hr38", "sr", "SKIP"), pt.size= 0.03, order=T)

DotPlot(KCab[[1]], feature=c("crb", "Dop1R2", "rn", "Hr38"))

FeaturePlot(TrioBrain.integrated_slim_labeled, feature=c("rn", "CG4829","mtt", "Dscam3","wake","Hr38"), pt.size= 0.01, order=T)

## DEG

KCab[[2]]

DEenrichRPlot(KCab[[1]], ident.1=1, ident.2=0, enrich.database='GO_Biological_Process_AutoRIF', max.genes=200, test.use = "MAST") ## doublet
DEenrichRPlot(KCab[[1]], ident.1=2, ident.2=0, enrich.database='GO_Biological_Process_AutoRIF', max.genes=200, test.use = "MAST")
DEenrichRPlot(KCab[[1]], ident.1=3, ident.2=0, enrich.database='GO_Biological_Process_AutoRIF', max.genes=200, test.use = "MAST")



DEenrichRPlot(KCab[[1]], ident.1=2, ident.2=0, enrich.database='GO_Molecular_Function_2018', max.genes=100)

DEenrichRPlot(KCab[[1]], ident.1=2, ident.2=0, enrich.database='Coexpression_Predicted_GO_Biological_Process_2018', max.genes=200)
DEenrichRPlot(KCab[[1]], ident.1=2, ident.2=0, enrich.database='GO_Biological_Process_GeneRIF_Predicted_zscore', max.genes=100)
DEenrichRPlot(KCab[[1]], ident.1=2, ident.2=0, enrich.database='GO_Biological_Process_AutoRIF', max.genes=100)

test <- DEenrichRPlot(KCab[[1]], ident.1=1, ident.2=0, 
                      enrich.database='GO_Biological_Process_AutoRIF', max.genes=200, test.use = "MAST",
                      return.gene.list=T)
View(test$pos)

## KCr ##

KCr <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 6, resolution = 0.2)

KCr[[1]]  <- SCTransform(KCr[[1]], method = "glmGamPoi", verbose = FALSE)

KCr_markers <- FindAllMarkers(KCr[[1]], only.pos = F, min.pct = 0.15, 
                               logfc.threshold = 0.25, test.use = "MAST")

KCr[[2]]

FeaturePlot(KCr[[1]], feature=c("rn", "Hr38", "sr", "SKIP"), pt.size= 0.03, order=T)

FeaturePlot(KCr[[1]], feature=c("crb", "Hr38", "sr", "rn"), pt.size= 0.03, order=T)
DotPlot(KCr[[1]], feature=c("sNPF", "Dop1R2", "sr", "Hr38", "mamo", "Imp", "Syp"))

FeaturePlot(TrioBrain.integrated_slim_labeled, feature=c("rn", "CG14989","SKIP", "crb"), pt.size= 0.01, order=T)

## KCa'b' ##

KCapbp <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 15, resolution = 0.2)

KCapbp[[1]]  <- SCTransform(KCapbp[[1]], method = "glmGamPoi", verbose = FALSE)

KCapbp_markers <- FindAllMarkers(KCapbp[[1]], only.pos = F, min.pct = 0.15, 
                              logfc.threshold = 0.25, test.use = "MAST")

KCapbp[[2]]

FeaturePlot(KCapbp[[1]], feature=c("rn", "Hr38", "sr", "SKIP"), pt.size= 0.03, order=T)


DotPlot(KCapbp[[1]], feature=c("sNPF", "Dop1R2", "sr", "Hr38", "mamo", "Imp", "Syp"))


FeaturePlot(KCapbp[[1]], feature=c("mamo", "Hr38", "sr"), pt.size= 0.1, order=F)



Poxn_markers12 <- FindMarkers(Poxn[[1]], only.pos = F, min.pct = 0.15, ident.1 = 0, ident.2= 1,
                               logfc.threshold = 0.25, test.use = "MAST")

FeaturePlot(Poxn[[1]], feature=c("Pde6", "caup", "IA-2", "fz2", "CG31191", "trol"), pt.size= 0.03, order=T)

FeaturePlot(TrioBrain.integrated_slim_labeled, feature=c("ct", "CG14989","SKIP", "crb"), pt.size= 0.01, order=T)

DotPlot(TrioBrain.integrated_slim_labeled, feature=c("Oaz", "otp"))

#KCs <- Subclustering(Seurat_object=TrioBrain.integrated, cluster = c(7,8,16), resolution = 0.1)

  
  DefaultAssay(TrioBrain.integrated_OPN0) <- "integrated"
  
  TrioBrain.integrated_OPN0<- ScaleData(TrioBrain.integrated_OPN0)
  
  TrioBrain.integrated_OPN0 <- RunPCA(TrioBrain.integrated_OPN0, npcs = 50, verbose = FALSE)
  # t-SNE and Clustering
  TrioBrain.integrated_OPN0 <- RunUMAP(TrioBrain.integrated_OPN0, reduction = "pca", dims = 1:50)
  TrioBrain.integrated_OPN0 <- FindNeighbors(TrioBrain.integrated_OPN0, reduction = "pca", dims = 1:50)
  TrioBrain.integrated_OPN0 <- FindClusters(TrioBrain.integrated_OPN0, resolution = 0.8)
  
  TrioBrain.integrated_OPN0_species <- TrioBrain.integrated_OPN0
  
  TrioBrain.integrated_OPN0_species@meta.data$orig.ident <- gsub("_rep1", "", TrioBrain.integrated_OPN0_species@meta.data$orig.ident)
  TrioBrain.integrated_OPN0_species@meta.data$orig.ident <- gsub("_rep2", "", TrioBrain.integrated_OPN0_species@meta.data$orig.ident)
  TrioBrain.integrated_OPN0_species@meta.data$orig.ident <- gsub("_rep3", "", TrioBrain.integrated_OPN0_species@meta.data$orig.ident)
  TrioBrain.integrated_OPN0_species@meta.data$orig.ident <- gsub("_rep4", "", TrioBrain.integrated_OPN0_species@meta.data$orig.ident)
  TrioBrain.integrated_OPN0_species@meta.data$orig.ident <- gsub("_rep5", "", TrioBrain.integrated_OPN0_species@meta.data$orig.ident)
  TrioBrain.integrated_OPN0_species@meta.data$orig.ident <- gsub("_rep6", "", TrioBrain.integrated_OPN0_species@meta.data$orig.ident)
  
  TrioBrain.integrated_OPN0_species@meta.data$orig.ident <- factor(TrioBrain.integrated_OPN0_species@meta.data$orig.ident, 
                                                                   levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))
  
  p1_OPN0 <- DimPlot(TrioBrain.integrated_OPN0_species, reduction = "umap",  label = TRUE, pt.size = 0.3)
  p2_OPN0 <- DimPlot(TrioBrain.integrated_OPN0_species, reduction = "umap",  label = F, pt.size = 0.3, group.by = "orig.ident")
  
  plot_grid(p1_OPN0, p2_OPN0)
  
  ggsave("Plots/TrioBrain_Integrated_UMAP_OPN0.png", width=12, height = 6)
  
  ###
  
  TrioBrain.integrated_OPN0.list <- SplitObject(TrioBrain.integrated_OPN0_species, split.by = "orig.ident")
  
  TrioBrain.OPN0.features <- SelectIntegrationFeatures(object.list = TrioBrain.integrated_OPN0.list, nfeatures = 3000)
  TrioBrain.integrated_OPN0.list <- PrepSCTIntegration(object.list = TrioBrain.integrated_OPN0.list, assay = "SCT", 
                                                       anchor.features = TrioBrain.OPN0.features,
                                                       verbose = TRUE)
  
  TrioBrain.integrated_OPN0.list <- lapply(X = TrioBrain.integrated_OPN0.list, FUN = RunPCA, features = TrioBrain.OPN0.features)
  
  TrioBrain.integrated_OPN0.anchors<- FindIntegrationAnchors(object.list = TrioBrain.integrated_OPN0.list, normalization.method = "SCT",
                                                             anchor.features = TrioBrain.OPN0.features, verbose = T, reduction = "rpca", reference=1)
  
  TrioBrain.reintegrated_OPN0_ <- IntegrateData(anchorset = TrioBrain.integrated_OPN0.anchors, normalization.method = "SCT", verbose = TRUE,
                                                features.to.integrate = shared_genes)
  
  TrioBrain.reintegrated_OPN0_ <- RunPCA(TrioBrain.reintegrated_OPN0_, npcs = 50, verbose = T)
  
  TrioBrain.reintegrated_OPN0_ <- RunUMAP(TrioBrain.reintegrated_OPN0_, dims = 1:50)
  TrioBrain.reintegrated_OPN0_ <- RunTSNE(TrioBrain.reintegrated_OPN0_, dims = 1:50)
  TrioBrain.reintegrated_OPN0_ <- FindNeighbors(TrioBrain.reintegrated_OPN0_, reduction = "pca", dims = 1:50)
  TrioBrain.reintegrated_OPN0_ <- FindClusters(TrioBrain.reintegrated_OPN0_, resolution = 0.8)
  
  p1_OPN0 <- DimPlot(UpdateSeuratObject(TrioBrain.reintegrated_OPN0_), reduction = "umap", 
                     group.by = "orig.ident", shuffle = T, pt.size = 0.01)
  p2_OPN0 <- DimPlot(TrioBrain.reintegrated_OPN0_, reduction = "umap", label = TRUE, raster=FALSE, pt.size = 0.01)
  plot_grid(p1_OPN0, p2_OPN0)
  
  
}


### Cell type classification using module score ###

df_KD_cluster_markers <- read.csv(file="Processed_Data/KD(2018)_cluster_markers.csv")
df_KD_cluster_annotation <- read.csv(file="Processed_Data/KD(2018)_cluster_annotation.csv") %>%
  dplyr::rename(cluster=Cluster_ID)

df_KD_cluster <- df_KD_cluster_markers %>%
  dplyr::left_join(., df_KD_cluster_annotation, by='cluster') %>%
  na.omit()

df_KD_cluster_select <- df_KD_cluster %>%
  dplyr::select(gene, celltype=Annotation)

glia_markers_KD <- df_KD_cluster_markers %>%
  dplyr::left_join(., df_KD_cluster_annotation, by='cluster') %>%
  na.omit() %>%
  dplyr::filter(grepl("glia", Annotation) | grepl("Glia", Annotation) | grepl("Astrocyte", Annotation)) %>%
  dplyr::group_by(Annotation) %>%
  top_n(n = 50, wt = -p.value) %>%
  top_n(n = 50, wt = mean_diff) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene, celltype=Annotation)

glia_markers_KD_common <-  df_KD_cluster_markers %>%
  dplyr::left_join(., df_KD_cluster_annotation, by='cluster') %>%
  na.omit() %>%
  dplyr::filter(grepl("glia", Annotation) | grepl("Glia", Annotation) | grepl("Astrocyte", Annotation)) %>% 
  dplyr::group_by(Annotation) %>%
  top_n(n = 100, wt = -p.value) %>%
  top_n(n = 100, wt = mean_diff) %>%
  dplyr::ungroup() %>%
  group_by(gene) %>%
  dplyr::summarise(n_celltype = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_celltype >=4)

glia_markers <- glia_markers_KD_common$gene
ENS_markers <- dplyr::filter(glia_markers_KD, grepl("Ensheathing", celltype))$gene
AST_markers <- dplyr::filter(glia_markers_KD, grepl("Astrocyte", celltype))$gene
PRN_markers <- dplyr::filter(glia_markers_KD, grepl("Perineu", celltype))$gene
SUB_markers <- dplyr::filter(glia_markers_KD, grepl("Subperi", celltype))$gene
CTX_markers <- dplyr::filter(glia_markers_KD, grepl("Cortex", celltype))$gene
CHI_markers <- dplyr::filter(glia_markers_KD, grepl("Chiasm", celltype))$gene

TrioBrain.integrated_slim_labeled <- AddModuleScore(object = TrioBrain.integrated_slim_labeled, features = list(glia_markers), name = 'Glia')
TrioBrain.integrated_slim_labeled <- AddModuleScore(object = TrioBrain.integrated_slim_labeled, features = list(ENS_markers), name = 'ens')
TrioBrain.integrated_slim_labeled <- AddModuleScore(object = TrioBrain.integrated_slim_labeled, features = list(AST_markers), name = 'ast')
TrioBrain.integrated_slim_labeled <- AddModuleScore(object = TrioBrain.integrated_slim_labeled, features = list(PRN_markers), name = 'prn')
TrioBrain.integrated_slim_labeled <- AddModuleScore(object = TrioBrain.integrated_slim_labeled, features = list(SUB_markers), name = 'sub')
TrioBrain.integrated_slim_labeled <- AddModuleScore(object = TrioBrain.integrated_slim_labeled, features = list(CTX_markers), name = 'ctx')
TrioBrain.integrated_slim_labeled <- AddModuleScore(object = TrioBrain.integrated_slim_labeled, features = list(CHI_markers), name = 'chi')

VlnPlot(TrioBrain.integrated_slim_labeled, features=c("Glia1"), pt.size=0) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank()) +
  NoLegend() 

ggsave("Plots/TrioBrain_glia_module_score.png", width=10, height = 4)

VlnPlot(TrioBrain.integrated_slim_labeled, features=c("ens1","ast1", "prn1","sub1","ctx1","chi1"), pt.size=0) *
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank()) +
  NoLegend() 

ggsave("Plots/TrioBrain_glia_subtype_module_scores.png", width=35, height = 15)

### re-clustering ###

TrioBrain.integrated_glia <- subset(TrioBrain.integrated, idents = c(10, 13, 23, 25, 32, 33, 35, 40))

DefaultAssay(TrioBrain.integrated_glia) <- "integrated"

TrioBrain.integrated_glia<- ScaleData(TrioBrain.integrated_glia)

TrioBrain.integrated_glia <- RunPCA(TrioBrain.integrated_glia, npcs = 50, verbose = FALSE)
# t-SNE and Clustering
TrioBrain.integrated_glia <- RunUMAP(TrioBrain.integrated_glia, reduction = "pca", dims = 1:50)
TrioBrain.integrated_glia <- FindNeighbors(TrioBrain.integrated_glia, reduction = "pca", dims = 1:50)
TrioBrain.integrated_glia <- FindClusters(TrioBrain.integrated_glia, resolution = 0.8)

TrioBrain.integrated_glia_species <- TrioBrain.integrated_glia

TrioBrain.integrated_glia_species@meta.data$orig.ident <- gsub("_rep1", "", TrioBrain.integrated_glia_species@meta.data$orig.ident)
TrioBrain.integrated_glia_species@meta.data$orig.ident <- gsub("_rep2", "", TrioBrain.integrated_glia_species@meta.data$orig.ident)
TrioBrain.integrated_glia_species@meta.data$orig.ident <- gsub("_rep3", "", TrioBrain.integrated_glia_species@meta.data$orig.ident)
TrioBrain.integrated_glia_species@meta.data$orig.ident <- gsub("_rep4", "", TrioBrain.integrated_glia_species@meta.data$orig.ident)
TrioBrain.integrated_glia_species@meta.data$orig.ident <- gsub("_rep5", "", TrioBrain.integrated_glia_species@meta.data$orig.ident)
TrioBrain.integrated_glia_species@meta.data$orig.ident <- gsub("_rep6", "", TrioBrain.integrated_glia_species@meta.data$orig.ident)

TrioBrain.integrated_glia_species@meta.data$orig.ident <- factor(TrioBrain.integrated_glia_species@meta.data$orig.ident, 
                                                          levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))

p1_glia <- DimPlot(TrioBrain.integrated_glia_species, reduction = "umap",  label = TRUE, pt.size = 0.3)
p2_glia <- DimPlot(TrioBrain.integrated_glia_species, reduction = "umap",  label = F, pt.size = 0.3, group.by = "orig.ident")

plot_grid(p1_glia, p2_glia)

ggsave("Plots/TrioBrain_Integrated_UMAP_glia.png", width=17, height = 7)

DefaultAssay(TrioBrain.integrated_glia_species) <- "SCT"

TrioBrain.integrated_glia_species <- AddModuleScore(object = TrioBrain.integrated_glia_species, features = list(glia_markers), name = 'Glia')
TrioBrain.integrated_glia_species <- AddModuleScore(object = TrioBrain.integrated_glia_species, features = list(ENS_markers), name = 'ens')
TrioBrain.integrated_glia_species <- AddModuleScore(object = TrioBrain.integrated_glia_species, features = list(AST_markers), name = 'ast')
TrioBrain.integrated_glia_species <- AddModuleScore(object = TrioBrain.integrated_glia_species, features = list(PRN_markers), name = 'prn')
TrioBrain.integrated_glia_species <- AddModuleScore(object = TrioBrain.integrated_glia_species, features = list(SUB_markers), name = 'sub')
TrioBrain.integrated_glia_species <- AddModuleScore(object = TrioBrain.integrated_glia_species, features = list(CTX_markers), name = 'ctx')
TrioBrain.integrated_glia_species <- AddModuleScore(object = TrioBrain.integrated_glia_species, features = list(CHI_markers), name = 'chi')

VlnPlot(TrioBrain.integrated_glia_species, features=c("Glia1"), pt.size=0) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank()) +
  NoLegend() 

ggsave("Plots/TrioBrain_glia_module_score_subclusters.png", width=10, height = 7)

VlnPlot(TrioBrain.integrated_glia_species, features=c("ens1","ast1", "prn1","sub1","ctx1","chi1"), pt.size=0) *
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank()) +
  NoLegend() 

ggsave("Plots/TrioBrain_glia_subtype_module_scores_subclusters.png", width=20, height = 10)

## subcluster marker genes ###

DefaultAssay(TrioBrain.integrated_glia) <- "SCT"

TrioBrain.integrated_glia <- PrepSCTFindMarkers(TrioBrain.integrated_glia, assay = "SCT", verbose = TRUE)

TrioBrain.integrated_glia.markers <- FindAllMarkers(TrioBrain.integrated_glia, only.pos = TRUE, min.pct = 0.15, 
                                               logfc.threshold = 0.25, test.use = "LR", latent.vars = "orig.ident")

saveRDS(TrioBrain.integrated_glia.markers, file = "Processed_Data/TrioBrain.integrated_glia.markers.rds")

TrioBrain.integrated_glia<- ScaleData(TrioBrain.integrated_glia)

top10_Brain_glia <- TrioBrain.integrated_glia.markers %>% group_by(cluster) %>% top_n(n = 10, wt = abs(avg_log2FC))

DoHeatmap(TrioBrain.integrated_glia, features = top10_Brain_glia$gene) + NoLegend()

ggsave("Plots/TrioBrain_glia_marekr_heatmap.png", width=25, height = 18)


## Dsec vs DsecNoni

Ast_DsecNoni.markers <- FindMarkers(TrioBrain.integrated_glia, ident.1=7, ident.2=6, only.pos = F, min.pct = 0.15, logfc.threshold = 0.25)
PRN_DsecNoni.markers <- FindMarkers(TrioBrain.integrated_glia, ident.1=17, ident.2=19, only.pos = F, min.pct = 0.15, logfc.threshold = 0.25)
ENS_DsecNoni.markers <- FindMarkers(TrioBrain.integrated_glia, ident.1=4, ident.2=2, only.pos = F, min.pct = 0.15, logfc.threshold = 0.25)

PRN_DsimDsec.markers <- FindMarkers(TrioBrain.integrated_glia, ident.1=17, ident.2=14, only.pos = F, min.pct = 0.15, logfc.threshold = 0.25) %>%
  dplyr::add_rownames(var = "gene") %>%
  dplyr::mutate(pair = "Dsec_Dsim")
PRN_DmelDsec.markers <- FindMarkers(TrioBrain.integrated_glia, ident.1=17, ident.2=9, only.pos = F, min.pct = 0.15, logfc.threshold = 0.25) %>%
  dplyr::add_rownames(var = "gene") %>%
  dplyr::mutate(pair = "Dsec_Dmel")

PRN_deg <- rbind(PRN_DsimDsec.markers, PRN_DmelDsec.markers) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(hits=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-hits,gene, p_val_adj)

### CAH

FeaturePlot(TrioBrain.integrated_glia, feature = c("CAH1", "CAH2","CAH3"), order=T, ncol=3)

ggsave("Plots/TrioBrain_glia_CAH123.png", width=16, height = 5)

FeaturePlot(TrioBrain.integrated_glia_species, feature = c("CAH2","CAH3"), order=T, split.by = "orig.ident")

ggsave("Plots/TrioBrain_glia_CAH23_species.png", width=15, height = 8)

FeaturePlot(TrioBrain.integrated_slim_labeled, feature = c("CAH1", "CAH2","CAH3"), order=T, ncol=3, raster=F)

ggsave("Plots/TrioBrain_CAH123_species.png", width=16, height = 5)


### Ast

AST_DsimDsec.markers <- FindMarkers(TrioBrain.integrated_glia, ident.1=7, ident.2=8, only.pos = F, min.pct = 0.15, logfc.threshold = 0.25) %>%
  dplyr::add_rownames(var = "gene") %>%
  dplyr::mutate(pair = "Dsec_Dsim")
AST_DmelDsec.markers <- FindMarkers(TrioBrain.integrated_glia, ident.1=7, ident.2=3, only.pos = F, min.pct = 0.15, logfc.threshold = 0.25) %>%
  dplyr::add_rownames(var = "gene") %>%
  dplyr::mutate(pair = "Dsec_Dmel")

AST_deg <- rbind(AST_DsimDsec.markers, AST_DmelDsec.markers) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(hits=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-hits,gene, p_val_adj)

FeaturePlot(TrioBrain.integrated_glia_species, feature = c("Ca-alpha1T"), order=T, split.by = "orig.ident")

### Others

FeaturePlot(TrioBrain.integrated_glia_species, feature = c("pk"), order=T, split.by = "orig.ident")

FeaturePlot(TrioBrain.integrated_glia_species, feature = c("Hr38", "sr"), order=T, split.by = "orig.ident")


### OPNs ###

OPN_markers_KD <- df_KD_cluster_markers %>%
  dplyr::left_join(., df_KD_cluster_annotation, by='cluster') %>%
  na.omit() %>%
  dplyr::filter(grepl("OPN", Annotation)) %>%
  top_n(n = 50, wt = -p.value) %>%
  top_n(n = 50, wt = mean_diff) %>%
  dplyr::select(gene, celltype=Annotation)

OPN_markers <- OPN_markers_KD$gene

TrioBrain.integrated_slim_labeled <- AddModuleScore(object = TrioBrain.integrated_slim_labeled, features = list(OPN_markers), name = 'OPN')

VlnPlot(TrioBrain.integrated_slim_labeled, features=c("OPN1"), pt.size=0) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank()) +
  NoLegend() 

ggsave("Plots/TrioBrain_OPN_module_score.png", width=10, height = 4)

### re-clustering ###

TrioBrain.integrated_OPN <- subset(TrioBrain.integrated, idents = c(21,43))

DefaultAssay(TrioBrain.integrated_OPN) <- "integrated"

TrioBrain.integrated_OPN<- ScaleData(TrioBrain.integrated_OPN)

TrioBrain.integrated_OPN <- RunPCA(TrioBrain.integrated_OPN, npcs = 50, verbose = FALSE)
# t-SNE and Clustering
TrioBrain.integrated_OPN <- RunUMAP(TrioBrain.integrated_OPN, reduction = "pca", dims = 1:50)
TrioBrain.integrated_OPN <- FindNeighbors(TrioBrain.integrated_OPN, reduction = "pca", dims = 1:50)
TrioBrain.integrated_OPN <- FindClusters(TrioBrain.integrated_OPN, resolution = 0.1)

TrioBrain.integrated_OPN_species <- TrioBrain.integrated_OPN

TrioBrain.integrated_OPN_species@meta.data$orig.ident <- gsub("_rep1", "", TrioBrain.integrated_OPN_species@meta.data$orig.ident)
TrioBrain.integrated_OPN_species@meta.data$orig.ident <- gsub("_rep2", "", TrioBrain.integrated_OPN_species@meta.data$orig.ident)
TrioBrain.integrated_OPN_species@meta.data$orig.ident <- gsub("_rep3", "", TrioBrain.integrated_OPN_species@meta.data$orig.ident)
TrioBrain.integrated_OPN_species@meta.data$orig.ident <- gsub("_rep4", "", TrioBrain.integrated_OPN_species@meta.data$orig.ident)
TrioBrain.integrated_OPN_species@meta.data$orig.ident <- gsub("_rep5", "", TrioBrain.integrated_OPN_species@meta.data$orig.ident)
TrioBrain.integrated_OPN_species@meta.data$orig.ident <- gsub("_rep6", "", TrioBrain.integrated_OPN_species@meta.data$orig.ident)

TrioBrain.integrated_OPN_species@meta.data$orig.ident <- factor(TrioBrain.integrated_OPN_species@meta.data$orig.ident, 
                                                                 levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))

p1_OPN <- DimPlot(TrioBrain.integrated_OPN_species, reduction = "umap",  label = TRUE, pt.size = 0.3)
p2_OPN <- DimPlot(TrioBrain.integrated_OPN_species, reduction = "umap",  label = F, pt.size = 0.3, group.by = "orig.ident")

plot_grid(p1_OPN, p2_OPN)

ggsave("Plots/TrioBrain_Integrated_UMAP_OPN.png", width=12, height = 6)

## subcluster marker genes ###

DefaultAssay(TrioBrain.integrated_OPN) <- "SCT"

TrioBrain.integrated_OPN <- PrepSCTFindMarkers(TrioBrain.integrated_OPN, assay = "SCT", verbose = TRUE)

TrioBrain.integrated_OPN.markers <- FindAllMarkers(TrioBrain.integrated_OPN, only.pos = TRUE, min.pct = 0.15, 
                                                    logfc.threshold = 0.25, test.use = "LR", latent.vars = "orig.ident")

saveRDS(TrioBrain.integrated_OPN.markers, file = "Processed_Data/TrioBrain.integrated_OPN.markers.rds")

TrioBrain.integrated_OPN<- ScaleData(TrioBrain.integrated_OPN)

top10_Brain_OPN <- TrioBrain.integrated_OPN.markers %>% group_by(cluster) %>% top_n(n = 10, wt = abs(avg_log2FC))

DoHeatmap(TrioBrain.integrated_OPN, features = top10_Brain_OPN$gene) + NoLegend()

ggsave("Plots/TrioBrain_OPN_marekr_heatmap.png", width=12, height = 8)

TrioBrain.integrated_OPN_species <- AddModuleScore(object = TrioBrain.integrated_OPN_species, features = list(OPN_markers), name = 'OPN')

VlnPlot(TrioBrain.integrated_OPN_species, features=c("OPN1"), pt.size=0) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank()) +
  NoLegend() 

ggsave("Plots/TrioBrain_OPN_module_score_subclusters.png", width=5, height = 2.5)


### Marker genes ###

FeaturePlot(TrioBrain.integrated_OPN, feature = c("VAChT", "VGlut","Gad1"), order=T, ncol=3)

ggsave("Plots/TrioBrain_OPN_neurotransmitter.png", width=16, height = 5)

FeaturePlot(TrioBrain.integrated_OPN, feature = c("otp", "SiaT", "hth", "Oaz"), order=T)

ggsave("Plots/TrioBrain_OPN_markers.png", width=12, height = 10)

FeaturePlot(TrioBrain.integrated_slim_labeled, feature = c("otp", "SiaT", "hth", "Oaz"), order=T, ncol=2, raster=F)

ggsave("Plots/TrioBrain_markers_OPN.png", width=14, height = 13)

FeaturePlot(TrioBrain.integrated_OPN, feature = c("nAChRalpha6", "CG42458", "Eaat1"), order=T, ncol=3)

ggsave("Plots/TrioBrain_OPN_markers_subtype.png", width=16, height = 5)

FeaturePlot(TrioBrain.integrated_slim_labeled, feature = c("nAChRalpha6", "CG42458", "Eaat1"), order=T, ncol=3, raster=F)

ggsave("Plots/TrioBrain_markers.png", width=16, height = 5)

### OPN_sub0 ####

TrioBrain.integrated_OPN0 <- subset(TrioBrain.integrated_OPN, idents = 0)


DefaultAssay(TrioBrain.integrated_OPN0) <- "integrated"

TrioBrain.integrated_OPN0<- ScaleData(TrioBrain.integrated_OPN0)

TrioBrain.integrated_OPN0 <- RunPCA(TrioBrain.integrated_OPN0, npcs = 50, verbose = FALSE)
# t-SNE and Clustering
TrioBrain.integrated_OPN0 <- RunUMAP(TrioBrain.integrated_OPN0, reduction = "pca", dims = 1:50)
TrioBrain.integrated_OPN0 <- FindNeighbors(TrioBrain.integrated_OPN0, reduction = "pca", dims = 1:50)
TrioBrain.integrated_OPN0 <- FindClusters(TrioBrain.integrated_OPN0, resolution = 0.8)

TrioBrain.integrated_OPN0_species <- TrioBrain.integrated_OPN0

TrioBrain.integrated_OPN0_species@meta.data$orig.ident <- gsub("_rep1", "", TrioBrain.integrated_OPN0_species@meta.data$orig.ident)
TrioBrain.integrated_OPN0_species@meta.data$orig.ident <- gsub("_rep2", "", TrioBrain.integrated_OPN0_species@meta.data$orig.ident)
TrioBrain.integrated_OPN0_species@meta.data$orig.ident <- gsub("_rep3", "", TrioBrain.integrated_OPN0_species@meta.data$orig.ident)
TrioBrain.integrated_OPN0_species@meta.data$orig.ident <- gsub("_rep4", "", TrioBrain.integrated_OPN0_species@meta.data$orig.ident)
TrioBrain.integrated_OPN0_species@meta.data$orig.ident <- gsub("_rep5", "", TrioBrain.integrated_OPN0_species@meta.data$orig.ident)
TrioBrain.integrated_OPN0_species@meta.data$orig.ident <- gsub("_rep6", "", TrioBrain.integrated_OPN0_species@meta.data$orig.ident)

TrioBrain.integrated_OPN0_species@meta.data$orig.ident <- factor(TrioBrain.integrated_OPN0_species@meta.data$orig.ident, 
                                                                levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))

p1_OPN0 <- DimPlot(TrioBrain.integrated_OPN0_species, reduction = "umap",  label = TRUE, pt.size = 0.3)
p2_OPN0 <- DimPlot(TrioBrain.integrated_OPN0_species, reduction = "umap",  label = F, pt.size = 0.3, group.by = "orig.ident")

plot_grid(p1_OPN0, p2_OPN0)

ggsave("Plots/TrioBrain_Integrated_UMAP_OPN0.png", width=12, height = 6)

###

TrioBrain.integrated_OPN0.list <- SplitObject(TrioBrain.integrated_OPN0_species, split.by = "orig.ident")

TrioBrain.OPN0.features <- SelectIntegrationFeatures(object.list = TrioBrain.integrated_OPN0.list, nfeatures = 3000)
TrioBrain.integrated_OPN0.list <- PrepSCTIntegration(object.list = TrioBrain.integrated_OPN0.list, assay = "SCT", 
                                                     anchor.features = TrioBrain.OPN0.features,
                                     verbose = TRUE)

TrioBrain.integrated_OPN0.list <- lapply(X = TrioBrain.integrated_OPN0.list, FUN = RunPCA, features = TrioBrain.OPN0.features)

TrioBrain.integrated_OPN0.anchors<- FindIntegrationAnchors(object.list = TrioBrain.integrated_OPN0.list, normalization.method = "SCT",
                                            anchor.features = TrioBrain.OPN0.features, verbose = T, reduction = "rpca", reference=1)

TrioBrain.reintegrated_OPN0_ <- IntegrateData(anchorset = TrioBrain.integrated_OPN0.anchors, normalization.method = "SCT", verbose = TRUE,
                                      features.to.integrate = shared_genes)

TrioBrain.reintegrated_OPN0_ <- RunPCA(TrioBrain.reintegrated_OPN0_, npcs = 50, verbose = T)

TrioBrain.reintegrated_OPN0_ <- RunUMAP(TrioBrain.reintegrated_OPN0_, dims = 1:50)
TrioBrain.reintegrated_OPN0_ <- RunTSNE(TrioBrain.reintegrated_OPN0_, dims = 1:50)
TrioBrain.reintegrated_OPN0_ <- FindNeighbors(TrioBrain.reintegrated_OPN0_, reduction = "pca", dims = 1:50)
TrioBrain.reintegrated_OPN0_ <- FindClusters(TrioBrain.reintegrated_OPN0_, resolution = 0.8)

p1_OPN0 <- DimPlot(UpdateSeuratObject(TrioBrain.reintegrated_OPN0_), reduction = "umap", 
                   group.by = "orig.ident", shuffle = T, pt.size = 0.01)
p2_OPN0 <- DimPlot(TrioBrain.reintegrated_OPN0_, reduction = "umap", label = TRUE, raster=FALSE, pt.size = 0.01)
plot_grid(p1_OPN0, p2_OPN0)

### DEG

TrioBrain.reintegrated_OPN0_SCT <- TrioBrain.reintegrated_OPN0_

DefaultAssay(TrioBrain.reintegrated_OPN0_SCT) <- "RNA"

TrioBrain.reintegrated_OPN0_SCT <- SCTransform(TrioBrain.reintegrated_OPN0_SCT, method = "glmGamPoi", verbose = FALSE)

TrioBrain.reintegrated_OPN0_SCT.markers <- FindAllMarkers(TrioBrain.reintegrated_OPN0_SCT, only.pos = TRUE, min.pct = 0.15, 
                                                   logfc.threshold = 0.25)

#saveRDS(TrioBrain.reintegrated_OPN0_.markers, file = "Processed_Data/TrioBrain.reintegrated_OPN0_.markers.rds")

top10_Brain_OPN0_SCT <- TrioBrain.reintegrated_OPN0_SCT.markers %>% group_by(cluster) %>% top_n(n = 10, wt = abs(avg_log2FC))

DoHeatmap(TrioBrain.reintegrated_OPN0_SCT, features = top10_Brain_OPN0_SCT$gene, assay="SCT") + NoLegend()

ggsave("Plots/TrioBrain_OPN_marekr_heatmap.png", width=12, height = 8)



### Let's write a function

DimPlot(subset(TrioBrain.integrated, idents=32), pt.size=0.01)


