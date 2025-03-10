library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(sctransform)
library(cowplot)
library(broom)
library(MAST)
library(enrichR)
library(pheatmap)
library(RColorBrewer)
setEnrichrSite("FlyEnrichr")

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

Dmel_decon_combined_DF <- readRDS(file = "Processed_Data/Dmel_combined_full_intron_decon_DF.rds")
Dsim_to_DmelRef_decon_combined_DF <- readRDS(file = "Processed_Data/Dsim_to_DmelRef_combined_full_intron_decon_DF.rds")
Dsec_to_DmelRef_decon_combined_DF <- readRDS(file = "Processed_Data/Dsec_to_DmelRef_combined_full_intron_decon_DF.rds")
DsecNoni_to_DmelRef_decon_combined_DF <- readRDS(file = "Processed_Data/DsecNoni_to_DmelRef_combined_full_intron_decon_DF.rds")

DefaultAssay(Dmel_decon_combined_DF) <- "SCT"
DefaultAssay(Dsim_to_DmelRef_decon_combined_DF) <- "SCT"
DefaultAssay(Dsec_to_DmelRef_decon_combined_DF) <- "SCT"
DefaultAssay(DsecNoni_to_DmelRef_decon_combined_DF) <- "SCT"

all.genes_Dmel <- rownames(x = Dmel_decon_combined_DF)
all.genes_Dsim_to_DmelRef <- rownames(x = Dsim_to_DmelRef_decon_combined_DF)
all.genes_Dsec_to_DmelRef <- rownames(x = Dsec_to_DmelRef_decon_combined_DF)
all.genes_DsecNoni_to_DmelRef <- rownames(x = DsecNoni_to_DmelRef_decon_combined_DF)

shared_genes_DmelRef <- all.genes_Dmel[(all.genes_Dmel %in% all.genes_Dsim_to_DmelRef) & (all.genes_Dmel %in% all.genes_Dsec_to_DmelRef) & (all.genes_Dmel %in% all.genes_DsecNoni_to_DmelRef)]
save(shared_genes_DmelRef,file="Processed_Data/shared_genes_DmelRef.RData")


DmelBrain_ortho_only_DF <- subset(Dmel_decon_combined_DF, features = shared_genes_DmelRef)
Dsim_to_DmelRef_Brain_ortho_only_DF <- subset(Dsim_to_DmelRef_decon_combined_DF, features = shared_genes_DmelRef)
Dsec_to_DmelRef_Brain_ortho_only_DF <- subset(Dsec_to_DmelRef_decon_combined_DF, features = shared_genes_DmelRef)
DsecNoni_to_DmelRef_Brain_ortho_only_DF <- subset(DsecNoni_to_DmelRef_decon_combined_DF, features = shared_genes_DmelRef)

DefaultAssay(DmelBrain_ortho_only_DF) <- "SCT"
DefaultAssay(Dsim_to_DmelRef_Brain_ortho_only_DF) <- "SCT"
DefaultAssay(Dsec_to_DmelRef_Brain_ortho_only_DF) <- "SCT"
DefaultAssay(DsecNoni_to_DmelRef_Brain_ortho_only_DF) <- "SCT"

TrioBrain_DmelRef_DF.list <- list(DmelBrain_ortho_only_DF, Dsim_to_DmelRef_Brain_ortho_only_DF, Dsec_to_DmelRef_Brain_ortho_only_DF, DsecNoni_to_DmelRef_Brain_ortho_only_DF)
TrioBrain_DmelRef_DF.list <- lapply(X = TrioBrain_DmelRef_DF.list, FUN = DietSeurat, scale.data = T)

#save(TrioBrain_DmelRef_DF.list, file="Processed_Data/TrioBrain_DmelRef_DF.list.RData")

### run on cluster ###

load("Processed_Data/TrioBrain_DmelRef_DF.list.RData")
load("Processed_Data/shared_genes_DmelRef.RData")

TrioBrain_DmelRef.features <- SelectIntegrationFeatures(object.list = TrioBrain_DmelRef_DF.list, nfeatures = 3000)
TrioBrain_DmelRef_DF.list <- PrepSCTIntegration(object.list = TrioBrain_DmelRef_DF.list, assay = "SCT", anchor.features = TrioBrain_DmelRef.features,
                                     verbose = TRUE)
TrioBrain_DmelRef_DF.list <- lapply(X = TrioBrain_DmelRef_DF.list, FUN = RunPCA, features = TrioBrain_DmelRef.features)

TrioBrain_DmelRef.anchors <- FindIntegrationAnchors(object.list = TrioBrain_DmelRef_DF.list, normalization.method = "SCT",
                                            anchor.features = TrioBrain_DmelRef.features, verbose = T, reduction = "rpca", reference=1)

TrioBrain_DmelRef.integrated <- IntegrateData(anchorset = TrioBrain_DmelRef.anchors, normalization.method = "SCT", verbose = TRUE,
                                      features.to.integrate = shared_genes_DmelRef)

TrioBrain_DmelRef.integrated <- RunPCA(TrioBrain_DmelRef.integrated, npcs = 50, verbose = T)

TrioBrain_DmelRef.integrated <- DietSeurat(TrioBrain_DmelRef.integrated, scale.data=F, dimreducs = "pca")

TrioBrain_DmelRef.integrated <- RunUMAP(TrioBrain_DmelRef.integrated, dims = 1:50)
TrioBrain_DmelRef.integrated <- RunTSNE(TrioBrain_DmelRef.integrated, dims = 1:50)
TrioBrain_DmelRef.integrated <- FindNeighbors(TrioBrain_DmelRef.integrated, reduction = "pca", dims = 1:50)
TrioBrain_DmelRef.integrated <- FindClusters(TrioBrain_DmelRef.integrated, resolution = 0.2)

DimPlot(TrioBrain_DmelRef.integrated, label=T)

saveRDS(TrioBrain_DmelRef.integrated, file = "Processed_Data/TrioBrain_DmelRef_combined_rpca.rds")

####

#TrioBrain_DmelRef.integrated <- readRDS(file = "Processed_Data/TrioBrain_DmelRef_combined_rpca.rds")

DefaultAssay(TrioBrain_DmelRef.integrated) <- "SCT"
FeaturePlot(TrioBrain_DmelRef.integrated, feature=c("Tret1-1"), raster=FALSE,order=T)

FeaturePlot(TrioBrain_DmelRef.integrated, feature=c("Gad1", "ChAT"), raster=FALSE,order=T)

### DEG analysis ###

DefaultAssay(TrioBrain_DmelRef.integrated) <- "RNA"

TrioBrain_DmelRef.integrated_slim <- DietSeurat(TrioBrain_DmelRef.integrated, assays = c("RNA"), dimreducs = c("umap", "tsne"))

TrioBrain_DmelRef.integrated_slim <- SCTransform(TrioBrain_DmelRef.integrated_slim, method = "glmGamPoi", verbose = FALSE)

saveRDS(TrioBrain_DmelRef.integrated_slim, file = "Processed_Data/TrioBrain_DmelRef_DF.integrated_slim.rds")
#TrioBrain_DmelRef.integrated_slim <- readRDS(file = "Processed_Data/TrioBrain_DmelRef_DF.integrated_slim.rds")

DimPlot(TrioBrain_DmelRef.integrated_slim, reduction = "tsne", label = TRUE, raster=FALSE, pt.size = 0.01)

p1 <- DimPlot(TrioBrain_DmelRef.integrated_slim, reduction = "tsne", label = TRUE, raster=FALSE, pt.size = 0.01) + NoLegend() + NoAxes()
p2 <- DimPlot(TrioBrain_DmelRef.integrated_slim, reduction = "umap", label = TRUE, raster=FALSE, pt.size = 0.01) + NoLegend() + NoAxes()
plot_grid(p1, p2)

ggsave("Plots/TrioBrain_DmelRef_DF_Integrated_tsne.png", width=17, height = 8)

FeaturePlot(TrioBrain_DmelRef.integrated_slim, feature=c("pros", "Imp", 
                                                 "ChAT", "VGlut", "Gad1"), raster=F, reduction = 'tsne', pt.size=0.01)

### plots ###

TrioBrain_DmelRef.integrated_slim_species <- TrioBrain_DmelRef.integrated_slim

TrioBrain_DmelRef.integrated_slim_species@meta.data$orig.ident <- gsub("_rep1", "", TrioBrain_DmelRef.integrated_slim_species@meta.data$orig.ident)
TrioBrain_DmelRef.integrated_slim_species@meta.data$orig.ident <- gsub("_rep2", "", TrioBrain_DmelRef.integrated_slim_species@meta.data$orig.ident)
TrioBrain_DmelRef.integrated_slim_species@meta.data$orig.ident <- gsub("_rep3", "", TrioBrain_DmelRef.integrated_slim_species@meta.data$orig.ident)
TrioBrain_DmelRef.integrated_slim_species@meta.data$orig.ident <- gsub("_rep4", "", TrioBrain_DmelRef.integrated_slim_species@meta.data$orig.ident)
TrioBrain_DmelRef.integrated_slim_species@meta.data$orig.ident <- gsub("_rep5", "", TrioBrain_DmelRef.integrated_slim_species@meta.data$orig.ident)
TrioBrain_DmelRef.integrated_slim_species@meta.data$orig.ident <- gsub("_rep6", "", TrioBrain_DmelRef.integrated_slim_species@meta.data$orig.ident)

TrioBrain_DmelRef.integrated_slim_species@meta.data$orig.ident <- factor(TrioBrain_DmelRef.integrated_slim_species@meta.data$orig.ident, 
                                                                 levels=c("Dmel", "Dsim_to_DmelRef", "Dsec", "DsecNoni_to_DmelRef"))

plot_dmel_tsne <- DimPlot(subset(TrioBrain_DmelRef.integrated_slim_species, orig.ident == "Dmel"), 
                          reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01,
                          cols=c("#E41A1C")) + NoLegend() + labs(title=NULL) + NoAxes()

plot_Dsim_to_DmelRef_tsne <- DimPlot(subset(TrioBrain_DmelRef.integrated_slim_species, orig.ident == "Dsim_to_DmelRef"), 
                          reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01,
                          cols=c("#4DAF4A")) + NoLegend() + labs(title=NULL) + NoAxes()

plot_Dsec_to_DmelRef_tsne <- DimPlot(subset(TrioBrain_DmelRef.integrated_slim_species, orig.ident == "Dsec"), 
                          reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01,
                          cols=c("#377EB8")) + NoLegend() + labs(title=NULL) + NoAxes()

plot_trio_merged_tsne <- DimPlot(subset(TrioBrain_DmelRef.integrated_slim_species, orig.ident != "DsecNoni_to_DmelRef"), 
                                 reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01,
                                 cols=c("#E41A1C", "#4DAF4A",  "#377EB8")) + NoLegend() + labs(title=NULL) + NoAxes()

plot_tsne_species <- cowplot::plot_grid(plot_dmel_tsne, plot_Dsim_to_DmelRef_tsne,
                                        plot_Dsec_to_DmelRef_tsne, plot_trio_merged_tsne, ncol=2)

plot_tsne_species

ggsave("Plots/Manuscript/Fig1_tsne_species.png", width=13, height = 12)

#DimPlot(TrioBrain_DmelRef.integrated_slim_species, reduction = "umap", shuffle = T, raster=FALSE, label=T, pt.size = 0.01) + NoLegend() + labs(title=NULL)
#ggsave("Plots/TrioBrain_DmelRef_DF_Integrated_umap_label.png", width=11, height = 10)

p1 <- DimPlot(TrioBrain_DmelRef.integrated_slim_species, reduction = "umap", 
              group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01, cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL)
p2 <- DimPlot(TrioBrain_DmelRef.integrated_slim_species, reduction = "umap", label = F, raster=F, pt.size = 0.01) + NoLegend()
plot_grid(p1, p2)

#ggsave("Plots/TrioBrain_DmelRef_DF_Integrated_umap.png", width=23, height = 10)

### marker expression ###

FeaturePlot(TrioBrain_DmelRef.integrated_slim, feature = 'repo', order=T, reduction = "tsne")
FeaturePlot(TrioBrain_DmelRef.integrated_slim, feature = 'SK', order=T, reduction = "tsne")

Glia <- WhichCells(TrioBrain_DmelRef.integrated_slim, expression = repo>0)
MAcells <- WhichCells(TrioBrain_DmelRef.integrated_slim, expression = Vmat>3)
DOPcells <- WhichCells(TrioBrain_DmelRef.integrated_slim, expression = ple>1)
AchCells <- WhichCells(TrioBrain_DmelRef.integrated_slim, expression = VAChT>1 & VGlut<=1.5 & Gad1<=1.2)
GluCells <- WhichCells(TrioBrain_DmelRef.integrated_slim, expression = VGlut>1.5 & Gad1<=1.2 & VAChT<=1)
GabaCells <- WhichCells(TrioBrain_DmelRef.integrated_slim, expression = Gad1>1.2 & VGlut<=1.5 & VAChT<=1)
MultiCells <- WhichCells(TrioBrain_DmelRef.integrated_slim, expression = Gad1>1.2 & VAChT>1 | VAChT>1 &VGlut>1.5 | VGlut>1.5 &Gad1>1.2)

p1_type <- DimPlot(TrioBrain_DmelRef.integrated_slim, label=F,shuffle=T,
                   cells.highlight= list(VAChT=AchCells, VGlut=GluCells, Gad1=GabaCells, Repo=Glia, Vmat=MAcells), 
                   reduction = "tsne", 
                   pt.size = 0.01,sizes.highlight = 0.01,
                   cols.highlight = c("darkorange", "darkgreen", "darkblue","black", "darkred"), cols= "grey40", raster=FALSE) + NoLegend() + NoAxes()

p2_type <- DimPlot(TrioBrain_DmelRef.integrated_slim, label=F,shuffle=T,
                   cells.highlight= list(VAChT=AchCells, VGlut=GluCells, Gad1=GabaCells, Repo=Glia, Vmat=MAcells), 
                   reduction = "umap", 
                   pt.size = 0.01,sizes.highlight = 0.01,
                   cols.highlight = c("darkorange", "darkgreen", "darkblue","black", "darkred"), cols= "grey40", raster=FALSE) + NoLegend() + NoAxes()

plot_grid(p1_type, p2_type)
ggsave("Plots/TrioBrain_DmelRef_DF_Integrated_tsne_UMAP_nt.png", width=17, height = 8)

DimPlot(TrioBrain_DmelRef.integrated_slim, label=F,shuffle=T,
        cells.highlight= list(VAChT=AchCells, VGlut=GluCells, Gad1=GabaCells, Repo=Glia, Vmat=MAcells), 
        reduction = "umap", 
        pt.size = 0.01,sizes.highlight = 0.01,
        cols.highlight = c("darkorange", "darkgreen", "darkblue","black", "darkred"), cols= "grey40", raster=T) + NoAxes()

ggsave("Plots/TrioBrain_DmelRef_DF_Integrated_umap_nt_for_label.pdf", width=9, height = 7)


p3 <- DimPlot(TrioBrain_DmelRef.integrated_slim, label=F,shuffle=T,
              cells.highlight= list(VAChT=AchCells, VGlut=GluCells, Gad1=GabaCells, Repo=Glia, Vmat=MAcells), 
              reduction = "umap", 
              pt.size = 0.01,sizes.highlight = 0.01, 
              cols.highlight = c("darkorange", "darkgreen", "darkblue","black", "darkred"), cols= "grey40", raster=F) +NoLegend()

plot_grid(p1, p3, p2, nrow=1)

ggsave("Plots/Manuscript/Fig1a_raster.pdf", width=30, height = 10)

p1_lab <- DimPlot(TrioBrain_DmelRef.integrated_slim_species, reduction = "umap", 
                  group.by = "orig.ident", shuffle = T, raster=T, pt.size = 0.01, cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + labs(title=NULL)

p2_lab <- DimPlot(TrioBrain_DmelRef.integrated_slim_species, reduction = "umap", label = F, raster=T, pt.size = 0.01)

p3_lab <- DimPlot(TrioBrain_DmelRef.integrated_slim, label=F,shuffle=T,
                  cells.highlight= list(VAChT=AchCells, VGlut=GluCells, Gad1=GabaCells, Repo=Glia, Vmat=MAcells), 
                  reduction = "umap", 
                  pt.size = 0.01,sizes.highlight = 0.01, 
                  cols.highlight = c("darkorange", "darkgreen", "darkblue","black", "darkred"), cols= "grey40", raster=T)

plot_grid(p1_lab, p3_lab, p2_lab, nrow=1)

ggsave("Plots/Manuscript/Fig1a_raster_lab.pdf", width=30, height = 10)

## pros, Imp (master regulator) ##

ImpCells <- WhichCells(TrioBrain_DmelRef.integrated_slim, expression = Imp>1.2)
Proscells <- WhichCells(TrioBrain_DmelRef.integrated_slim, expression = pros>2)

p1_mr <- DimPlot(TrioBrain_DmelRef.integrated_slim, label=F,shuffle=T,
                 cells.highlight= list(Pros=Proscells, Imp=ImpCells), 
                 reduction = "tsne", 
                 pt.size = 0.01,sizes.highlight = 0.01,
                 cols.highlight = c("darkcyan", "darkmagenta"), cols= "grey40", raster=FALSE) + NoLegend() + NoAxes()

p2_mr <- DimPlot(TrioBrain_DmelRef.integrated_slim, label=F,shuffle=T,
                 cells.highlight= list(Pros=Proscells, Imp=ImpCells), 
                 reduction = "umap", 
                 pt.size = 0.01,sizes.highlight = 0.01,
                 cols.highlight = c("darkcyan", "darkmagenta"), cols= "grey40", raster=FALSE) + NoLegend() + NoAxes()

plot_grid(p1_mr, p2_mr)
ggsave("Plots/TrioBrain_DmelRef_DF_Integrated_tsne_UMAP_mr.png", width=17, height = 8)


DimPlot(TrioBrain_DmelRef.integrated_slim, label=F,shuffle=T,
        cells.highlight= list(Pros=Proscells, Imp=ImpCells), 
        reduction = "umap", 
        pt.size = 0.01,sizes.highlight = 0.01,
        cols.highlight = c("darkcyan", "darkmagenta"), cols= "grey40", raster=T) + NoAxes()

ggsave("Plots/TrioBrain_DmelRef_DF_Integrated_umap_mr_for_label.pdf", width=9, height = 7)

### Iterative subclustering - see Interative_subclustering.r ###

TrioBrain_DmelRef.integrated_slim_ISub_DF <- readRDS(file = "Processed_Data/TrioBrain_DmelRef.integrated_slim_ISub_DF.rds")

ClusterIdents <- unique(Idents(TrioBrain_DmelRef.integrated_slim_ISub_DF))

DimPlot(TrioBrain_DmelRef.integrated_slim_ISub_DF, reduction='umap', label=T) + NoLegend()

## QC and additional doublet filtering ##

clusters_15 <- ClusterIdents[grep(15, ClusterIdents)] ## PRN
VlnPlot(subset(TrioBrain_DmelRef.integrated_slim_ISub_DF, idents = clusters_15), features = c("nFeature_RNA")) + NoLegend()

clusters_7 <- as.character(ClusterIdents[grep('cluster7', ClusterIdents)]) ## ENS
VlnPlot(subset(TrioBrain_DmelRef.integrated_slim_ISub_DF, idents = c(7,clusters_7)), features = c("nFeature_RNA")) + NoLegend()

clusters_0 <- as.character(ClusterIdents[grep('cluster0', ClusterIdents)])
VlnPlot(subset(TrioBrain_DmelRef.integrated_slim_ISub_DF, idents = c(0,clusters_0)), features = c("nFeature_RNA")) + NoLegend()

## nFeature_RNA

cluster_features_nFeature_RNA <- TrioBrain_DmelRef.integrated_slim_ISub_DF@meta.data %>%
  group_by(seurat_clusters) %>%
  mutate(total_cell_count=n(), 
         total_mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
         total_median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
         total_sd_nFeature_RNA = sd(nFeature_RNA, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  group_by(seurat_clusters, CellType) %>%
  mutate(
    mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    sd_nFeature_RNA = sd(nFeature_RNA, na.rm = TRUE),
    cell_count = n(),
    fraction_within_cluter = n()/total_cell_count
  ) %>%
  dplyr::distinct(seurat_clusters, CellType, total_cell_count, total_mean_nFeature_RNA, total_median_nFeature_RNA, total_sd_nFeature_RNA,
                  mean_nFeature_RNA, median_nFeature_RNA, sd_nFeature_RNA, cell_count, fraction_within_cluter) %>%
  dplyr::mutate(upper = total_mean_nFeature_RNA+0.5*total_sd_nFeature_RNA,
                lower = total_mean_nFeature_RNA-0.5*total_sd_nFeature_RNA)


cluster_features_nFeature_RNA %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::filter(n()>=2) %>%
  dplyr::ungroup() %>%
  ggplot(.) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(y=total_mean_nFeature_RNA), linetype=2, alpha=0.5, color='red') +
  geom_point() +
  geom_line(alpha=0.5, linewidth=0.3) +
  aes(x=fraction_within_cluter, y=mean_nFeature_RNA) +
  theme_bw() +
  facet_wrap(~seurat_clusters, ncol=5, scales='free')

## nCount_RNA

cluster_features_nCount_RNA <- TrioBrain_DmelRef.integrated_slim_ISub_DF@meta.data %>%
  group_by(seurat_clusters) %>%
  mutate(total_cell_count=n(), 
         total_mean_nCount_RNA = mean(nCount_RNA, na.rm = TRUE),
         total_median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
         total_sd_nCount_RNA = sd(nCount_RNA, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  group_by(seurat_clusters, CellType) %>%
  mutate(
    mean_nCount_RNA = mean(nCount_RNA, na.rm = TRUE),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    sd_nCount_RNA = sd(nCount_RNA, na.rm = TRUE),
    cell_count = n(),
    fraction_within_cluter = n()/total_cell_count
  ) %>%
  dplyr::distinct(seurat_clusters, CellType, total_cell_count, total_mean_nCount_RNA, total_median_nCount_RNA, total_sd_nCount_RNA,
                  mean_nCount_RNA, median_nCount_RNA, sd_nCount_RNA, cell_count, fraction_within_cluter) %>%
  dplyr::mutate(upper = total_mean_nCount_RNA+total_sd_nCount_RNA,
                lower = total_mean_nCount_RNA-total_sd_nCount_RNA)

cluster_features_nCount_RNA %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::filter(n()>=2) %>%
  dplyr::ungroup() %>%
  ggplot(.) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(y=total_mean_nCount_RNA), linetype=2, alpha=0.5, color='red') +
  geom_point() +
  geom_line(alpha=0.5, linewidth=0.3) +
  aes(x=fraction_within_cluter, y=mean_nCount_RNA) +
  theme_bw() +
  facet_wrap(~seurat_clusters, ncol=5, scales='free')




## Annotation ##

### get markers ##

library(future)
plan(multisession, workers = 16)

#TrioBrain_DmelRef_ISub_DF.integrated.markers <- FindAllMarkers(TrioBrain_DmelRef.integrated_slim_ISub_DF, only.pos = TRUE, min.pct = 0.15, 
#                                               logfc.threshold = 0.25, test.use = "MAST")

plan(sequential)  # Reset to default plan
gc()  # Cleanup unused memory

#saveRDS(TrioBrain_DmelRef_ISub_DF.integrated.markers, file = "Processed_Data/TrioBrain_DmelRef_ISub_DF.integrated.markers.rds")
TrioBrain_DmelRef_ISub_DF.integrated.markers <- readRDS(file = "Processed_Data/TrioBrain_DmelRef_ISub_DF.integrated.markers.rds")

### annotate based on markeres ###

TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled <- TrioBrain_DmelRef.integrated_slim_ISub_DF

total_cell_number <- nrow(TrioBrain_DmelRef.integrated_slim_ISub_DF@meta.data)

ISub_DF_metadata_summary <- TrioBrain_DmelRef.integrated_slim_ISub_DF@meta.data %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-n) %>%
  dplyr::mutate(percent = n/total_cell_number*100)

anno_complete <- NULL

annotate_celltype <- function(seurat_object=TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled,
                              cell_type, df_markers, exclude='NA', 
                              metadata=ISub_DF_metadata_summary) {
  
  anno_cluster <- NULL
  
  celltype_clusters <- df_markers %>%
    dplyr::filter(celltype==cell_type & !cluster %in% exclude)
  
  cluster_ids <- as.character(celltype_clusters$cluster)
  
  if(length(cluster_ids)>1) {
    
    cluster_size <- metadata %>%
      dplyr::filter(CellType %in% cluster_ids) %>%
      dplyr::arrange(-percent) %>%
      tibble::rownames_to_column() %>%
      dplyr::mutate(cluster_name=paste(cell_type,rowname,sep="_"))
    
    for (i in 1:nrow(cluster_size)) {
      
      Idents(seurat_object, 
             cells = WhichCells(seurat_object, idents=cluster_size$CellType[i])) <- cluster_size$cluster_name[i]
      
      anno_cluster <- c(anno_cluster, as.character(cluster_size$CellType[i]))
      
    }
    
  }
  
  else if(length(cluster_ids)==1) {
    
    Idents(seurat_object, 
           cells = WhichCells(seurat_object, idents=cluster_ids)) <- cell_type
    
    anno_cluster <- cluster_ids
    
  }
  
  return_list <- list(seurat_object, anno_cluster)
  
  return(return_list)
  
}


## annotation workflow - define non-neuronal clusters -> define KC clusters -> define neuronal clusters by neurotransmitter markers 

Markers_sigs <- TrioBrain_DmelRef_ISub_DF.integrated.markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::select(gene, cluster,p_val_adj) %>%
  tidyr::spread(key=gene, value=p_val_adj) %>%
  tidyr::gather(-cluster, key='gene', value='p_val_adj')

## glial population ## (non-neuronal, no blood, no fat from new pipeline)

glial <- c("alrm", "Eaat1", "Gat", "axo","Gs2","Tret1-1","dve", "sog", "baz","CG40470")

df_glial_markers <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene %in% c(glial)) %>%
  dplyr::mutate(gene = ifelse(gene == "Tret1-1", "Tret1_1", gene)) %>%
  tidyr::spread(gene, p_val_adj) %>%
  dplyr::mutate(celltype=ifelse(!is.na(axo) & !is.na(Gs2) & is.na(sog), "ENS",
                                       ifelse(is.na(axo) & !is.na(Gs2) & !is.na(alrm), "AST",
                                              ifelse(!is.na(Tret1_1) & !is.na(dve) & !is.na(sog), "PRN",
                                                     ifelse(!is.na(Gs2) & !is.na(baz), "SUB",
                                                            ifelse(!is.na(axo) & !is.na(CG40470), "CTX", "others")))))) %>%
  dplyr::arrange(celltype)

clu14_markers <- TrioBrain_DmelRef_ISub_DF.integrated.markers %>%
  dplyr::filter(cluster ==14)

## annotate glial clusters ##

glial_celltypes <- c("ENS","AST","PRN","SUB", "CTX")

anno_complete <- NULL

for (ct in glial_celltypes) {
  
  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_glial_markers)
  
  TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
  anno_complete <- c(anno_complete, out_annotate_celltype[[2]])
}

#DimPlot(subset(TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled, idents = c("ENS_1","ENS_2","ENS_3", "PRN", "SUB", "AST", "CTX")))

## KC populations ##

KC_markers <- c("crb","CG32204","Pka-C1","mamo","sNPF")

df_KC_markers <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene %in% c(KC_markers)) %>%
  dplyr::mutate(gene = ifelse(gene == "Pka-C1", "Pka_C1", gene)) %>%
  tidyr::spread(gene, p_val_adj) %>%
  dplyr::mutate(celltype=ifelse(!is.na(crb) & !is.na(CG32204) & !is.na(Pka_C1) & !is.na(sNPF) & is.na(mamo), "αβ-KC",
                                ifelse(!is.na(mamo) & !is.na(CG32204) & !is.na(Pka_C1) & !is.na(sNPF) & is.na(crb), "γ-KC",
                                       ifelse(!is.na(mamo) & !is.na(CG32204) & !is.na(Pka_C1) & is.na(sNPF) & is.na(crb), "α'β'-KC", "others")))) %>%
  dplyr::arrange(celltype)

## annotate KC clusters ##

KC_celltypes <- c("αβ-KC","γ-KC","α'β'-KC")

for (ct in KC_celltypes) {
  
  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_KC_markers, exclude=anno_complete)
  
  TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
  anno_complete <- c(anno_complete, out_annotate_celltype[[2]])
  
}

DimPlot(subset(TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled, idents = c("αβ-KC_1","αβ-KC_2","αβ-KC_3","γ-KC_1","γ-KC_2","γ-KC_3","γ-KC_4","α'β'-KC")))

FeaturePlot(subset(TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled, idents = c("αβ-KC_1","αβ-KC_2","αβ-KC_3","γ-KC_1","γ-KC_2","γ-KC_3","γ-KC_4","α'β'-KC")),
            feature='CG33098')

## Monoaminergic neurons ##

MA_markers <- c("Vmat", "ple", "DAT","Trhn","SerT","Tdc2","Tbh")

df_MA_markers <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene %in% c(MA_markers)) %>%
  tidyr::spread(gene, p_val_adj) %>%
  dplyr::filter(!is.na(Vmat)) %>%
  dplyr::mutate(celltype=ifelse(!is.na(Tdc2) &is.na(SerT) & is.na(Tbh), "TY",
                                ifelse(!is.na(Tdc2) & !is.na(Tbh) & is.na(SerT), "OCTY", 
                                       ifelse(!is.na(Trhn), "SER",
                                              ifelse(!is.na(ple) & is.na(Trhn) & is.na(Tdc2), "DOP",
                                                     ifelse(!is.na(ple) | !is.na(Trhn) | !is.na(Tdc2), "MON","others")))))) %>%
  dplyr::arrange(celltype)

## annotate MA clusters ##

MA_celltypes <- c("MON","TY","SER","DOP", "OCTY")

for (ct in MA_celltypes) {

  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_MA_markers, exclude=anno_complete)
  
  TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
  anno_complete <- c(anno_complete, out_annotate_celltype[[2]])
  
}

#DimPlot(subset(TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled, CellType != 'Doublets'), label=T, reduction = 'tsne') + NoLegend()


## intersect with KD ##


## KD intersection for annotation

## KD (2018) annotation

df_KD_cluster_markers <- read.csv(file="Processed_Data/KD(2018)_cluster_markers.csv")
df_KD_cluster_annotation <- read.csv(file="Processed_Data/KD(2018)_cluster_annotation.csv") %>%
  dplyr::rename(cluster=Cluster_ID)
df_KD_cluster <- df_KD_cluster_markers %>%
  dplyr::left_join(., df_KD_cluster_annotation, by='cluster') %>%
  na.omit()
df_KD_cluster_select <- df_KD_cluster %>%
  dplyr::select(gene, celltype=Annotation)

TrioBrain_DmelRef_ISub_DF.integrated.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC) %>%
  dplyr::left_join(., df_KD_cluster_select, by='gene') %>%
  dplyr::ungroup() %>%
  na.omit() %>%
  dplyr::group_by(celltype, cluster) %>%
  dplyr::summarise(n=n()) %>%
  top_n(n=15, wt=n) %>%
  dplyr::ungroup() %>%
  #dplyr::arrange(celltype,-n) %>%
  ggplot(.) +
  geom_col() +
  aes(x=reorder(cluster,-n) , y=n) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=9, angle=90, hjust = 1, vjust = 0.5, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        axis.title.y=element_text(size=11, color='black'),
        axis.title.x=element_blank()) +
  facet_wrap(~celltype, scale='free', ncol=6)

#ggsave(glue::glue("Plots/TrioBrain_DmelRef_ISub_DF_KD_cluster_intersect_K50.pdf"), width=35, height =40)

## annotate Fru, clock, Poxn, OPN

fru_bru3 <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene %in% c('bru3','fru','tei','Ca-alpha1T','CG2269','jeb')) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n=n())

fru_markers <- c('bru3','fru','tei','Ca-alpha1T','CG2269','jeb')

glial_clusters <- df_glial_markers %>%
  dplyr::filter(celltype != 'others') %>%
  dplyr::pull(cluster) %>%
  as.character()
  
clock_neurons <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene %in% c('tim','Clk'), !cluster %in% glial_clusters)  %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n=n())

clock_markers <- c('tim','Dh44')

OPN_neurons <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene %in% c('SiaT','otp','acj6'), !cluster %in% glial_clusters)  %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n=n())

OPN_markers <- c('SiaT','otp','acj6')

Poxn_neurons <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene %in% c('Poxn','CG14687'), !cluster %in% glial_clusters)  %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::arrange(-n)

Poxn_markers <- c('Poxn','CG14687')

df_KD_markers <- Markers_sigs %>%
  na.omit() %>%
  filter(gene %in% c(fru_markers,OPN_markers,Poxn_markers,clock_markers), 
         !cluster %in% glial_clusters) %>%
  tidyr::spread(key = gene, value = p_val_adj) %>%
  rowwise() %>%  # Operate on each row individually
  mutate(celltype = case_when(
    sum(!is.na(c_across(all_of(fru_markers)))) == length(fru_markers) ~ "Fru",  # At least 3 markers must be non-NA
    sum(!is.na(c_across(all_of(OPN_markers)))) == length(OPN_markers) ~ "OPN",
    sum(!is.na(c_across(all_of(Poxn_markers)))) == length(Poxn_markers) ~ "Poxn",
    sum(!is.na(c_across(all_of(clock_markers)))) == length(clock_markers) ~ "Clock",
    TRUE ~ "others"
  )) %>%
  ungroup() %>%
  dplyr::filter(celltype != "others") %>%
  arrange(celltype) 

## annotate KD clusters ##

KD_celltypes <- c("Fru","OPN","Poxn","Clock")

for (ct in KD_celltypes) {
  
  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_KD_markers, exclude=anno_complete)
  
  TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
  anno_complete <- c(anno_complete, out_annotate_celltype[[2]])
  
}

## annotate based on neuropeptides ##

# load neuropeptide genes #

Dmel_gene_name <- read.table(file="genelist/Dmel_ID_genename.tsv") %>% dplyr::distinct(V1,V2) %>%
  dplyr::rename(FBID=V1, gene=V2)

df_neuropep <- data.table::fread("genelist/list_neuropeptide_FBgg0000179.txt", header = F) %>%
  dplyr::mutate(class="Neuropeptide") %>%
  dplyr::select(FBID=V1, class) %>%
  dplyr::left_join(., Dmel_gene_name, by="FBID") %>% ### few genes without matching gene_name
  na.omit()

neuropep_genes <- df_neuropep$gene

df_np_markers_specificity <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene %in% neuropep_genes) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(n_cluster=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(n_cluster)

specificity_threshold = nrow(ISub_DF_metadata_summary)*0.05

np_specific <- df_np_markers_specificity %>%
  dplyr::filter(n_cluster <= specificity_threshold) %>%
  dplyr::pull(gene) 

df_np_markers <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene %in% np_specific) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(n_sp_np = n()) %>%
  dplyr::mutate(celltype=paste(gene, collapse = "_")) %>%
  dplyr::distinct(cluster, celltype)

## annotate NP clusters ##

NP_celltypes <- unique(df_np_markers$celltype)

for (ct in NP_celltypes) {
  
  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_np_markers, exclude=anno_complete)
  
  TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
  anno_complete <- c(anno_complete, out_annotate_celltype[[2]])
  
}

## annotate remianing clusters by neurotransmitter ##

NT_genes <- c("ChAT", "VAChT","VGlut", "Gad1", "VGAT")

df_NT_markers <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene %in% c(NT_genes), !cluster %in% c(anno_complete, "Doublets")) %>%
  tidyr::spread(gene, p_val_adj) %>%
  dplyr::mutate(celltype=ifelse(!is.na(ChAT) | !is.na(VAChT), "Ach",
                                ifelse(!is.na(VGlut), "Glu",
                                       ifelse(!is.na(Gad1) | !is.na(VGAT), "GABA", "others")))) %>%
  dplyr::arrange(celltype)

## annotate NT clusters ##

NT_celltypes <- unique(df_NT_markers$celltype)

for (ct in NT_celltypes) {
  
  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_NT_markers, exclude=anno_complete)
  
  TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
  anno_complete <- c(anno_complete, out_annotate_celltype[[2]])
  
}

TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled[["ClusterLabel"]] <- Idents(object = TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled)

DimPlot(subset(TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled, CellType == 'Doublets'), label=T, reduction = 'tsne', repel=T, group.by='ClusterLabel') + NoLegend() + NoAxes()

#saveRDS(TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled, file = "Processed_Data/TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled.rds")
#TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled <- readRDS(file = "Processed_Data/TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled.rds")

DimPlot(subset(TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled, CellType != 'Doublets'), 
        label=T, reduction = 'tsne', repel=T, group.by='ClusterLabel', raster=F) + NoLegend() + NoAxes() +
  ggtitle("")
TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled

ggsave("Plots/TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled.png", width=14, height = 12)


## cluster size quantification and cell type order by cluster size ##

#TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled <- StashIdent(object = TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled, save.name = "ClusterLabel")

df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata <- TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_to_DmelRef_rep1", "Dsim_to_DmelRef_rep2", "Dsim_to_DmelRef_rep3", "Dsim_to_DmelRef_rep4", "Dsim_to_DmelRef_rep5", "Dsim_to_DmelRef_rep6"), "Dsim_to_DmelRef",
                                      ifelse(orig.ident %in% c("Dsec_to_DmelRef_rep1", "Dsec_to_DmelRef_rep2", "Dsec_to_DmelRef_rep3", "Dsec_to_DmelRef_rep4", "Dsec_to_DmelRef_rep5", "Dsec_to_DmelRef_rep6"), "Dsec",
                                             "DsecNoni_to_DmelRef")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim_to_DmelRef", "Dsec","DsecNoni_to_DmelRef")))

#final_cluster_id <- unique(c(new.cluster.ids, anno_cluster))
#unique(df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata$ClusterLabel)[!unique(df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata$ClusterLabel) %in% final_cluster_id]
#final_cluster_id %in% unique(df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata$ClusterLabel)

df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary_rep <- df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata %>%
  dplyr::group_by(species, ClusterLabel, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary <- df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary_rep %>%
  dplyr::group_by(species, ClusterLabel) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

ClusterLabel_order <- df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary %>%
  dplyr::group_by(ClusterLabel) %>%
  dplyr::summarise(mean_freq=mean(percent_combined)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-mean_freq) %>%
  dplyr::select(ClusterLabel, mean_freq) %>%
  as.vector()

ClusterLabel_order_anno <- df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary %>%
  dplyr::filter(!grepl("cluster", ClusterLabel)) %>%
  dplyr::filter(!ClusterLabel %in%c(9,19,'Doublets')) %>%
  dplyr::group_by(ClusterLabel) %>%
  dplyr::summarise(mean_freq=mean(percent_combined)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-mean_freq) %>%
  dplyr::select(ClusterLabel, mean_freq) %>%
  as.vector()

anno_cluster <- as.character(unique(ClusterLabel_order_anno$ClusterLabel))

#save(anno_cluster, ClusterLabel_order, ClusterLabel_order_anno, file="Processed_Data/celltype_order.RData")
load("Processed_Data/celltype_order.RData")

df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata$ClusterLabel <- factor(df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata$ClusterLabel, levels=ClusterLabel_order$ClusterLabel,
                                                                 labels = ClusterLabel_order$ClusterLabel)

df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary_rep <- df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata %>%
  dplyr::group_by(species, ClusterLabel, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary <- df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary_rep %>%
  dplyr::group_by(species, ClusterLabel) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary %>%
  ggplot(.) +
  geom_col(aes(y=ClusterLabel, x=percent_combined), alpha = 0.8) + 
  geom_errorbar(aes(y=ClusterLabel, x=percent_combined,
                    xmin=percent_combined-sem, xmax=percent_combined+sem), width=.2, position=position_dodge(.9)) +
  geom_point(data=df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary_rep, size = 2, alpha =0.8, shape=21,  aes(y=ClusterLabel, x=percent, fill=species)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=ClusterLabel, y=n/total*100)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11, color='black'),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=10, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        panel.grid = element_blank()) +
  labs(y="Cluster", x="Percent of cell type (%)", fill="")+
  facet_grid(~species)

#ggsave("Plots/Trio_combined_full_intron/TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_cluster_size_percent_facet.png", width=10, height = 12)

df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary %>%
  ggplot(.) +
  geom_col(aes(x=ClusterLabel, y=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.9)) + 
  geom_errorbar(aes(x=ClusterLabel, y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,0,percent_combined-sem), ymax=percent_combined+sem, group=species), width=0.9, position=position_dodge(width=0.9)) +
  geom_point(data=df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary_rep,
             size = 1, alpha =0.8,  
             aes(x=ClusterLabel, y=percent, group=species,fill=species), shape=21, position=position_dodge(width=0.9)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=ClusterLabel, y=n/total*100)) +
  theme_bw() +
  theme(axis.title.y = element_text(size=13, color='black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black', angle=90, hjust=1, vjust=0.5),
        legend.text = element_text(size=13, face='italic', color='black'),
        panel.grid = element_blank(),
        legend.position=c(0.93,0.80)) +
  labs(x="Cluster", y="Percent of cluster (%)", fill="")  +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) +
  scale_y_continuous(expand=c(0,0.1))

ggsave("Plots/TrioBrain_DmelRef_DF.integrated_cluster_size_percent_alhor.png", width=20, height = 10)

df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary %>%
  dplyr::filter(species != "DsecNoni_to_DmelRef") %>%
  ggplot(.) +
  geom_col(aes(x=species, y=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.8)) + 
  geom_errorbar(aes(x=species, y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,0,percent_combined-sem), ymax=percent_combined+sem, group=species), width=0.8, position=position_dodge(width=0.8)) +
  geom_point(data=dplyr::filter(df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary_rep, species != "DsecNoni_to_DmelRef"),
             size = 1, alpha =0.8,  
             aes(x=species, y=percent, group=species,fill=species), shape=21, position=position_dodge(width=0.8)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=ClusterLabel, y=n/total*100)) +
  theme_bw() +
  theme(axis.title.y = element_text(size=14, color='black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x=element_text(size=11, color='black', angle=45, hjust=1, vjust=1, face='italic'),
        legend.text = element_text(size=13, face='italic', color='black'),
        panel.grid = element_blank(),
        strip.text = element_text(size=12, color='black'),
        legend.position='none') +
  labs(x="Cluster", y="Frequency (% of cells)", fill="")  +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) +
  facet_wrap(~ClusterLabel, scales='free', ncol=8)

ggsave("Plots/FigS6_celltype_frequencies_all.pdf", width=13, height = 16)

for (cluster in unique(df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary$ClusterLabel)) {
  
  clustesr_namefix <-  gsub("/", "_", cluster)
  
  df_trio_reps <- df_TrioBrain_DmelRef.integrated_slim_ISub_DF_labeled_metadata_summary_rep %>%
    dplyr::filter(ClusterLabel == cluster) %>%
    dplyr::group_by(species, orig.ident) %>%
    dplyr::summarise(percent = sum(percent)) %>%
    dplyr::ungroup()
  
  df_trio_summary <- df_trio_reps %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(mean=mean(percent), sem=sd(percent)/sqrt(6)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim_to_DmelRef", "Dsec","DsecNoni_to_DmelRef")))
  
  df_trio_summary %>%
    ggplot(.) +
    geom_col(aes(y=species, x=mean, fill=species), alpha = 1, position = "dodge") + 
    geom_point(data=df_trio_reps, size = 2, alpha =0.8,  
               aes(y=species, x=percent, group=species), fill='grey', shape=21) +
    geom_errorbar(aes(y=species, x=mean,
                      xmin=ifelse(mean-sem<0,0,mean-sem), 
                      xmax=mean+sem), 
                  alpha=0.7, width=0.5) +
    #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=ClusterLabel, y=n/total*100)) +
    theme_bw() +
    theme(axis.title.x = element_text(size=12, color='black'),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=10, color='black'),
          axis.text.y=element_text(size=10, face='italic', color='black'),
          panel.grid = element_blank(),
          legend.position = 'none') +
    labs(y="Species", x=glue::glue("Percent of {cluster} cells (%)")) +
    scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3"))
  
  ggsave(glue::glue("Plots/comparison/frequencies_SCT_rpca_ISub_DF/Trio_merged_percent_{clustesr_namefix}.png"), width=3.5, height = 3.5)
  
  
}



