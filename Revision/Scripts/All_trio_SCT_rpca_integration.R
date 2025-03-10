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
Dsim_decon_combined_DF <- readRDS(file = "Processed_Data/Dsim_combined_full_intron_decon_DF.rds")
Dsec_decon_combined_DF <- readRDS(file = "Processed_Data/Dsec_combined_full_intron_decon_DF.rds")
DsecNoni_decon_combined_DF <- readRDS(file = "Processed_Data/DsecNoni_combined_full_intron_decon_DF.rds")

DefaultAssay(Dmel_decon_combined_DF) <- "SCT"
DefaultAssay(Dsim_decon_combined_DF) <- "SCT"
DefaultAssay(Dsec_decon_combined_DF) <- "SCT"
DefaultAssay(DsecNoni_decon_combined_DF) <- "SCT"

all.genes_Dmel <- rownames(x = Dmel_decon_combined_DF)
all.genes_Dsim <- rownames(x = Dsim_decon_combined_DF)
all.genes_Dsec <- rownames(x = Dsec_decon_combined_DF)
all.genes_DsecNoni <- rownames(x = DsecNoni_decon_combined_DF)

shared_genes <- all.genes_Dmel[(all.genes_Dmel %in% all.genes_Dsim) & (all.genes_Dmel %in% all.genes_Dsec) & (all.genes_Dmel %in% all.genes_DsecNoni)]
save(shared_genes,file="Processed_Data/shared_genes.RData")

DmelBrain_ortho_only_DF <- subset(Dmel_decon_combined_DF, features = shared_genes)
DsimBrain_ortho_only_DF <- subset(Dsim_decon_combined_DF, features = shared_genes)
DsecBrain_ortho_only_DF <- subset(Dsec_decon_combined_DF, features = shared_genes)
DsecNoniBrain_ortho_only_DF <- subset(DsecNoni_decon_combined_DF, features = shared_genes)

DefaultAssay(DmelBrain_ortho_only_DF) <- "SCT"
DefaultAssay(DsimBrain_ortho_only_DF) <- "SCT"
DefaultAssay(DsecBrain_ortho_only_DF) <- "SCT"
DefaultAssay(DsecNoniBrain_ortho_only_DF) <- "SCT"

TrioBrain_DF.list <- list(DmelBrain_ortho_only_DF, DsimBrain_ortho_only_DF, DsecBrain_ortho_only_DF, DsecNoniBrain_ortho_only_DF)
TrioBrain_DF.list <- lapply(X = TrioBrain_DF.list, FUN = DietSeurat, scale.data = T)

#save(TrioBrain_DF.list, file="Processed_Data/TrioBrain_DF.list.RData")

### run on cluster ###

load("Processed_Data/TrioBrain_DF.list.RData")
load("Processed_Data/shared_genes.RData")

TrioBrain.features <- SelectIntegrationFeatures(object.list = TrioBrain_DF.list, nfeatures = 3000)
TrioBrain_DF.list <- PrepSCTIntegration(object.list = TrioBrain_DF.list, assay = "SCT", anchor.features = TrioBrain.features,
                                     verbose = TRUE)
TrioBrain_DF.list <- lapply(X = TrioBrain_DF.list, FUN = RunPCA, features = TrioBrain.features)

TrioBrain.anchors <- FindIntegrationAnchors(object.list = TrioBrain_DF.list, normalization.method = "SCT",
                                            anchor.features = TrioBrain.features, verbose = T, reduction = "rpca", reference=1)

TrioBrain.integrated <- IntegrateData(anchorset = TrioBrain.anchors, normalization.method = "SCT", verbose = TRUE,
                                      features.to.integrate = shared_genes)

TrioBrain.integrated <- RunPCA(TrioBrain.integrated, npcs = 50, verbose = T)

TrioBrain.integrated <- DietSeurat(TrioBrain.integrated, scale.data=F, dimreducs = "pca")

TrioBrain.integrated <- RunUMAP(TrioBrain.integrated, dims = 1:50)
TrioBrain.integrated <- RunTSNE(TrioBrain.integrated, dims = 1:50)
TrioBrain.integrated <- FindNeighbors(TrioBrain.integrated, reduction = "pca", dims = 1:50)
TrioBrain.integrated <- FindClusters(TrioBrain.integrated, resolution = 0.2)

DimPlot(TrioBrain.integrated, label=T)

#saveRDS(TrioBrain.integrated, file = "Processed_Data/TrioBrain_combined_rpca.rds")

####

#TrioBrain.integrated <- readRDS(file = "Processed_Data/TrioBrain_combined_rpca.rds")

DefaultAssay(TrioBrain.integrated) <- "SCT"
FeaturePlot(TrioBrain.integrated, feature=c("Tret1-1"), raster=FALSE,order=T)

FeaturePlot(TrioBrain.integrated, feature=c("Gad1", "ChAT"), raster=FALSE,order=T)

### DEG analysis ###

DefaultAssay(TrioBrain.integrated) <- "RNA"

TrioBrain.integrated_slim <- DietSeurat(TrioBrain.integrated, assays = c("RNA"), dimreducs = c("umap", "tsne"))

TrioBrain.integrated_slim <- SCTransform(TrioBrain.integrated_slim, method = "glmGamPoi", verbose = FALSE)

saveRDS(TrioBrain.integrated_slim, file = "Processed_Data/TrioBrain_DF.integrated_slim.rds")
#TrioBrain.integrated_slim <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim.rds")

DimPlot(TrioBrain.integrated_slim, reduction = "tsne", label = TRUE, raster=FALSE, pt.size = 0.01)

p1 <- DimPlot(TrioBrain.integrated_slim, reduction = "tsne", label = TRUE, raster=FALSE, pt.size = 0.01) + NoLegend() + NoAxes()
p2 <- DimPlot(TrioBrain.integrated_slim, reduction = "umap", label = TRUE, raster=FALSE, pt.size = 0.01) + NoLegend() + NoAxes()
plot_grid(p1, p2)

ggsave("Plots/TrioBrain_DF_Integrated_tsne.png", width=17, height = 8)

FeaturePlot(TrioBrain.integrated_slim, feature=c("pros", "Imp", 
                                                 "ChAT", "VGlut", "Gad1"), raster=F, reduction = 'tsne', pt.size=0.01)

### plots ###

TrioBrain.integrated_slim_species <- TrioBrain.integrated_slim

TrioBrain.integrated_slim_species@meta.data$orig.ident <- gsub("_rep1", "", TrioBrain.integrated_slim_species@meta.data$orig.ident)
TrioBrain.integrated_slim_species@meta.data$orig.ident <- gsub("_rep2", "", TrioBrain.integrated_slim_species@meta.data$orig.ident)
TrioBrain.integrated_slim_species@meta.data$orig.ident <- gsub("_rep3", "", TrioBrain.integrated_slim_species@meta.data$orig.ident)
TrioBrain.integrated_slim_species@meta.data$orig.ident <- gsub("_rep4", "", TrioBrain.integrated_slim_species@meta.data$orig.ident)
TrioBrain.integrated_slim_species@meta.data$orig.ident <- gsub("_rep5", "", TrioBrain.integrated_slim_species@meta.data$orig.ident)
TrioBrain.integrated_slim_species@meta.data$orig.ident <- gsub("_rep6", "", TrioBrain.integrated_slim_species@meta.data$orig.ident)

TrioBrain.integrated_slim_species@meta.data$orig.ident <- factor(TrioBrain.integrated_slim_species@meta.data$orig.ident, 
                                                                 levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))

plot_dmel_tsne <- DimPlot(subset(TrioBrain.integrated_slim_species, orig.ident == "Dmel"), 
                          reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01,
                          cols=c("#E41A1C")) + NoLegend() + labs(title=NULL) + NoAxes()

plot_dsim_tsne <- DimPlot(subset(TrioBrain.integrated_slim_species, orig.ident == "Dsim"), 
                          reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01,
                          cols=c("#4DAF4A")) + NoLegend() + labs(title=NULL) + NoAxes()

plot_dsec_tsne <- DimPlot(subset(TrioBrain.integrated_slim_species, orig.ident == "Dsec"), 
                          reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01,
                          cols=c("#377EB8")) + NoLegend() + labs(title=NULL) + NoAxes()

plot_trio_merged_tsne <- DimPlot(subset(TrioBrain.integrated_slim_species, orig.ident != "DsecNoni"), 
                                 reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01,
                                 cols=c("#E41A1C", "#4DAF4A",  "#377EB8")) + NoLegend() + labs(title=NULL) + NoAxes()

plot_tsne_species <- cowplot::plot_grid(plot_dmel_tsne, plot_dsim_tsne,
                                        plot_dsec_tsne, plot_trio_merged_tsne, ncol=2)

plot_tsne_species

ggsave("Plots/Manuscript/Fig1_tsne_species.png", width=13, height = 12)

#DimPlot(TrioBrain.integrated_slim_species, reduction = "umap", shuffle = T, raster=FALSE, label=T, pt.size = 0.01) + NoLegend() + labs(title=NULL)
#ggsave("Plots/TrioBrain_DF_Integrated_umap_label.png", width=11, height = 10)

p1 <- DimPlot(TrioBrain.integrated_slim_species, reduction = "umap", 
              group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01, cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL)
p2 <- DimPlot(TrioBrain.integrated_slim_species, reduction = "umap", label = F, raster=F, pt.size = 0.01) + NoLegend()
plot_grid(p1, p2)

#ggsave("Plots/TrioBrain_DF_Integrated_umap.png", width=23, height = 10)

### marker expression ###

FeaturePlot(TrioBrain.integrated_slim, feature = 'repo', order=T, reduction = "tsne")
FeaturePlot(TrioBrain.integrated_slim, feature = 'SK', order=T, reduction = "tsne")

Glia <- WhichCells(TrioBrain.integrated_slim, expression = repo>0)
MAcells <- WhichCells(TrioBrain.integrated_slim, expression = Vmat>3)
DOPcells <- WhichCells(TrioBrain.integrated_slim, expression = ple>1)
AchCells <- WhichCells(TrioBrain.integrated_slim, expression = VAChT>1 & VGlut<=1.5 & Gad1<=1.2)
GluCells <- WhichCells(TrioBrain.integrated_slim, expression = VGlut>1.5 & Gad1<=1.2 & VAChT<=1)
GabaCells <- WhichCells(TrioBrain.integrated_slim, expression = Gad1>1.2 & VGlut<=1.5 & VAChT<=1)
MultiCells <- WhichCells(TrioBrain.integrated_slim, expression = Gad1>1.2 & VAChT>1 | VAChT>1 &VGlut>1.5 | VGlut>1.5 &Gad1>1.2)

p1_type <- DimPlot(TrioBrain.integrated_slim, label=F,shuffle=T,
                   cells.highlight= list(VAChT=AchCells, VGlut=GluCells, Gad1=GabaCells, Repo=Glia, Vmat=MAcells), 
                   reduction = "tsne", 
                   pt.size = 0.01,sizes.highlight = 0.01,
                   cols.highlight = c("darkorange", "darkgreen", "darkblue","pink", "darkred"), cols= "grey40", raster=FALSE) + NoLegend() + NoAxes()

p2_type <- DimPlot(TrioBrain.integrated_slim, label=F,shuffle=T,
                   cells.highlight= list(VAChT=AchCells, VGlut=GluCells, Gad1=GabaCells, Repo=Glia, Vmat=MAcells), 
                   reduction = "umap", 
                   pt.size = 0.01,sizes.highlight = 0.01,
                   cols.highlight = c("darkorange", "darkgreen", "darkblue","pink", "darkred"), cols= "grey40", raster=FALSE) + NoLegend() + NoAxes()

plot_grid(p1_type, p2_type)
ggsave("Plots/Manuscript/TrioBrain_DF_Integrated_tsne_UMAP_nt.png", width=17, height = 8)

DimPlot(TrioBrain.integrated_slim, label=F,shuffle=T,
        cells.highlight= list(VAChT=AchCells, VGlut=GluCells, Gad1=GabaCells, Repo=Glia, Vmat=MAcells), 
        reduction = "umap", 
        pt.size = 0.01,sizes.highlight = 0.01,
        cols.highlight = c("darkorange", "darkgreen", "darkblue","black", "darkred"), cols= "grey40", raster=T) + NoAxes()

ggsave("Plots/TrioBrain_DF_Integrated_umap_nt_for_label.pdf", width=9, height = 7)


p3 <- DimPlot(TrioBrain.integrated_slim, label=F,shuffle=T,
              cells.highlight= list(VAChT=AchCells, VGlut=GluCells, Gad1=GabaCells, Repo=Glia, Vmat=MAcells), 
              reduction = "umap", 
              pt.size = 0.01,sizes.highlight = 0.01, 
              cols.highlight = c("darkorange", "darkgreen", "darkblue","black", "darkred"), cols= "grey40", raster=F) +NoLegend()

plot_grid(p1, p3, p2, nrow=1)

ggsave("Plots/Manuscript/Fig1a_raster.pdf", width=30, height = 10)

p1_lab <- DimPlot(TrioBrain.integrated_slim_species, reduction = "umap", 
                  group.by = "orig.ident", shuffle = T, raster=T, pt.size = 0.01, cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + labs(title=NULL)

p2_lab <- DimPlot(TrioBrain.integrated_slim_species, reduction = "umap", label = F, raster=T, pt.size = 0.01)

p3_lab <- DimPlot(TrioBrain.integrated_slim, label=F,shuffle=T,
                  cells.highlight= list(VAChT=AchCells, VGlut=GluCells, Gad1=GabaCells, Repo=Glia, Vmat=MAcells), 
                  reduction = "umap", 
                  pt.size = 0.01,sizes.highlight = 0.01, 
                  cols.highlight = c("darkorange", "darkgreen", "darkblue","black", "darkred"), cols= "grey40", raster=T)

plot_grid(p1_lab, p3_lab, p2_lab, nrow=1)

ggsave("Plots/Manuscript/Fig1a_raster_lab.pdf", width=30, height = 10)

## pros, Imp (master regulator) ##

ImpCells <- WhichCells(TrioBrain.integrated_slim, expression = Imp>1.2)
Proscells <- WhichCells(TrioBrain.integrated_slim, expression = pros>2)

p1_mr <- DimPlot(TrioBrain.integrated_slim, label=F,shuffle=T,
                 cells.highlight= list(Pros=Proscells, Imp=ImpCells), 
                 reduction = "tsne", 
                 pt.size = 0.01,sizes.highlight = 0.01,
                 cols.highlight = c("darkcyan", "darkmagenta"), cols= "grey40", raster=FALSE) + NoLegend() + NoAxes()

p2_mr <- DimPlot(TrioBrain.integrated_slim, label=F,shuffle=T,
                 cells.highlight= list(Pros=Proscells, Imp=ImpCells), 
                 reduction = "umap", 
                 pt.size = 0.01,sizes.highlight = 0.01,
                 cols.highlight = c("darkcyan", "darkmagenta"), cols= "grey40", raster=FALSE) + NoLegend() + NoAxes()

plot_grid(p1_mr, p2_mr)
ggsave("Plots/Manuscript/TrioBrain_DF_Integrated_tsne_UMAP_mr.png", width=17, height = 8)

DimPlot(TrioBrain.integrated_slim, label=F,shuffle=T,
        cells.highlight= list(Pros=Proscells, Imp=ImpCells), 
        reduction = "umap", 
        pt.size = 0.01,sizes.highlight = 0.01,
        cols.highlight = c("darkcyan", "darkmagenta"), cols= "grey40", raster=T) + NoAxes()

ggsave("Plots/TrioBrain_DF_Integrated_umap_mr_for_label.pdf", width=9, height = 7)

### Iterative subclustering - see Interative_subclustering.r ###

TrioBrain.integrated_slim_ISub_DF <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_ISub_DF.rds")

ClusterIdents <- unique(Idents(TrioBrain.integrated_slim_ISub_DF))

DimPlot(TrioBrain.integrated_slim_ISub_DF, reduction='umap', label=T) + NoLegend()

## QC and additional doublet filtering ##

clusters_15 <- ClusterIdents[grep(15, ClusterIdents)] ## PRN
VlnPlot(subset(TrioBrain.integrated_slim_ISub_DF, idents = clusters_15), features = c("nFeature_RNA")) + NoLegend()

clusters_7 <- as.character(ClusterIdents[grep('cluster7', ClusterIdents)]) ## ENS
VlnPlot(subset(TrioBrain.integrated_slim_ISub_DF, idents = c(7,clusters_7)), features = c("nFeature_RNA")) + NoLegend()

clusters_0 <- as.character(ClusterIdents[grep('cluster0', ClusterIdents)])
VlnPlot(subset(TrioBrain.integrated_slim_ISub_DF, idents = c(0,clusters_0)), features = c("nFeature_RNA")) + NoLegend()

## nFeature_RNA

cluster_features_nFeature_RNA <- TrioBrain.integrated_slim_ISub_DF@meta.data %>%
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

cluster_features_nCount_RNA <- TrioBrain.integrated_slim_ISub_DF@meta.data %>%
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

#TrioBrain_ISub_DF.integrated.markers <- FindAllMarkers(TrioBrain.integrated_slim_ISub_DF, only.pos = TRUE, min.pct = 0.15, 
#                                               logfc.threshold = 0.25, test.use = "MAST")

plan(sequential)  # Reset to default plan
gc()  # Cleanup unused memory

#saveRDS(TrioBrain_ISub_DF.integrated.markers, file = "Processed_Data/TrioBrain_ISub_DF.integrated.markers.rds")
TrioBrain_ISub_DF.integrated.markers <- readRDS(file = "Processed_Data/TrioBrain_ISub_DF.integrated.markers.rds")

### annotate based on markeres ###

TrioBrain.integrated_slim_ISub_DF_labeled <- TrioBrain.integrated_slim_ISub_DF

total_cell_number <- nrow(TrioBrain.integrated_slim_ISub_DF@meta.data)

ISub_DF_metadata_summary <- TrioBrain.integrated_slim_ISub_DF@meta.data %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-n) %>%
  dplyr::mutate(percent = n/total_cell_number*100)

anno_complete <- NULL

annotate_celltype <- function(seurat_object=TrioBrain.integrated_slim_ISub_DF_labeled,
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

Markers_sigs <- TrioBrain_ISub_DF.integrated.markers %>%
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

#clu14_markers <- TrioBrain_ISub_DF.integrated.markers %>%  dplyr::filter(cluster ==14)

## annotate glial clusters ##

glial_celltypes <- c("ENS","AST","PRN","SUB", "CTX")

anno_complete <- NULL

for (ct in glial_celltypes) {
  
  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_glial_markers)
  
  TrioBrain.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
  anno_complete <- c(anno_complete, out_annotate_celltype[[2]])
}

#DimPlot(subset(TrioBrain.integrated_slim_ISub_DF_labeled, idents = c("ENS_1","ENS_2","ENS_3", "PRN", "SUB", "AST", "CTX")))

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
  
  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_KC_markers, exclude=anno_complete)
  
  TrioBrain.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
  anno_complete <- c(anno_complete, out_annotate_celltype[[2]])
  
}

#DimPlot(subset(TrioBrain.integrated_slim_ISub_DF_labeled, idents = c("αβ-KC_1","αβ-KC_2","αβ-KC_3","γ-KC_1","γ-KC_2","γ-KC_3","γ-KC_4","α'β'-KC")))

#FeaturePlot(subset(TrioBrain.integrated_slim_ISub_DF_labeled, idents = c("αβ-KC_1","αβ-KC_2","αβ-KC_3","γ-KC_1","γ-KC_2","γ-KC_3","γ-KC_4","α'β'-KC")),feature='CG33098')

## Monoaminergic neurons ##

#MA_markers <- c("Vmat", "ple", "DAT","Trhn","SerT","Tdc2","Tbh")
MA_marekrs <- "Vmat"

df_MA_markers <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene == "Vmat") %>%
  tidyr::spread(gene, p_val_adj) %>%
  dplyr::filter(!is.na(Vmat)) %>%
  dplyr::mutate(celltype=ifelse(!is.na(Vmat), "MON","others")) %>%
  dplyr::arrange(celltype)

## annotate MA clusters ##

#MA_celltypes <- c("MON","TY","SER","DOP", "OCTY")
MA_celltypes <- "MON"
 
for (ct in MA_celltypes) {

  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_MA_markers, exclude=anno_complete)
  
  TrioBrain.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
  anno_complete <- c(anno_complete, out_annotate_celltype[[2]])
  
}

#DimPlot(subset(TrioBrain.integrated_slim_ISub_DF_labeled, CellType != 'Doublets'), label=T, reduction = 'tsne') + NoLegend()


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

TrioBrain_ISub_DF.integrated.markers %>% 
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

#ggsave(glue::glue("Plots/TrioBrain_ISub_DF_KD_cluster_intersect_K50.pdf"), width=35, height =40)

## Clock neurons ##

Clock_marekrs <- c("Clk","tim")

df_Clock_markers <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene %in% Clock_marekrs) %>%
  tidyr::spread(gene, p_val_adj) %>%
  dplyr::filter(!is.na(Clk) | !is.na(tim)) %>%
  dplyr::mutate(celltype=ifelse(!is.na(Clk) | !is.na(tim), "Clock","others")) %>%
  dplyr::arrange(celltype)

## annotate Clock clusters ##

Clock_celltypes <- "Clock"

for (ct in Clock_celltypes) {
  
  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_Clock_markers, exclude=anno_complete)
  
  TrioBrain.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
  anno_complete <- c(anno_complete, out_annotate_celltype[[2]])
  
}

#DimPlot(subset(TrioBrain.integrated_slim_ISub_DF_labeled, CellType != 'Doublets'), label=T, reduction = 'tsne') + NoLegend() + NoAxes()


## OPN neurons ##

OPN_markers <- c('otp','acj6')

df_OPN_markers <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene %in% OPN_markers) %>%
  tidyr::spread(gene, p_val_adj) %>%
  dplyr::filter(!is.na(acj6) | !is.na(otp)) %>%
  dplyr::mutate(celltype=ifelse(!is.na(otp) & !is.na(acj6), "OPN","others")) %>%
  dplyr::arrange(celltype)

## annotate OPN clusters ##

OPN_celltypes <- "OPN"

for (ct in OPN_celltypes) {
  
  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_OPN_markers, exclude=anno_complete)
  
  TrioBrain.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
  anno_complete <- c(anno_complete, out_annotate_celltype[[2]])
  
}

#DimPlot(subset(TrioBrain.integrated_slim_ISub_DF_labeled, CellType != 'Doublets'), label=T, reduction = 'tsne') + NoLegend() + NoAxes()


## Poxn neurons ##

Poxn_markers <- 'Poxn'

df_Poxn_markers <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene == "Poxn") %>%
  tidyr::spread(gene, p_val_adj) %>%
  dplyr::filter(!is.na(Poxn)) %>%
  dplyr::mutate(celltype=ifelse(!is.na(Poxn), "Poxn","others")) %>%
  dplyr::arrange(celltype)

## annotate Poxn clusters ##

Poxn_celltypes <- "Poxn"

for (ct in Poxn_celltypes) {
  
  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_Poxn_markers, exclude=anno_complete)
  
  TrioBrain.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
  anno_complete <- c(anno_complete, out_annotate_celltype[[2]])
  
}

#DimPlot(subset(TrioBrain.integrated_slim_ISub_DF_labeled, CellType != 'Doublets'), label=T, reduction = 'tsne') + NoLegend() + NoAxes()


## fru neurons ##

fru_markers <- 'fru'

df_fru_markers <- Markers_sigs %>%
  na.omit() %>%
  dplyr::filter(gene == "fru") %>%
  tidyr::spread(gene, p_val_adj) %>%
  dplyr::filter(!is.na(fru)) %>%
  dplyr::mutate(celltype=ifelse(!is.na(fru), "fru","others")) %>%
  dplyr::arrange(celltype)

## annotate fru clusters ##

fru_celltypes <- "fru"

for (ct in fru_celltypes) {
  
  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_fru_markers, exclude=anno_complete)
  
  TrioBrain.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
  anno_complete <- c(anno_complete, out_annotate_celltype[[2]])
  
}

#DimPlot(subset(TrioBrain.integrated_slim_ISub_DF_labeled, CellType != 'Doublets'), label=T, reduction = 'tsne') + NoLegend() + NoAxes()

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
  dplyr::filter(n_cluster < specificity_threshold) %>%
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

save(np_specific, NP_celltypes, file="Processed_Data/NP_celltypes.RData")

for (ct in NP_celltypes) {
  
  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_np_markers, exclude=anno_complete)
  
  TrioBrain.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
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
  
  out_annotate_celltype <- annotate_celltype(seurat_object=TrioBrain.integrated_slim_ISub_DF_labeled,
                                             cell_type=ct, df_markers = df_NT_markers, exclude=anno_complete)
  
  TrioBrain.integrated_slim_ISub_DF_labeled <- out_annotate_celltype[[1]]
  
  anno_complete <- c(anno_complete, out_annotate_celltype[[2]])
  
}

TrioBrain.integrated_slim_ISub_DF_labeled[["ClusterLabel"]] <- Idents(object = TrioBrain.integrated_slim_ISub_DF_labeled)

DimPlot(subset(TrioBrain.integrated_slim_ISub_DF_labeled, CellType != 'Doublets'), label=T, reduction = 'tsne', repel=T, group.by='ClusterLabel') + NoLegend() + NoAxes()

saveRDS(TrioBrain.integrated_slim_ISub_DF_labeled, file = "Processed_Data/TrioBrain.integrated_slim_ISub_DF_labeled.rds")
#TrioBrain.integrated_slim_ISub_DF_labeled <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_ISub_DF_labeled.rds")

DimPlot(subset(TrioBrain.integrated_slim_ISub_DF_labeled, CellType != 'Doublets'), 
        label=T, reduction = 'tsne', repel=T, group.by='ClusterLabel', raster=F) + NoLegend() + NoAxes() +
  ggtitle("")
TrioBrain.integrated_slim_ISub_DF_labeled

ggsave("Plots/TrioBrain.integrated_slim_ISub_DF_labeled.png", width=14, height = 12)
