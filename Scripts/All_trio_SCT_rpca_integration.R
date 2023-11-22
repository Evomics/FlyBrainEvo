library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(sctransform)
library(cowplot)
library(broom)
library(MAST)
library(enrichR)
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

save(TrioBrain_DF.list, file="Processed_Data/TrioBrain_DF.list.RData")

### run on cluster ###

load("/work/FAC/FBM/CIG/rbenton/neuroflies/Daehan/FlyBrainEvo/R/TrioBrain.list.RData")
load("/work/FAC/FBM/CIG/rbenton/neuroflies/Daehan/FlyBrainEvo/R/shared_genes.RData")

TrioBrain.features <- SelectIntegrationFeatures(object.list = TrioBrain.list, nfeatures = 3000)
TrioBrain.list <- PrepSCTIntegration(object.list = TrioBrain.list, assay = "SCT", anchor.features = TrioBrain.features,
                                     verbose = TRUE)
TrioBrain.list <- lapply(X = TrioBrain.list, FUN = RunPCA, features = TrioBrain.features)

TrioBrain.anchors <- FindIntegrationAnchors(object.list = TrioBrain.list, normalization.method = "SCT",
                                            anchor.features = TrioBrain.features, verbose = T, reduction = "rpca", reference=1)

TrioBrain.integrated <- IntegrateData(anchorset = TrioBrain.anchors, normalization.method = "SCT", verbose = TRUE,
                                      features.to.integrate = shared_genes)

TrioBrain.integrated <- RunPCA(TrioBrain.integrated, npcs = 50, verbose = T)

TrioBrain.integrated <- DietSeurat(TrioBrain.integrated, scale.data=F, dimreducs = "pca")

TrioBrain.integrated <- RunUMAP(TrioBrain.integrated, dims = 1:50)
TrioBrain.integrated <- RunTSNE(TrioBrain.integrated, dims = 1:50)
TrioBrain.integrated <- FindNeighbors(TrioBrain.integrated, reduction = "pca", dims = 1:50)
TrioBrain.integrated <- FindClusters(TrioBrain.integrated, resolution = 0.2)

saveRDS(TrioBrain.integrated, file = "/work/FAC/FBM/CIG/rbenton/neuroflies/Daehan/FlyBrainEvo/R/TrioBrain_combined_rpca.rds")

####

#saveRDS(TrioBrain.integrated, file = "Processed_Data/TrioBrain_DF_combined_rpca.rds")

#TrioBrain.integrated <- readRDS(file = "Processed_Data/TrioBrain_combined_rpca.rds")
TrioBrain.integrated <- readRDS(file = "Processed_Data/TrioBrain_DF_combined_rpca.rds")

p1 <- DimPlot(UpdateSeuratObject(TrioBrain.integrated), reduction = "umap", group.by = "orig.ident", shuffle = T, raster=FALSE, pt.size = 0.01)
p2 <- DimPlot(TrioBrain.integrated, reduction = "umap", label = TRUE, raster=FALSE, pt.size = 0.01)
plot_grid(p1, p2)

ggsave("Plots/TrioBrain_DF_Integrated_UMAP.png", width=17, height = 7)

DefaultAssay(TrioBrain.integrated) <- "SCT"
FeaturePlot(TrioBrain.integrated, feature=c("Tret1-1"), raster=FALSE,order=T)

FeaturePlot(TrioBrain.integrated, feature=c("Gad1", "ChAT"), raster=FALSE,order=T)

### DEG analysis ###

DefaultAssay(TrioBrain.integrated) <- "RNA"

TrioBrain.integrated_slim <- DietSeurat(TrioBrain.integrated, assays = c("RNA"), dimreducs = c("umap", "tsne"))

#saveRDS(TrioBrain.integrated_slim, file = "Processed_Data/TrioBrain_DF.integrated_slim.rds")

TrioBrain.integrated_slim <- SCTransform(TrioBrain.integrated_slim, method = "glmGamPoi", verbose = FALSE)

saveRDS(TrioBrain.integrated_slim, file = "Processed_Data/TrioBrain_DF.integrated_slim.rds")
TrioBrain.integrated_slim <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim.rds")

DimPlot(TrioBrain.integrated_slim, reduction = "tsne", label = TRUE, raster=FALSE, pt.size = 0.01)

p1 <- DimPlot(TrioBrain.integrated_slim, reduction = "tsne", label = TRUE, raster=FALSE, pt.size = 0.01) + NoLegend() + NoAxes()
p2 <- DimPlot(TrioBrain.integrated_slim, reduction = "umap", label = TRUE, raster=FALSE, pt.size = 0.01) + NoLegend() + NoAxes()
plot_grid(p1, p2)

ggsave("Plots/TrioBrain_DF_Integrated_tsne_UMAP.png", width=17, height = 8)


FeaturePlot(TrioBrain.integrated_slim, feature="dpr6", raster=F, reduction = 'tsne')

FeaturePlot(TrioBrain.integrated_slim, feature=c("CAH2", "CAH3"), raster=F, reduction = 'tsne', order=T)



FeaturePlot(TrioBrain.integrated_slim, feature=c("pros", "Imp", 
                                                 "ChAT", "VGlut", "Gad1"), raster=F, reduction = 'tsne', pt.size=0.01)

##


TrioBrain.integrated.markers <- FindAllMarkers(TrioBrain.integrated_slim, only.pos = TRUE, min.pct = 0.15, 
                                               logfc.threshold = 0.25, test.use = "MAST")

saveRDS(TrioBrain.integrated.markers, file = "Processed_Data/TrioBrain_DF.integrated.markers.rds")
TrioBrain.integrated.markers <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated.markers.rds")

top10_Brain <- TrioBrain.integrated.markers %>% group_by(cluster) %>% top_n(n = 10, wt = abs(avg_log2FC))

#DoHeatmap(TrioBrain.integrated_slim, features = top10_Brain$gene, slot="data") + NoLegend()
#ggsave("Plots/TrioBrain_DF_marekr_heatmap.png", width=25, height = 18)

### 24 markers for all clusters ###

top24_TrioBrain.integrated <- TrioBrain.integrated.markers %>% 
  group_by(cluster) %>% 
  #dplyr::filter(pct.2 < 0.2) %>%
  top_n(n = 24, wt = -p_val) %>% top_n(n = 24, wt = avg_log2FC)

saveRDS(top24_TrioBrain.integrated, file = "Processed_Data/top24_TrioBrain_DF.integrated.rds")
#top24_TrioBrain.integrated <- readRDS(file = "Processed_Data/top24_TrioBrain_DF.integrated.rds")

for (i in unique(top24_TrioBrain.integrated$cluster)) {
  
  top16_TrioBrain.integrated_cluster <- top24_TrioBrain.integrated %>%
    dplyr::filter(cluster == i)
  
  VlnPlot(TrioBrain.integrated_slim, features = top16_TrioBrain.integrated_cluster$gene, 
          pt.size = 0, ncol=3, assay = "SCT")
  
  ggsave(glue::glue("Plots/Markers/All_DF02/TrioBrain_markers_clusters{i}.png"), width=24, height = 24)
  
}


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

TrioBrain.integrated.markers %>% 
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

ggsave(glue::glue("Plots/TrioBrain_DF_KD_cluster_intersect_K50.png"), width=25, height =40)

#### Module scores ###

## KC ##

KC_markers_KD <- df_KD_cluster_markers %>%
  dplyr::left_join(., df_KD_cluster_annotation, by='cluster') %>%
  na.omit() %>%
  dplyr::filter(grepl("KC", Annotation)) %>%
  #dplyr::filter(gene %in% Trio.combined.markers_filter$gene)  %>%
  dplyr::group_by(Annotation) %>%
  top_n(n = 50, wt = -p.value) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene, celltype=Annotation)

Kab_markers <- dplyr::filter(KC_markers_KD, grepl("a/b", celltype))$gene
Kapbp_markers <- dplyr::filter(KC_markers_KD, grepl("a'/b'", celltype))$gene
Kg_markers <- dplyr::filter(KC_markers_KD, grepl("G-KC", celltype))$gene

TrioBrain.integrated_slim <- AddModuleScore(object = TrioBrain.integrated_slim, features = list(Kab_markers), name = 'ab')
TrioBrain.integrated_slim <- AddModuleScore(object = TrioBrain.integrated_slim, features = list(Kg_markers), name = 'g')
TrioBrain.integrated_slim <- AddModuleScore(object = TrioBrain.integrated_slim, features = list(Kapbp_markers), name = 'apbp')

VlnPlot(TrioBrain.integrated_slim, features=c("ab1","g1", "apbp1"), pt.size=0) + NoLegend()

ggsave("Plots/TrioBrain_DF_KC_module_score.png", width=25, height=9)


## glia ##

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

TrioBrain.integrated_slim <- AddModuleScore(object = TrioBrain.integrated_slim, features = list(glia_markers), name = 'Glia')
TrioBrain.integrated_slim <- AddModuleScore(object = TrioBrain.integrated_slim, features = list(ENS_markers), name = 'ens')
TrioBrain.integrated_slim <- AddModuleScore(object = TrioBrain.integrated_slim, features = list(AST_markers), name = 'ast')
TrioBrain.integrated_slim <- AddModuleScore(object = TrioBrain.integrated_slim, features = list(PRN_markers), name = 'prn')
TrioBrain.integrated_slim <- AddModuleScore(object = TrioBrain.integrated_slim, features = list(SUB_markers), name = 'sub')
TrioBrain.integrated_slim <- AddModuleScore(object = TrioBrain.integrated_slim, features = list(CTX_markers), name = 'ctx')
TrioBrain.integrated_slim <- AddModuleScore(object = TrioBrain.integrated_slim, features = list(CHI_markers), name = 'chi')

VlnPlot(TrioBrain.integrated_slim, features=c("Glia1"), pt.size=0) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank()) +
  NoLegend() 

ggsave("Plots/TrioBrain_DF_glia_module_score.png", width=10, height = 4)

VlnPlot(TrioBrain.integrated_slim, features=c("ens1","ast1", "prn1","sub1","ctx1","chi1"), pt.size=0) *
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_blank()) +
  NoLegend() 

ggsave("Plots/TrioBrain_glia_subtype_module_scores.png", width=35, height = 15)


###

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


DimPlot(TrioBrain.integrated_slim_species, reduction = "umap", shuffle = T, raster=FALSE, label=T, pt.size = 0.01) + NoLegend() + labs(title=NULL)
ggsave("Plots/TrioBrain_DF_Integrated_umap_label.png", width=11, height = 10)

p1 <- DimPlot(TrioBrain.integrated_slim_species, reduction = "umap", 
              group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01, cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL)
p2 <- DimPlot(TrioBrain.integrated_slim_species, reduction = "umap", label = F, raster=F, pt.size = 0.01) + NoLegend()
plot_grid(p1, p2)

ggsave("Plots/TrioBrain_DF_Integrated_umap.png", width=23, height = 10)

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
        cols.highlight = c("darkorange", "darkgreen", "darkblue","black", "darkred"), cols= "grey40", raster=FALSE) + NoLegend() + NoAxes()

p2_type <- DimPlot(TrioBrain.integrated_slim, label=F,shuffle=T,
                   cells.highlight= list(VAChT=AchCells, VGlut=GluCells, Gad1=GabaCells, Repo=Glia, Vmat=MAcells), 
                   reduction = "umap", 
                   pt.size = 0.01,sizes.highlight = 0.01,
                   cols.highlight = c("darkorange", "darkgreen", "darkblue","black", "darkred"), cols= "grey40", raster=FALSE) + NoLegend() + NoAxes()

plot_grid(p1_type, p2_type)
ggsave("Plots/TrioBrain_DF_Integrated_tsne_UMAP_nt.png", width=17, height = 8)

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
ggsave("Plots/TrioBrain_DF_Integrated_tsne_UMAP_mr.png", width=17, height = 8)


DimPlot(TrioBrain.integrated_slim, label=F,shuffle=T,
        cells.highlight= list(Pros=Proscells, Imp=ImpCells), 
        reduction = "umap", 
        pt.size = 0.01,sizes.highlight = 0.01,
        cols.highlight = c("darkcyan", "darkmagenta"), cols= "grey40", raster=T) + NoAxes()

ggsave("Plots/TrioBrain_DF_Integrated_umap_mr_for_label.pdf", width=9, height = 7)


## label cluster ###

FeaturePlot(TrioBrain.integrated_slim, feature=c("AstA", "CCHa2", "tim", "Clk", "cry"), order=T)

FeaturePlot(TrioBrain.integrated_slim, feature=c("Dh31", "sNPF", "fru", "Lmx1a"))

new.cluster.ids <- c("Cluster0", "Cluster1", "GABA", 
                     "αβ-KC", "γ-KC",
                     "ENS", "Cluster6", "Cluster7", "AST",
                     "MON", "α'β'-KC",
                     "fru(E)", "sNPF(E)", "Dh31(E)", "Cluster14",
                     "Cluster15", "Da6(E)","Cluster17","Poxn",
                     "Cluster19", "Clock", "CCHa2(1)", "Cluster22", 
                     "Cluster23","CTX", "Da5/7","Dh44","Cluster27",
                     "Cluster28","Cluster29","Dh31(I2)","ort", "AstA","Cluster33",
                     "RYa", "Cluster35", "SIFa", "Cluster37")

names(new.cluster.ids) <- levels(TrioBrain.integrated_slim)
TrioBrain.integrated_slim_labeled <- RenameIdents(TrioBrain.integrated_slim, new.cluster.ids)

DimPlot(TrioBrain.integrated_slim_labeled, reduction = "umap", label = TRUE, pt.size = 0.01) + NoLegend()

### Add subclustser identity ###

load("Processed_Data/Cluster0.RData")
load("Processed_Data/Cluster1.RData")
load("Processed_Data/Cluster2.RData")
load("Processed_Data/Cluster7.RData")
load("Processed_Data/Cluster9.RData")
load("Processed_Data/Cluster11.RData")
load("Processed_Data/Cluster12.RData")
load("Processed_Data/Cluster14.RData")
load("Processed_Data/Cluster15.RData")
load("Processed_Data/Cluster17.RData")

### Cluster example ###

DimPlot(Cluster9_seurat, reduction = 'umap', label=T) + NoLegend() + NoAxes()
ggsave("Plots/Manuscript/FigS1_MON_subclustering.pdf", width=6, height = 5)

FeaturePlot(Cluster9_seurat, feature=c("ple", "SerT", "Tdc2", "Tbh"), ncol=2, pt.size=0.2) * NoAxes()
ggsave("Plots/Manuscript/FigS1_MON_subclustering_exp.pdf", width=8, height = 5)


DimPlot(Cluster17_seurat, reduction = 'umap', label=T) + NoLegend() + NoAxes()
ggsave("Plots/Manuscript/FigS1_PEP_subclustering.pdf", width=5.5, height = 5)

FeaturePlot(Cluster17_seurat, feature=c("Hug", "Ilp2", "Ilp3", "Tk", "Mip", "Crz"), ncol=3, pt.size=0.2) * NoAxes()
ggsave("Plots/Manuscript/FigS1_PEP_subclustering_exp.pdf", width=8.5, height = 5)


FMRFa1_cellid <- WhichCells(Cluster0_seurat, idents=5)

OPN_cellid <- WhichCells(Cluster1_seurat, idents=3)
Proc_cellid <- WhichCells(Cluster1_seurat, idents=4)

Dh31I1_cellid <- WhichCells(Cluster2_seurat, idents=3)
Da6I_cellid <- WhichCells(Cluster2_seurat, idents=8)
fruI_cellid <- WhichCells(Cluster2_seurat, idents=7)
sNPFI_cellid <- WhichCells(Cluster2_seurat, idents=11)

Mip1_cellid <- WhichCells(Cluster7_seurat, idents=2)

DOP_cellid <- WhichCells(Cluster9_seurat, idents=c(0,1))
SER_cellid <- WhichCells(Cluster9_seurat, idents=2)
OCTY_cellid <- WhichCells(Cluster9_seurat, idents=c(3,6))
TY_cellid <- WhichCells(Cluster9_seurat, idents=c(4,5))
Tbh_cellid <- WhichCells(Cluster9_seurat, idents=c(8))

fruGlu_cellid <- WhichCells(Cluster11_seurat, idents=1)
fruMs_cellid <- WhichCells(Cluster11_seurat, idents=3)

CCHa22_cellid <- WhichCells(Cluster12_seurat, idents=2)
FMRFa2_cellid <- WhichCells(Cluster12_seurat, idents=4)

PRN_cellid <- WhichCells(Cluster14_seurat, idents=c(0,1))
SUB_cellid <- WhichCells(Cluster14_seurat, idents=c(4))
Blood_cellid <- WhichCells(Cluster14_seurat, idents=c(6))

Fat_cellid <- WhichCells(Cluster15_seurat, idents=c(3))

Hug_cellid <- WhichCells(Cluster17_seurat, idents=c(2))
IPC_cellid <- WhichCells(Cluster17_seurat, idents=c(3))
Tk_cellid <- WhichCells(Cluster17_seurat, idents=c(4))
Mip2_cellid <- WhichCells(Cluster17_seurat, idents=c(5))
Crz_cellid <- WhichCells(Cluster17_seurat, idents=c(6))

Idents(TrioBrain.integrated_slim_labeled, cells = FMRFa1_cellid) <- "FMRFa(1)"

Idents(TrioBrain.integrated_slim_labeled, cells = OPN_cellid) <- "OPN"
Idents(TrioBrain.integrated_slim_labeled, cells = Proc_cellid) <- "Proc"

Idents(TrioBrain.integrated_slim_labeled, cells = Dh31I1_cellid) <- "Dh31(I1)"
Idents(TrioBrain.integrated_slim_labeled, cells = Da6I_cellid) <- "Da6(I)"
Idents(TrioBrain.integrated_slim_labeled, cells = fruI_cellid) <- "fru(I)"
Idents(TrioBrain.integrated_slim_labeled, cells = sNPFI_cellid) <- "sNPF(I)"

Idents(TrioBrain.integrated_slim_labeled, cells = Mip1_cellid) <- "Mip(1)"

Idents(TrioBrain.integrated_slim_labeled, cells = DOP_cellid) <- "DOP"
Idents(TrioBrain.integrated_slim_labeled, cells = SER_cellid) <- "SER"
Idents(TrioBrain.integrated_slim_labeled, cells = OCTY_cellid) <- "OCTY"
Idents(TrioBrain.integrated_slim_labeled, cells = TY_cellid) <- "TY"
Idents(TrioBrain.integrated_slim_labeled, cells = Tbh_cellid) <- "Tbh"

Idents(TrioBrain.integrated_slim_labeled, cells = fruGlu_cellid) <- "fru(Glu)"
Idents(TrioBrain.integrated_slim_labeled, cells = fruMs_cellid) <- "fru(Ms)"

Idents(TrioBrain.integrated_slim_labeled, cells = CCHa22_cellid) <- "CCHa2(2)"
Idents(TrioBrain.integrated_slim_labeled, cells = FMRFa2_cellid) <- "FMRFa(2)"

Idents(TrioBrain.integrated_slim_labeled, cells = PRN_cellid) <- "PRN"
Idents(TrioBrain.integrated_slim_labeled, cells = SUB_cellid) <- "SUB"
Idents(TrioBrain.integrated_slim_labeled, cells = Blood_cellid) <- "Blood"

Idents(TrioBrain.integrated_slim_labeled, cells = Fat_cellid) <- "Fat"

Idents(TrioBrain.integrated_slim_labeled, cells = Hug_cellid) <- "Hug"
Idents(TrioBrain.integrated_slim_labeled, cells = IPC_cellid) <- "IPC"
Idents(TrioBrain.integrated_slim_labeled, cells = Tk_cellid) <- "Tk"
Idents(TrioBrain.integrated_slim_labeled, cells = Mip2_cellid) <- "Mip(2)"
Idents(TrioBrain.integrated_slim_labeled, cells = Crz_cellid) <- "Crz"

anno_cluster <- c("GABA", "DOP", "SER", "OCTY", "TY", "Tbh", "MON", "Da5/7","Da6(E)", "Da6(I)","ort",
                  "ENS", "AST", "PRN", "SUB", "CTX",
                  "αβ-KC", "γ-KC", "α'β'-KC",
                  "sNPF(E)", "sNPF(I)", "CCHa2(1)","CCHa2(2)", "Proc", "Mip(1)", "Mip(2)", "FMRFa(1)",  "FMRFa(2)",
                  "AstA", "RYa", "Hug", "IPC", "Tk", "SIFa", "Crz", 
                  "Dh31(E)","Dh31(I1)","Dh31(I2)","Dh44",
                  "OPN", "fru(E)", "fru(Glu)","fru(Ms)","fru(I)","Clock",  "Poxn", "Fat", "Blood")

## cluster size quantification and cell type order by cluster size ##

TrioBrain.integrated_slim_labeled[["CellType"]] <- Idents(object = TrioBrain.integrated_slim_labeled)
#TrioBrain.integrated_slim_labeled <- StashIdent(object = TrioBrain.integrated_slim_labeled, save.name = "CellType")

df_TrioBrain.integrated_slim_labeled_metadata <- TrioBrain.integrated_slim_labeled@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))

#final_cluster_id <- unique(c(new.cluster.ids, anno_cluster))
#unique(df_TrioBrain.integrated_slim_labeled_metadata$CellType)[!unique(df_TrioBrain.integrated_slim_labeled_metadata$CellType) %in% final_cluster_id]
#final_cluster_id %in% unique(df_TrioBrain.integrated_slim_labeled_metadata$CellType)

df_TrioBrain.integrated_slim_labeled_metadata_summary_rep <- df_TrioBrain.integrated_slim_labeled_metadata %>%
  dplyr::group_by(species, CellType, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

df_TrioBrain.integrated_slim_labeled_metadata_summary <- df_TrioBrain.integrated_slim_labeled_metadata_summary_rep %>%
  dplyr::group_by(species, CellType) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

CellType_order <- df_TrioBrain.integrated_slim_labeled_metadata_summary %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(mean_freq=mean(percent_combined)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-mean_freq) %>%
  dplyr::select(CellType, mean_freq) %>%
  as.vector()

CellType_order_anno <- df_TrioBrain.integrated_slim_labeled_metadata_summary %>%
  dplyr::filter(!grepl("Cluster", CellType)) %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(mean_freq=mean(percent_combined)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-mean_freq) %>%
  dplyr::select(CellType, mean_freq) %>%
  as.vector()

#save(anno_cluster, CellType_order, CellType_order_anno, file="Processed_Data/celltype_order.RData")
load("Processed_Data/celltype_order.RData")

df_TrioBrain.integrated_slim_labeled_metadata$CellType <- factor(df_TrioBrain.integrated_slim_labeled_metadata$CellType, levels=CellType_order$CellType,
                                                                 labels = CellType_order$CellType)

df_TrioBrain.integrated_slim_labeled_metadata_summary_rep <- df_TrioBrain.integrated_slim_labeled_metadata %>%
  dplyr::group_by(species, CellType, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

df_TrioBrain.integrated_slim_labeled_metadata_summary <- df_TrioBrain.integrated_slim_labeled_metadata_summary_rep %>%
  dplyr::group_by(species, CellType) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

df_TrioBrain.integrated_slim_labeled_metadata_summary %>%
  ggplot(.) +
  geom_col(aes(y=CellType, x=percent_combined), alpha = 0.8) + 
  geom_errorbar(aes(y=CellType, x=percent_combined,
                    xmin=percent_combined-sem, xmax=percent_combined+sem), width=.2, position=position_dodge(.9)) +
  geom_point(data=df_TrioBrain.integrated_slim_labeled_metadata_summary_rep, size = 2, alpha =0.8, shape=21,  aes(y=CellType, x=percent, fill=species)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=CellType, y=n/total*100)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11, color='black'),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=10, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        panel.grid = element_blank()) +
  labs(y="Cluster", x="Percent of cell type (%)", fill="")+
  facet_grid(~species)

#ggsave("Plots/Trio_combined_full_intron/TrioBrain.integrated_slim_labeled_cluster_size_percent_facet.png", width=10, height = 12)

df_TrioBrain.integrated_slim_labeled_metadata_summary %>%
  ggplot(.) +
  geom_col(aes(x=CellType, y=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.9)) + 
  geom_errorbar(aes(x=CellType, y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,0,percent_combined-sem), ymax=percent_combined+sem, group=species), width=0.9, position=position_dodge(width=0.9)) +
  geom_point(data=df_TrioBrain.integrated_slim_labeled_metadata_summary_rep,
             size = 1, alpha =0.8,  
             aes(x=CellType, y=percent, group=species,fill=species), shape=21, position=position_dodge(width=0.9)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=CellType, y=n/total*100)) +
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

ggsave("Plots/TrioBrain_DF.integrated_cluster_size_percent_alhor.png", width=20, height = 10)

df_TrioBrain.integrated_slim_labeled_metadata_summary %>%
  dplyr::filter(species != "DsecNoni") %>%
  ggplot(.) +
  geom_col(aes(x=species, y=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.8)) + 
  geom_errorbar(aes(x=species, y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,0,percent_combined-sem), ymax=percent_combined+sem, group=species), width=0.8, position=position_dodge(width=0.8)) +
  geom_point(data=dplyr::filter(df_TrioBrain.integrated_slim_labeled_metadata_summary_rep, species != "DsecNoni"),
             size = 1, alpha =0.8,  
             aes(x=species, y=percent, group=species,fill=species), shape=21, position=position_dodge(width=0.8)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=CellType, y=n/total*100)) +
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
  facet_wrap(~CellType, scales='free', ncol=8)

ggsave("Plots/FigS6_celltype_frequencies_all.pdf", width=13, height = 16)

for (cluster in unique(df_TrioBrain.integrated_slim_labeled_metadata_summary$CellType)) {
  
  clustesr_namefix <-  gsub("/", "_", cluster)
  
  df_trio_reps <- df_TrioBrain.integrated_slim_labeled_metadata_summary_rep %>%
    dplyr::filter(CellType == cluster) %>%
    dplyr::group_by(species, orig.ident) %>%
    dplyr::summarise(percent = sum(percent)) %>%
    dplyr::ungroup()
  
  df_trio_summary <- df_trio_reps %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(mean=mean(percent), sem=sd(percent)/sqrt(6)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))
  
  df_trio_summary %>%
    ggplot(.) +
    geom_col(aes(y=species, x=mean, fill=species), alpha = 1, position = "dodge") + 
    geom_point(data=df_trio_reps, size = 2, alpha =0.8,  
               aes(y=species, x=percent, group=species), fill='grey', shape=21) +
    geom_errorbar(aes(y=species, x=mean,
                      xmin=ifelse(mean-sem<0,0,mean-sem), 
                      xmax=mean+sem), 
                  alpha=0.7, width=0.5) +
    #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=CellType, y=n/total*100)) +
    theme_bw() +
    theme(axis.title.x = element_text(size=12, color='black'),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=10, color='black'),
          axis.text.y=element_text(size=10, face='italic', color='black'),
          panel.grid = element_blank(),
          legend.position = 'none') +
    labs(y="Species", x=glue::glue("Percent of {cluster} cells (%)")) +
    scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3"))
  
  ggsave(glue::glue("Plots/comparison/frequencies_SCT_rpca_DF02_48annos/Trio_merged_percent_{clustesr_namefix}.png"), width=3.5, height = 3.5)
  
  
}


DotPlot(TrioBrain.integrated_slim_labeled, features = c("axo", "Ilp2", "Crz", "Mip"))

TrioBrain.integrated_slim_labeled@active.ident <- factor(TrioBrain.integrated_slim_labeled@active.ident, levels=CellType_order$CellType,
                                                                 labels = CellType_order$CellType)

DimPlot(TrioBrain.integrated_slim_labeled, reduction = "umap", label = TRUE, pt.size = 0.01, raster=F) + NoLegend()
DimPlot(TrioBrain.integrated_slim_labeled, reduction = "tsne", label = TRUE, pt.size = 0.01, raster=F) + NoLegend() + NoAxes()
FeaturePlot(subset(TrioBrain.integrated_slim_labeled, orig.ident != "DsecNoni", ident='Clock'), reduction = "tsne", feature=c("Imp","pros"))

IPCcells <- WhichCells(TrioBrain.integrated_slim_labeled, idents = "IPC")
Fatcells <- WhichCells(TrioBrain.integrated_slim_labeled, idents = "Fat")
Crzcells <- WhichCells(TrioBrain.integrated_slim_labeled, idents = "Crz")
FruIcells <- WhichCells(TrioBrain.integrated_slim_labeled, idents = "fru(I)")

DimPlot(TrioBrain.integrated_slim_labeled, reduction = "tsne", 
        cells.highlight = FruIcells, shuffle = T, raster=F, sizes.highlight = 0.01, pt.size = 0.01) + NoLegend() + labs(title=NULL)

DimPlot(subset(TrioBrain.integrated_slim_labeled, orig.ident != "DsecNoni"), 
        reduction = "tsne", label = F, pt.size = 0.01, raster=F) + NoLegend() + NoAxes()

ggsave("Plots/Manuscript/Fig1a_clustering.png", width=11, height = 10)


DimPlot(subset(TrioBrain.integrated_slim_labeled, orig.ident != "DsecNoni"), 
        reduction = "tsne", label = T, pt.size = 0.01, raster=F) + NoLegend() + NoAxes()

ggsave("Plots/Manuscript/Fig1a_clustering_labeled.png", width=11, height = 10)

FeaturePlot(subset(TrioBrain.integrated_slim_labeled, orig.ident != "DsecNoni"), 
        reduction = "tsne", feature="Tret1-1", order=T, raster=F) + NoAxes()

ggsave("Plots/tsne_tret1_1.png", width=11, height = 10)

FeaturePlot(subset(TrioBrain.integrated_slim_labeled, orig.ident != "DsecNoni"), 
            reduction = "tsne", feature="crb", order=T, raster=F) + NoAxes()

ggsave("Plots/tsne_crb.png", width=11, height = 10)

DimPlot(subset(TrioBrain.integrated_slim_labeled, orig.ident != "DsecNoni", idents = c('Clock', 'Fat', 'IPC', 'Poxn')), 
        reduction = "tsne", label = F, pt.size = 0.01, raster=F) + NoLegend() + NoAxes()

#saveRDS(TrioBrain.integrated_slim_labeled, file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled.rds")
TrioBrain.integrated_slim_labeled <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled.rds")

DefaultAssay(TrioBrain.integrated_slim_labeled)

TrioBrain.integrated_slim_labeled.markers <- FindAllMarkers(TrioBrain.integrated_slim_labeled, only.pos = TRUE, min.pct = 0.15, 
                                               logfc.threshold = 0.25, test.use = "MAST")
save(TrioBrain.integrated_slim_labeled.markers, file="Processed_Data/cell_type_markers.RData")

# Calculate the number of transcripts (UMIs) per cell
nUMI <- colSums(TrioBrain.integrated_slim_labeled@assays$RNA@counts)
names(nUMI) <- colnames(TrioBrain.integrated_slim_labeled)
TrioBrain.integrated_slim_labeled <- TrioBrain.integrated_slim_labeled %>% AddMetaData(object = ., metadata = nUMI, col.name = "nUMI")

# Calculate the number of genes detected per cell
nGene <- colSums(TrioBrain.integrated_slim_labeled@assays$RNA@counts > 0)
names(nGene) <- colnames(TrioBrain.integrated_slim_labeled)
TrioBrain.integrated_slim_labeled <- TrioBrain.integrated_slim_labeled %>% AddMetaData(object = ., metadata = nGene, col.name = "nGene")

# Extract metadata containing nUMI and nGene
TrioBrain_metadata <- TrioBrain.integrated_slim_labeled@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))

TrioBrain_metadata_summary <- TrioBrain_metadata %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(mean_nGene=mean(nGene), mean_nUMI=mean(nUMI))


# Display the first few rows of metadata
head(TrioBrain_metadata)



backgroud_DEG <- VariableFeatures(TrioBrain.integrated_slim_labeled)
save(backgroud_DEG, file="Processed_Data/background_vari_features.RData")

## unannotated clusters ##

DotPlot(TrioBrain.integrated_slim_labeled, feature=c("nSyb", "SK")) 
all_clusters <- levels(Idents(TrioBrain.integrated_slim_labeled))
unannotated <- setdiff(all_clusters,anno_cluster)

subset_unanno <- subset(TrioBrain.integrated_slim_labeled, idents = unannotated)

DotPlot(subset(TrioBrain.integrated_slim_labeled, idents = unannotated), feature=c("nSyb", "SK")) 

## Marker for each plot ##

DotPlot(TrioBrain.integrated_slim_labeled, feature=c("axo","Gs2","Tret1-1","CG40470")) ## Glial cell types

DotPlot(TrioBrain.integrated_slim_labeled, feature=c("crb","CG32204","Pka-C1", "mamo")) ## Kenyon cell types

DotPlot(TrioBrain.integrated_slim_labeled, feature=c("Gad1", "Vmat", "nAChRalpha5","nAChRalpha6","nAChRalpha7", "ort",
                                                    "Proc","Mip", "sNPF","CCHa2","AstA","RYa","Dh31", "Dh44", "SIFa")) ## Chemoconnectome (Neurotransmitter, neuropeptides, hormone)

DotPlot(TrioBrain.integrated_slim_labeled, feature=c("Crz", "Dsk","Hug","Ilp2","Ilp3","Ilp5")) 


DotPlot(TrioBrain.integrated_slim_labeled, feature=c("fru", "Oaz", "tim", "Clk", "Poxn")) ## others

DotPlot(TrioBrain.integrated_slim_labeled, feature=c("VAChT","VGlut", "Gad1", "Vmat",
                                                     "axo","Gs2","Tret1-1","CG40470",
                                                     "crb","CG32204","Pka-C1", "mamo",
                                                     "sNPF","CCHa2","AstA","RYa","Dh31", "Dh44", "SIFa",
                                                     "fru", "tim", "Clk", "Poxn"))

TrioBrain.integrated_slim_labeled_filter <- subset(TrioBrain.integrated_slim_labeled, idents = anno_cluster)

TrioBrain.integrated_slim_labeled_filter@active.ident <- factor(TrioBrain.integrated_slim_labeled_filter@active.ident,
                                                                levels= anno_cluster)

DotPlot(TrioBrain.integrated_slim_labeled_filter, 
        dot.scale = 4.5,
        feature=c("pros","Imp", 
                  "VAChT","VGlut", "Gad1", "Vmat", "DAT", "Trh", "Tdc2", "Tbh", 
                  "nAChRalpha5","nAChRalpha6","nAChRalpha7", "ort",
                                                            "axo","Gs2","Tret1-1", "baz", "CG40470",
                                                            "crb","CG32204","Pka-C1", "mamo",
                                                            "sNPF", "CCHa2",  "Proc", "Mip", "FMRFa", "AstA", "RYa", "Hug", "Ilp2","Ilp3","Ilp5","Tk", "SIFa", "Crz", "Ms",
                                                            "Dh31", "Dh44",
                                                            "Oaz", "SiaT", "fru", "tim", "Clk", "Poxn", "fit", "Hml")) +
  coord_flip() + labs(x=NULL, y="Clusters") +
  theme(legend.position='bottom',
        axis.title = element_blank(),
        axis.text.x = element_text(color='black', angle=45, hjust=1, vjust=1))

ggsave("Plots/Manuscript/Fig1_cluster_markers.pdf", width=15, height = 10)

## check each species ##

TrioBrain.integrated_slim_labeled_filter_species <- TrioBrain.integrated_slim_labeled_filter

TrioBrain.integrated_slim_labeled_filter_species@meta.data$orig.ident <- gsub("_rep1", "", TrioBrain.integrated_slim_labeled_filter_species@meta.data$orig.ident)
TrioBrain.integrated_slim_labeled_filter_species@meta.data$orig.ident <- gsub("_rep2", "", TrioBrain.integrated_slim_labeled_filter_species@meta.data$orig.ident)
TrioBrain.integrated_slim_labeled_filter_species@meta.data$orig.ident <- gsub("_rep3", "", TrioBrain.integrated_slim_labeled_filter_species@meta.data$orig.ident)
TrioBrain.integrated_slim_labeled_filter_species@meta.data$orig.ident <- gsub("_rep4", "", TrioBrain.integrated_slim_labeled_filter_species@meta.data$orig.ident)
TrioBrain.integrated_slim_labeled_filter_species@meta.data$orig.ident <- gsub("_rep5", "", TrioBrain.integrated_slim_labeled_filter_species@meta.data$orig.ident)
TrioBrain.integrated_slim_labeled_filter_species@meta.data$orig.ident <- gsub("_rep6", "", TrioBrain.integrated_slim_labeled_filter_species@meta.data$orig.ident)

TrioBrain.integrated_slim_labeled_filter_species@meta.data$orig.ident <- factor(TrioBrain.integrated_slim_labeled_filter_species@meta.data$orig.ident, 
                                                                         levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))

# Dmel #

DotPlot(subset(TrioBrain.integrated_slim_labeled_filter_species, orig.ident=="Dmel"), 
        dot.scale = 4.5,
        feature=c("pros","Imp", 
                  "VAChT","VGlut", "Gad1", "Vmat", "DAT", "Trh", "Tdc2", "Tbh", 
                  "nAChRalpha5","nAChRalpha6","nAChRalpha7", "ort",
                  "axo","Gs2","Tret1-1", "baz", "CG40470",
                  "crb","CG32204","Pka-C1", "mamo",
                  "sNPF", "CCHa2",  "Proc", "Mip", "FMRFa", "AstA", "RYa", "Hug", "Ilp2","Ilp3","Ilp5","Tk", "SIFa", "Crz", "Ms",
                  "Dh31", "Dh44",
                  "Oaz", "SiaT", "fru", "tim", "Clk", "Poxn", "fit", "Hml")) +
  coord_flip() + labs(x=NULL, y="Clusters") +
  theme(legend.position='bottom',
        axis.title = element_blank(),
        axis.text.x = element_text(color='black', angle=90, hjust=1, vjust=0.5))

ggsave("Plots/Manuscript/FigS2_cluster_markers_Dmel.pdf", width=15, height = 10)


# Dsim #

DotPlot(subset(TrioBrain.integrated_slim_labeled_filter_species, orig.ident=="Dsim"), 
        dot.scale = 4.5,
        feature=c("pros","Imp", 
                  "VAChT","VGlut", "Gad1", "Vmat", "DAT", "Trh", "Tdc2", "Tbh", 
                  "nAChRalpha5","nAChRalpha6","nAChRalpha7", "ort",
                  "axo","Gs2","Tret1-1", "baz", "CG40470",
                  "crb","CG32204","Pka-C1", "mamo",
                  "sNPF", "CCHa2",  "Proc", "Mip", "FMRFa", "AstA", "RYa", "Hug", "Ilp2","Ilp3","Ilp5","Tk", "SIFa", "Crz", "Ms",
                  "Dh31", "Dh44",
                  "Oaz", "SiaT", "fru", "tim", "Clk", "Poxn", "fit", "Hml")) +
  coord_flip() + labs(x=NULL, y="Clusters") +
  theme(legend.position='bottom',
        axis.title = element_blank(),
        axis.text.x = element_text(color='black', angle=90, hjust=1, vjust=0.5))

ggsave("Plots/Manuscript/FigS2_cluster_markers_Dsim.pdf", width=15, height = 10)

# Dsec #

DotPlot(subset(TrioBrain.integrated_slim_labeled_filter_species, orig.ident=="Dsec"), 
        dot.scale = 4.5,
        feature=c("pros","Imp", 
                  "VAChT","VGlut", "Gad1", "Vmat", "DAT", "Trh", "Tdc2", "Tbh", 
                  "nAChRalpha5","nAChRalpha6","nAChRalpha7", "ort",
                  "axo","Gs2","Tret1-1", "baz", "CG40470",
                  "crb","CG32204","Pka-C1", "mamo",
                  "sNPF", "CCHa2",  "Proc", "Mip", "FMRFa", "AstA", "RYa", "Hug", "Ilp2","Ilp3","Ilp5","Tk", "SIFa", "Crz", "Ms",
                  "Dh31", "Dh44",
                  "Oaz", "SiaT", "fru", "tim", "Clk", "Poxn", "fit", "Hml")) +
  coord_flip() + labs(x=NULL, y="Clusters") +
  theme(legend.position='bottom',
        axis.title = element_blank(),
        axis.text.x = element_text(color='black', angle=90, hjust=1, vjust=0.5))

ggsave("Plots/Manuscript/FigS2_cluster_markers_Dsec.pdf", width=15, height = 10)


## Correlation

df_cell_comp <- df_TrioBrain.integrated_slim_labeled_metadata_summary %>%
  dplyr::select(-sem) %>%
  tidyr::spread(key="species", value="percent_combined") %>%
  dplyr::mutate(Dsec_Dmel=Dsec/Dmel, Dsec_Dsim=Dsec/Dsim, Dsim_Dmel = Dsim/Dmel) 

df_cell_comp %>%
  ggplot(.) +
  geom_point() +
  geom_smooth() +
  geom_label_repel(data=dplyr::filter(df_cell_comp, (abs(log2(Dsec_Dmel)) > 0.5 & abs(log2(Dsec_Dsim)) > 0.5)), 
                   force=5,
                   aes(label=CellType)) +
  aes(x=log2(Dsec_Dmel), y=log2(Dsec_Dsim)) +
  theme_bw()+
  lims(x=c(-1,1),y=c(-1,1))

ggsave("Plots/TrioBrain_cell_comp_rho.png", width = 10, height = 8)

cor(df_cell_comp$Dsec_Dmel, df_cell_comp$Dsec_Dsim, method = 'spearman')

### ANOVA ###

df_TrioBrain_anova <- df_TrioBrain.integrated_slim_labeled_metadata_summary_rep %>%
  dplyr::filter(species != "DsecNoni") %>%
  dplyr::select(CellType, species, percent) %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(p_value=tidy(aov(percent ~ species))$p.value[1]) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_value_adjust = p.adjust(p_value, method = "fdr")) %>%
  dplyr::mutate(sig = ifelse(p_value_adjust < 0.05, T, F))

df_TrioBrain_stat_split <- df_TrioBrain.integrated_slim_labeled_metadata_summary_rep %>%
  dplyr::filter(species != "DsecNoni") %>%
  dplyr::select(CellType, species, percent) %>%
  group_split(CellType)

clusters <- unique(df_TrioBrain.integrated_slim_labeled_metadata_summary_rep$CellType)

df_comp_padj <- NULL

for (i in c(1:length(unique(clusters)))) {
  
  cl=as.character(clusters[i])
  
  padj <- t(TukeyHSD(aov(data=df_TrioBrain_stat_split[[i]], percent ~ species), conf.level=.95)$species)[4,]
  df_padj <- data.frame(cluster=cl, 
                             Dsim_Dmel=as.numeric(padj[1]), 
                             Dsec_Dmel=as.numeric(padj[2]), 
                             Dsec_Dsim=as.numeric(padj[3]))
  
  df_comp_padj <- rbind(df_comp_padj, df_padj)
  
}

df_comp_padj_gather <- df_comp_padj %>%
  tidyr::gather(-cluster, key="pair", value="padj")

df_cell_comp_padj <- df_cell_comp %>%
  dplyr::select(cluster=CellType, Dsim_Dmel, Dsec_Dmel, Dsec_Dsim) %>%
  tidyr::gather(-cluster, key="pair", value="fold_diff") %>%
  dplyr::left_join(., df_comp_padj_gather, by=c('cluster', 'pair')) %>%
  dplyr::mutate(log2fc = log2(fold_diff)) %>%
  dplyr::group_by(pair) %>%
  dplyr::mutate(p_value_adjust = p.adjust(padj, method = "fdr")) %>%
  dplyr::mutate(sig = ifelse(p_value_adjust < 0.05, T, F)) %>%
  dplyr::ungroup()

df_cell_comp_padj %>%
  #dplyr::filter(pair != "Dsim_Dmel") %>%
  ggplot(.) +
  geom_vline(xintercept=0, size=0.2, alpha=0.5) +
  geom_hline(yintercept=-log10(0.05/nrow(CellType_order)), size=0.2, alpha=0.5, linetype=2) +
  geom_point(size=2, alpha=0.8) +
  geom_text_repel(data=dplyr::filter(df_cell_comp_padj, cluster %in% c("PRN", "sNPF(E)", "Dh44")), 
                   force=5, nudge_y=0.2, nudge_x=0.2,
                   aes(label=cluster)) +
  aes(x=log2fc, y=-log10(padj), color=pair) +
  theme_bw() +
  scale_color_manual(values=c("orange", "navy",  "brown")) +
  labs(x="log2(fold difference)", y="-log10(adjusted P-value)") +
  theme(panel.grid=element_blank(),
        axis.title = element_text(size=12, color='black'),
        axis.text = element_text(size=11, color='black'),
        legend.title=element_blank(),
        legend.text = element_text(size=11, color='black')) +
  scale_y_continuous(limit=c(-0.05,5), expand=c(0,0))

ggsave("Plots/Manuscript/Fig2_volcano.pdf", width=7.5, height = 5)

df_cell_comp_padj %>%
  #dplyr::filter(pair != "Dsim_Dmel") %>%
  ggplot(.) +
  geom_vline(xintercept=0, size=0.2, alpha=0.5) +
  geom_hline(yintercept=-log10(0.05/64), size=0.2, alpha=0.5, linetype=2) +
  geom_point(size=2, alpha=0.8) +
  geom_text_repel(data=dplyr::filter(df_cell_comp_padj, padj<0.01), 
                  force=5,
                  aes(label=cluster)) +
  aes(x=log2fc, y=-log10(padj), color=pair) +
  theme_bw() +
  scale_color_manual(values= c("darkblue", "brown", "orange")) +
  labs(x="log2(fold difference)", y="-log10(adjusted P-value)") +
  theme(panel.grid=element_blank(),
        axis.title = element_text(size=12, color='black'),
        axis.text = element_text(size=11, color='black'),
        legend.title=element_blank(),
        legend.text = element_text(size=11, color='black')) +
  scale_y_continuous(limit=c(-0.05,5), expand=c(0,0))

df_cell_comp_padj %>%
  ggplot(.) +
  geom_boxplot() +
  aes(x=pair, y=-log10(padj)) +
  theme_bw()

df_cell_comp_padj %>%
  ggplot(.) +
  geom_boxplot() +
  aes(x=pair, y=abs(log2fc)) +
  theme_bw()

View(as_tibble(list_padj))

df_TrioBrain_anova_posthoc <- df_TrioBrain.integrated_slim_labeled_metadata_summary_rep %>%
  dplyr::filter(species != "DsecNoni") %>%
  dplyr::select(CellType, species, percent) %>%
  dplyr::group_by(CellType) %>%
  do(multitst = tidy(TukeyHSD(aov(percent ~ species, data = .))))

clu_sel <- c("α'β'-KC", "fru(E)", "Clock", "sNPF(I)", "PRN", "sNPF(E)", "Dh44")

df_TrioBrain.integrated_slim_labeled_metadata_summary %>%
  dplyr::filter(species != "DsecNoni", CellType %in% clu_sel) %>%
  ggplot(.) +
  geom_col(aes(x=factor(CellType, levels=clu_sel), y=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.9)) + 
  geom_errorbar(aes(x=factor(CellType, levels=clu_sel), y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,0,percent_combined-sem), ymax=percent_combined+sem, group=species), width=0.9, position=position_dodge(width=0.9)) +
  geom_point(data=dplyr::filter(df_TrioBrain.integrated_slim_labeled_metadata_summary_rep,
                                species != "DsecNoni", CellType %in% clu_sel),
             size = 1.8, alpha =0.8,  
             aes(x=CellType, y=percent, group=species,fill=species), shape=21, position=position_dodge(width=0.9)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=CellType, y=n/total*100)) +
  theme_bw() +
  theme(axis.title.y = element_text(size=13, color='black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        legend.text = element_text(size=13, face='italic', color='black'),
        panel.grid = element_blank(),
        legend.position=c(0.9,0.85),
        legend.background = element_blank()) +
  labs(x="Cluster", y="Frequency (% of cells)", fill="")  +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) +
  scale_y_continuous(limits = c(-0.05, 2.5), expand=c(0,0))

ggsave("Plots/Manuscript/Fig2_barplots.pdf", width=6, height = 4.2)


### Cell frequency divergence ###

TrioBrain.integrated_slim_labeled_species <- TrioBrain.integrated_slim_labeled

TrioBrain.integrated_slim_labeled_species@meta.data$orig.ident <- gsub("_rep1", "", TrioBrain.integrated_slim_labeled_species@meta.data$orig.ident)
TrioBrain.integrated_slim_labeled_species@meta.data$orig.ident <- gsub("_rep2", "", TrioBrain.integrated_slim_labeled_species@meta.data$orig.ident)
TrioBrain.integrated_slim_labeled_species@meta.data$orig.ident <- gsub("_rep3", "", TrioBrain.integrated_slim_labeled_species@meta.data$orig.ident)
TrioBrain.integrated_slim_labeled_species@meta.data$orig.ident <- gsub("_rep4", "", TrioBrain.integrated_slim_labeled_species@meta.data$orig.ident)
TrioBrain.integrated_slim_labeled_species@meta.data$orig.ident <- gsub("_rep5", "", TrioBrain.integrated_slim_labeled_species@meta.data$orig.ident)
TrioBrain.integrated_slim_labeled_species@meta.data$orig.ident <- gsub("_rep6", "", TrioBrain.integrated_slim_labeled_species@meta.data$orig.ident)

TrioBrain.integrated_slim_labeled_species@meta.data$orig.ident <- factor(TrioBrain.integrated_slim_labeled_species@meta.data$orig.ident, 
                                                                                levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))

#saveRDS(TrioBrain.integrated_slim_labeled_species, file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled_species.rds")
TrioBrain.integrated_slim_labeled_species <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled_species.rds")


DotPlot(TrioBrain.integrated_slim_labeled_species, feature=c("CAH2","CAH3"))
DotPlot(TrioBrain.integrated_slim_labeled_species, feature=c("CG9394","CG18135"))

DimPlot(TrioBrain.integrated_slim_labeled_species, reduction = "umap", 
              group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01, 
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL)

sNPFcells <- WhichCells(TrioBrain.integrated_slim_labeled_species, idents = "sNPF(E)")
Dh44cells <- WhichCells(TrioBrain.integrated_slim_labeled_species, idents = "Dh44")
PRNcells <- WhichCells(TrioBrain.integrated_slim_labeled_species, idents = "PRN")

DimPlot(TrioBrain.integrated_slim_labeled_species, reduction = "tsne", 
        cells.highlight = sNPFcells, shuffle = T, raster=F, pt.size = 0.01) + NoLegend() + labs(title=NULL)
DimPlot(TrioBrain.integrated_slim_labeled_species, reduction = "tsne", 
        cells.highlight = Dh44cells, shuffle = T, raster=F, pt.size = 0.01) + NoLegend() + labs(title=NULL)

DimPlot(TrioBrain.integrated_slim_labeled_species, label=F,shuffle=T,
        cells.highlight= list(Dh44=Dh44cells, PRN=PRNcells), 
        reduction = "tsne", 
        pt.size = 0.01,sizes.highlight = 0.01,
        cols.highlight = c("darkorange","blue"), cols= "grey40", raster=FALSE)

ggsave("Plots/Manuscript/Fig2_highlight.png", width=10, height = 8)

FeaturePlot(TrioBrain.integrated_slim_labeled_species, feature=c("sNPF", "CCHa2"), reduction = "tsne")

### PRN ###

DimPlot(subset(TrioBrain.integrated_slim_labeled_species, ident = "PRN", orig.ident != "DsecNoni"), reduction = "tsne", 
        group.by = "orig.ident", shuffle = T, raster=F, pt.size = 1, 
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL) +
  coord_cartesian(xlim=c(2.5,5.5), ylim=c(-9.5,-12))+ NoAxes()

ggsave("Plots/Manuscript/Fig2_PRN.pdf", width=8, height = 8)

DimPlot(subset(TrioBrain.integrated_slim_labeled_species, ident = "PRN", orig.ident != "DsecNoni"), reduction = "tsne", 
        group.by = "orig.ident", split.by= "orig.ident", shuffle = T, raster=F, pt.size = 1, 
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL) +
  coord_cartesian(xlim=c(2.5,5.5), ylim=c(-9.5,-12)) + NoAxes()

ggsave("Plots/Manuscript/Fig2_PRN_species.pdf", width=16, height = 6)

### Dh44 ###

DimPlot(subset(TrioBrain.integrated_slim_labeled_species, ident = "Dh44", orig.ident != "DsecNoni"), reduction = "tsne", 
        group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.5, 
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL) +
  coord_cartesian(xlim=c(-11.2,-10), ylim=c(13,13.8))+ NoAxes()

ggsave("Plots/Manuscript/Fig2_Dh44.pdf", width=8, height = 8)

DimPlot(subset(TrioBrain.integrated_slim_labeled_species, ident = "Dh44", orig.ident != "DsecNoni"), reduction = "tsne", 
        group.by = "orig.ident", split.by= "orig.ident", shuffle = T, raster=F, pt.size = 0.5, 
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL) +
  coord_cartesian(xlim=c(-11.2,-10), ylim=c(13,13.8)) + NoAxes()

ggsave("Plots/Manuscript/Fig2_Dh44_species.pdf", width=16, height = 6)

### sNPF ###

DimPlot(subset(TrioBrain.integrated_slim_labeled_species, ident = c("sNPF(E)"), orig.ident != "DsecNoni"), reduction = "tsne", 
        group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.8, 
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL) +
  coord_cartesian(xlim=c(-14,-9), ylim=c(8,12)) + NoAxes()

ggsave("Plots/Manuscript/Fig2_sNPF.pdf", width=8, height = 8)

DimPlot(subset(TrioBrain.integrated_slim_labeled_species, ident = "sNPF(E)", orig.ident != "DsecNoni"), reduction = "tsne", 
        group.by = "orig.ident", split.by= "orig.ident", shuffle = T, raster=F, pt.size = 0.8, 
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL)+
  coord_cartesian(xlim=c(-14,-9), ylim=c(8,12)) + NoAxes()

ggsave("Plots/Manuscript/Fig2_sNPF_species.pdf", width=16, height = 6)


### Cluster 33 - GO enrichment ###

cluster33_markers <- TrioBrain.integrated.markers %>%
  dplyr::filter(cluster==33) %>%
  dplyr::select(gene)

cluster33_markers_GO <- enrichr(genes=cluster33_markers$gene, database='GO_Biological_Process_AutoRIF')

View(cluster33_markers_GO[[1]])

plotEnrich(cluster33_markers_GO[[1]])

background_markers_GO <- enrichr(genes=unique(TrioBrain.integrated.markers$gene), database='GO_Biological_Process_AutoRIF')
plotEnrich(background_markers_GO[[1]])

View(background_markers_GO[[1]])


