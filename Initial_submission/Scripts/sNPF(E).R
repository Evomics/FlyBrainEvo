library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(sctransform)
library(cowplot)
library(broom)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/shared_genes.RData")

TrioBrain.integrated_slim_labeled<- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled.rds")
TrioBrain.integrated.markers <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated.markers.rds")

df_TrioBrain.integrated_slim_labeled_metadata <- TrioBrain.integrated_slim_labeled@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))

df_TrioBrain.integrated_slim_labeled_metadata_summary_rep <- df_TrioBrain.integrated_slim_labeled_metadata %>%
  dplyr::group_by(species, CellType, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

load("Processed_Data/df_TrioBrain.integrated_slim_labeled_metadata_summary_rep.RData")

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

sNPF_E <- Subclustering(Seurat_object=TrioBrain.integrated_slim_labeled, cluster = 'sNPF(E)', resolution = 0.5, npc=20)

#load("Processed_Data/sNPF_E_subcluster.RData")

sNPF_E_seurat <- sNPF_E[[1]]

sNPF_E_seurat_species <- sNPF_E_seurat

sNPF_E_seurat_species@meta.data$orig.ident <- gsub("_rep1", "", sNPF_E_seurat_species@meta.data$orig.ident)
sNPF_E_seurat_species@meta.data$orig.ident <- gsub("_rep2", "", sNPF_E_seurat_species@meta.data$orig.ident)
sNPF_E_seurat_species@meta.data$orig.ident <- gsub("_rep3", "", sNPF_E_seurat_species@meta.data$orig.ident)
sNPF_E_seurat_species@meta.data$orig.ident <- gsub("_rep4", "", sNPF_E_seurat_species@meta.data$orig.ident)
sNPF_E_seurat_species@meta.data$orig.ident <- gsub("_rep5", "", sNPF_E_seurat_species@meta.data$orig.ident)
sNPF_E_seurat_species@meta.data$orig.ident <- gsub("_rep6", "", sNPF_E_seurat_species@meta.data$orig.ident)

sNPF_E_seurat_species <- subset(sNPF_E_seurat_species, orig.ident != "DsecNoni")

sNPF_E_seurat_species@meta.data$orig.ident <- factor(sNPF_E_seurat_species@meta.data$orig.ident, 
                                                     levels=c("Dmel", "Dsim", "Dsec"))

### Plots for MS ###

Plot_sNPFe_cluster <- DimPlot(sNPF_E_seurat, reduction = 'tsne', label=T, label.box=T, repel=T, pt.size=0.3) + NoAxes() + NoLegend()
DefaultAssay(sNPF_E_seurat) <- 'SCT'
Plot_sNPFe_sNPF <- FeaturePlot(sNPF_E_seurat, reduction = 'tsne', feature='sNPF', pt.size=0.3) + NoAxes() + labs(title=NULL) + NoLegend()

Plot_sNPFe_species <-DimPlot(sNPF_E_seurat_species, reduction = "tsne", 
        group.by = "orig.ident", shuffle = T, raster=FALSE, cols=c("#E41A1C", "#4DAF4A",  "#377EB8"), pt.size=0.3) + labs(title=NULL) + NoAxes() + NoLegend()

Plot_sNPFe_Dmel <-DimPlot(subset(sNPF_E_seurat_species, orig.ident == 'Dmel') , reduction = "tsne", 
                             group.by = "orig.ident", shuffle = T, raster=FALSE, cols="#E41A1C", pt.size=0.5) + labs(title=NULL) + NoAxes() + NoLegend()

Plot_sNPFe_Dsim <-DimPlot(subset(sNPF_E_seurat_species, orig.ident == 'Dsim') , reduction = "tsne", 
                          group.by = "orig.ident", shuffle = T, raster=FALSE, cols="#4DAF4A", pt.size=0.5) + labs(title=NULL) + NoAxes() + NoLegend()

Plot_sNPFe_Dsec <-DimPlot(subset(sNPF_E_seurat_species, orig.ident == 'Dsec') , reduction = "tsne", 
                          group.by = "orig.ident", shuffle = T, raster=FALSE, cols="#377EB8", pt.size=0.5) + labs(title=NULL) + NoAxes() + NoLegend()

plot_sNPFe1 <- cowplot::plot_grid(Plot_sNPFe_cluster, Plot_sNPFe_species, Plot_sNPFe_Dmel, 
                                  Plot_sNPFe_sNPF, Plot_sNPFe_Dsim, Plot_sNPFe_Dsec, nrow=2)
plot_sNPFe1

ggsave(plot_sNPFe1, file="Plots/sNPFe1.pdf", width=14, height=9)

sNPF_E_metadata <- sNPF_E[[5]]

df_exp_rep_size <- df_TrioBrain.integrated_slim_labeled_metadata_summary_rep %>%
  dplyr::distinct(orig.ident, total)

sNPF_E_metadata_total_rep <- sNPF_E_metadata %>%
  dplyr::filter(species != "DsecNoni", seurat_clusters %in% 0:9) %>%
  dplyr::left_join(., df_exp_rep_size, by='orig.ident') %>%
  dplyr::mutate(pct_total = n/total.y*100) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec")))

sNPF_E_metadata_total <- sNPF_E_metadata_total_rep %>%
  dplyr::group_by(species, seurat_clusters) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total.y)*100, sem=sd(pct_total)/sqrt(6)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec")))

plot_sNPFe_freq_bars <- sNPF_E_metadata_total %>%
  ggplot(.) +
  geom_col(aes(y=as.factor(seurat_clusters), x=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.9)) + 
  geom_errorbar(aes(y=as.factor(seurat_clusters), x=percent_combined,
                    xmin=ifelse(percent_combined-sem<0,0,percent_combined-sem), xmax=percent_combined+sem, group=species), width=0.9, position=position_dodge(width=0.9),linewidth=0.2) +
  geom_point(data=sNPF_E_metadata_total_rep, size = 1.3, alpha =0.8,  
             aes(y=as.factor(seurat_clusters), x=pct_total, group=species, fill=species), shape=21, position=position_dodge(width=0.9), stroke = 0.2) +
   theme_bw() +
  theme(axis.title = element_text(size=14, color='black'),
        axis.text.x = element_text(size=13, color='black'),
        axis.text.y=element_text(size=13, color='black'),
        panel.grid = element_blank(),
        legend.position = c(0.95,0.85),
        legend.text=element_text(size=13, color='black', face='italic')) +
  labs(y="Subcluster ID", x="Frequency (% of cells)", fill="") +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8")) +
  scale_x_continuous(expand=c(0,0.003))

plot_sNPFe_freq_bars

ggsave(plot_sNPFe_freq_bars, file="Plots/sNPFe_freq.pdf", width=14, height=5.3)

## DEG cluster marker

sNPF_E_seurat <- PrepSCTFindMarkers(sNPF_E_seurat, assay = "SCT", verbose = TRUE)

sNPF_E_markers <- FindAllMarkers(sNPF_E_seurat, only.pos = F, min.pct = 0.10, 
                                    logfc.threshold = 0.25)

save(sNPF_E_markers, file="Processed_Data/sNPF_E_markers.RData")
load("Processed_Data/sNPF_E_markers.RData")

top3_sNPFe <- sNPF_E_markers %>% 
  dplyr::filter(cluster %in% 0:9) %>%
  group_by(cluster) %>% top_n(n =4, wt = -p_val_adj)

top3_sNPFe_markers <- unique(top3_sNPFe$gene)

plot_sNPFe_subtype_markers <- DotPlot(subset(sNPF_E_seurat, idents=0:9), feature = top3_sNPFe_markers)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=13, color='black'),
        axis.text.y = element_text(color='black', size=12),
        axis.text.x = element_text(color='black', size=12, angle=90, hjust=1, vjust=0.5),
        legend.title=element_text(size=13, color='black'),
        legend.text=element_text(size=12, color='black'))

plot_sNPFe_subtype_markers

plot_sNPFe2 <- cowplot::plot_grid(plot_sNPFe_freq_bars, plot_sNPFe_subtype_markers, nrow=2, rel_heights = c(4,3))
plot_sNPFe2

ggsave(plot_sNPFe2, file="Plots/Manuscript/raw_plots/sNPF/sNPFe2.pdf", width=14, height=9)


DotPlot(sNPF_E_seurat, feature = "Scp2")

plot_sNPFe <- cowplot::plot_grid(plot_sNPFe1,plot_sNPFe2, nrow=2, rel_heights = c(1,1))

ggsave(plot_sNPFe, file="Plots/Manuscript/raw_plots/sNPF/sNPFe.pdf", width=14, height=18)

#save(sNPF_E, df_exp_rep_size, df_TrioBrain.integrated_slim_labeled_metadata_summary_rep, file="Processed_Data/sNPF_E_subcluster.RData")

## statistical significance ##

df_sNPFE_anova <- sNPF_E_metadata_total_rep %>%
  dplyr::filter(species != "DsecNoni") %>%
  dplyr::select(CellType=seurat_clusters, species, percent) %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(p_value=tidy(aov(percent ~ species))$p.value[1]) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_value_adjust = p.adjust(p_value, method = "fdr")) %>%
  dplyr::mutate(sig = ifelse(p_value_adjust < 0.05, T, F))

df_sNPFE_stat_split <- sNPF_E_metadata_total_rep %>%
  dplyr::filter(species != "DsecNoni") %>%
  dplyr::select(CellType=seurat_clusters, species, percent) %>%
  group_split(CellType)

clusters <- unique(sNPF_E_metadata_total_rep$seurat_clusters)

df_comp_padj_sNPFE <- NULL

for (i in c(1:length(unique(clusters)))) {
  
  cl=as.character(clusters[i])
  
  padj <- t(TukeyHSD(aov(data=df_sNPFE_stat_split[[i]], percent ~ species), conf.level=.95)$species)[4,]
  df_padj_sNPFE <- data.frame(cluster=cl, 
                        Dsim_Dmel=as.numeric(padj[1]), 
                        Dsec_Dmel=as.numeric(padj[2]), 
                        Dsec_Dsim=as.numeric(padj[3]))
  
  df_comp_padj_sNPFE <- rbind(df_comp_padj_sNPFE, df_padj_sNPFE)
  
}

df_comp_padj_sNPFE_gather <- df_comp_padj_sNPFE %>%
  tidyr::gather(-cluster, key="pair", value="padj")

df_cell_comp_sNPFE <- sNPF_E_metadata_total %>%
  dplyr::select(-sem) %>%
  tidyr::spread(key="species", value="percent_combined") %>%
  dplyr::mutate(Dsec_Dmel=Dsec/Dmel, Dsec_Dsim=Dsec/Dsim, Dsim_Dmel = Dsim/Dmel) 

df_cell_comp_padj_sNPFE <- df_cell_comp_sNPFE %>%
  dplyr::select(cluster=seurat_clusters, Dsim_Dmel, Dsec_Dmel, Dsec_Dsim) %>%
  tidyr::gather(-cluster, key="pair", value="fold_diff") %>%
  dplyr::left_join(., df_comp_padj_sNPFE_gather, by=c('cluster', 'pair')) %>%
  dplyr::mutate(log2fc = log2(fold_diff)) %>%
  dplyr::group_by(pair) %>%
  dplyr::mutate(p_value_adjust = p.adjust(padj, method = "fdr")) %>%
  dplyr::mutate(sig = ifelse(p_value_adjust < 0.05, T, F)) %>%
  dplyr::ungroup()
