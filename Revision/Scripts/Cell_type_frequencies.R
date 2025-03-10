library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(broom)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

TrioBrain.integrated_slim_labeled_final <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_labeled_final.rds")

TrioBrain.integrated_slim_labeled_final[["ClusterLabel"]] <- Idents(object = TrioBrain.integrated_slim_labeled_final)

df_TrioBrain.integrated_slim_labeled_final_metadata <- TrioBrain.integrated_slim_labeled_final@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))

df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep <- df_TrioBrain.integrated_slim_labeled_final_metadata %>%
  dplyr::group_by(species, ClusterLabel, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

df_TrioBrain.integrated_slim_labeled_final_metadata_summary <- df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep %>%
  dplyr::group_by(species, ClusterLabel) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

ClusterLabel_order <- df_TrioBrain.integrated_slim_labeled_final_metadata_summary %>%
  dplyr::group_by(ClusterLabel) %>%
  dplyr::summarise(mean_freq=mean(percent_combined)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-mean_freq) %>%
  dplyr::select(ClusterLabel, mean_freq) %>%
  as.vector()

ClusterLabel_order_anno <- df_TrioBrain.integrated_slim_labeled_final_metadata_summary %>%
  dplyr::filter(!grepl("cluster", ClusterLabel)) %>%
  #dplyr::filter(!ClusterLabel %in%c(9,19,'Doublets')) %>%
  dplyr::group_by(ClusterLabel) %>%
  dplyr::summarise(mean_freq=mean(percent_combined)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-mean_freq) %>%
  dplyr::select(ClusterLabel, mean_freq) %>%
  as.vector()

Final_idents <- as.character(unique(Idents(TrioBrain.integrated_slim_labeled_final)))
exclude <- c(Final_idents[grep('Unanno', Final_idents)], 0:100)

Anno_idents <- setdiff(Final_idents,exclude)

save(Anno_idents, ClusterLabel_order, ClusterLabel_order_anno, file="Processed_Data/celltype_order.RData")
load("Processed_Data/celltype_order.RData")

df_TrioBrain.integrated_slim_labeled_final_metadata$ClusterLabel <- factor(df_TrioBrain.integrated_slim_labeled_final_metadata$ClusterLabel, levels=ClusterLabel_order$ClusterLabel,
                                                                             labels = ClusterLabel_order$ClusterLabel)

df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep <- df_TrioBrain.integrated_slim_labeled_final_metadata %>%
  dplyr::group_by(species, ClusterLabel, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

df_TrioBrain.integrated_slim_labeled_final_metadata_summary <- df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep %>%
  dplyr::group_by(species, ClusterLabel) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

df_TrioBrain.integrated_slim_labeled_final_metadata_summary %>%
  ggplot(.) +
  geom_col(aes(y=ClusterLabel, x=percent_combined), alpha = 0.8) + 
  geom_errorbar(aes(y=ClusterLabel, x=percent_combined,
                    xmin=percent_combined-sem, xmax=percent_combined+sem), width=.2, position=position_dodge(.9)) +
  geom_point(data=df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep, size = 2, alpha =0.8, shape=21,  aes(y=ClusterLabel, x=percent, fill=species)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=ClusterLabel, y=n/total*100)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11, color='black'),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=10, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        panel.grid = element_blank()) +
  labs(y="Cluster", x="Percent of cell type (%)", fill="")+
  facet_grid(~species)

#ggsave("Plots/Trio_combined_full_intron/TrioBrain.integrated_slim_labeled_final_cluster_size_percent_facet.png", width=10, height = 12)

df_TrioBrain.integrated_slim_labeled_final_metadata_summary %>%
  dplyr::filter(species != "DsecNoni",ClusterLabel %in% Anno_idents) %>%
  ggplot(.) +
  geom_col(aes(x=ClusterLabel, y=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.9)) + 
  geom_errorbar(aes(x=ClusterLabel, y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,0,percent_combined-sem), ymax=percent_combined+sem, group=species), width=0.9, position=position_dodge(width=0.9)) +
  geom_point(data=dplyr::filter(df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep,species != "DsecNoni",ClusterLabel %in% Anno_idents),
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

ggsave("Plots/TrioBrain.final_cluster_size_percent_alhor.png", width=25, height = 10)

df_TrioBrain.integrated_slim_labeled_final_metadata_summary %>%
  dplyr::filter(species != "DsecNoni",ClusterLabel %in% Anno_idents) %>%
  ggplot(.) +
  geom_col(aes(x=species, y=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.8)) + 
  geom_errorbar(aes(x=species, y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,0,percent_combined-sem), ymax=percent_combined+sem, group=species), width=0.8, position=position_dodge(width=0.8)) +
  geom_point(data=dplyr::filter(df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep, species != "DsecNoni",ClusterLabel %in% Anno_idents),
             size = 1, alpha =0.8,  
             aes(x=species, y=percent, group=species,fill=species), shape=21, position=position_dodge(width=0.8)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=ClusterLabel, y=n/total*100)) +
  theme_bw() +
  theme(axis.title.y = element_text(size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=9, color='black'),
        axis.text.x=element_text(size=9, color='black', angle=45, hjust=1, vjust=1, face='italic'),
        legend.text = element_text(size=10, face='italic', color='black'),
        panel.grid = element_blank(),
        strip.text = element_text(size=8.5, color='black'),
        legend.position='none') +
  labs(x="Cluster", y="Frequency (% of cells)", fill="")  +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) +
  facet_wrap(~ClusterLabel, scales='free', ncol=9)

ggsave("Plots/Manuscript/FigS7_celltype_frequencies_all.pdf", width=15, height = 20)

for (cluster in unique(df_TrioBrain.integrated_slim_labeled_final_metadata_summary$ClusterLabel)) {
  
  clustesr_namefix <-  gsub("/", "_", cluster)
  
  df_trio_reps <- df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep %>%
    dplyr::filter(ClusterLabel == cluster) %>%
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
  
  ggsave(glue::glue("Plots/comparison/frequencies_final/Trio_merged_percent_{clustesr_namefix}.png"), width=3.5, height = 3.5)
  
  
}

## Correlation

df_cell_comp <- df_TrioBrain.integrated_slim_labeled_final_metadata_summary %>%
  dplyr::filter(species != "DsecNoni",ClusterLabel %in% Anno_idents) %>%
  dplyr::select(-sem) %>%
  tidyr::spread(key="species", value="percent_combined") %>%
  dplyr::mutate(Dsec_Dmel=Dsec/Dmel, Dsec_Dsim=Dsec/Dsim, Dsim_Dmel = Dsim/Dmel) 

df_cell_comp %>%
  ggplot(.) +
  geom_point() +
  geom_smooth() +
  geom_label_repel(data=dplyr::filter(df_cell_comp, (abs(log2(Dsec_Dmel)) > 0.5 & abs(log2(Dsec_Dsim)) > 0.5)), 
                   force=5,
                   aes(label=ClusterLabel)) +
  aes(x=log2(Dsec_Dmel), y=log2(Dsec_Dsim)) +
  theme_bw()+
  lims(x=c(-1,1),y=c(-1,1))

#ggsave("Plots/TrioBrain_cell_comp_rho.png", width = 10, height = 8)

cor(df_cell_comp$Dsec_Dmel, df_cell_comp$Dsec_Dsim, method = 'spearman')


### glial freqeuncy ###

Seurat_obj_trio <- TrioBrain.integrated_slim_labeled_final
Idents(Seurat_obj_trio, cells = Glial_cells) <- 'Glia(all)'
Seurat_obj_trio[["ClusterLabel"]] <- Idents(object = Seurat_obj_trio)

Seurat_obj_trio_metadata <- Seurat_obj_trio@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))

Seurat_obj_trio_metadata_summary_rep <- Seurat_obj_trio_metadata %>%
  dplyr::group_by(species, ClusterLabel, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100) %>%
  dplyr::filter(ClusterLabel=='Glia(all)')

Seurat_obj_trio_metadata_summary <- Seurat_obj_trio_metadata_summary_rep %>%
  dplyr::group_by(species, ClusterLabel) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()%>%
  dplyr::filter(ClusterLabel=='Glia(all)')

df_freq_barplot <- rbind(df_TrioBrain.integrated_slim_labeled_final_metadata_summary, Seurat_obj_trio_metadata_summary)
df_freq_barplot_rep <- rbind(df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep, Seurat_obj_trio_metadata_summary_rep)

clu_sel <- c("α'β'-KC", "Clock_1", "Ach_2", "PRN", "Ach_1","Poxn_1", "Tk_1", "Glia(all)")

df_freq_barplot %>%
  dplyr::filter(species != "DsecNoni", ClusterLabel %in% clu_sel) %>%
  ggplot(.) +
  geom_col(aes(x=factor(ClusterLabel, levels=clu_sel), y=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.9)) + 
  geom_errorbar(aes(x=factor(ClusterLabel, levels=clu_sel), y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,0,percent_combined-sem), ymax=percent_combined+sem, group=species), 
                width=0.5, position=position_dodge(width=0.9), linewidth=0.2) +
  geom_point(data=dplyr::filter(df_freq_barplot_rep,
                                species != "DsecNoni", ClusterLabel %in% clu_sel),
             size = 1.6, alpha =0.8,  
             aes(x=ClusterLabel, y=percent, group=species,fill=species), shape=21, position=position_dodge(width=0.9)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=ClusterLabel, y=n/total*100)) +
  theme_bw() +
  theme(axis.title.y = element_text(size=13, color='black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x=element_text(size=12, color='black'),
        legend.text = element_text(size=13, face='italic', color='black'),
        panel.grid = element_blank(),
        legend.position=c(0.1,0.88),
        legend.background = element_blank()) +
  labs(x="Cluster", y="Frequency (% of cells)", fill="")  +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) +
  scale_y_continuous(limits = c(-0.05, 10), expand=c(0,0))

ggsave("Plots/Manuscript/Fig2_barplots.pdf", width=6, height = 4.2)

### ANOVA ###

df_TrioBrain_anova <- df_freq_barplot_rep %>%
  dplyr::filter(species != "DsecNoni",ClusterLabel %in% c(Anno_idents, "Glia(all)")) %>%
  dplyr::select(ClusterLabel, species, percent) %>%
  dplyr::group_by(ClusterLabel) %>%
  dplyr::summarise(p_value=tidy(aov(percent ~ species))$p.value[1]) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_value_adjust = p.adjust(p_value, method = "fdr")) %>%
  dplyr::mutate(sig = ifelse(p_value_adjust < 0.05, T, F))

df_TrioBrain_stat_split <- df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep %>%
  dplyr::filter(species != "DsecNoni",ClusterLabel %in% c(Anno_idents)) %>%
  dplyr::select(ClusterLabel, species, percent) %>%
  group_split(ClusterLabel)

clusters <- intersect(unique(df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep$ClusterLabel), Anno_idents)

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
  dplyr::select(cluster=ClusterLabel, Dsim_Dmel, Dsec_Dmel, Dsec_Dsim) %>%
  tidyr::gather(-cluster, key="pair", value="fold_diff") %>%
  dplyr::left_join(., df_comp_padj_gather, by=c('cluster', 'pair')) %>%
  dplyr::mutate(log2fc = log2(fold_diff)) %>%
  dplyr::group_by(pair) %>%
  dplyr::mutate(p_value_adjust = p.adjust(padj, method = "fdr")) %>%
  dplyr::mutate(sig = ifelse(p_value_adjust < 0.05, T, F)) %>%
  dplyr::ungroup()

df_cell_comp_padj %>%
  dplyr::filter(cluster %in% Anno_idents) %>%
  ggplot(.) +
  geom_vline(xintercept=0, size=0.2, alpha=0.5) +
  geom_hline(yintercept=-log10(0.05/107), size=0.2, alpha=0.5, linetype=2) +
  geom_point(size=2, alpha=0.8) +
  geom_text_repel(data=dplyr::filter(df_cell_comp_padj, cluster %in% c("PRN", "Tk_1", "Ach_1", 'Poxn_1')),
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
  scale_y_continuous(limit=c(-0.05,6), expand=c(0,0))

ggsave("Plots/Manuscript/Fig2_volcano.pdf", width=7.5, height = 5)

df_TrioBrain_anova_posthoc <- TrioBrain.integrated_slim_labeled_final_metadata_summary_rep %>%
  dplyr::filter(species != "DsecNoni") %>%
  dplyr::select(ClusterLabel, species, percent) %>%
  dplyr::group_by(ClusterLabel) %>%
  do(multitst = tidy(TukeyHSD(aov(percent ~ species, data = .))))


## UMAP plot for PRN, Ach_1, Poxn_1

### PRN ###

DimPlot(subset(TrioBrain.integrated_slim_labeled_final_species, ident = "PRN", orig.ident != "DsecNoni"), reduction = "tsne", 
        group.by = "orig.ident", shuffle = T, raster=F, pt.size = 1, 
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL) +
  coord_cartesian(xlim=c(12,15.7), ylim=c(-7.5,-10))+ NoAxes()

ggsave("Plots/Manuscript/Fig2_PRN.pdf", width=8, height = 8)

DimPlot(subset(TrioBrain.integrated_slim_labeled_final_species, ident = "PRN", orig.ident != "DsecNoni"), reduction = "tsne", 
        group.by = "orig.ident", split.by= "orig.ident", shuffle = T, raster=F, pt.size = 1, 
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL) +
  coord_cartesian(xlim=c(12,15.7), ylim=c(-7.5,-10))+ NoAxes()

ggsave("Plots/Manuscript/Fig2_PRN_species.pdf", width=16, height = 6)

### Ach_1 ###

DimPlot(subset(TrioBrain.integrated_slim_labeled_final_species, ident = "Ach_1", orig.ident != "DsecNoni"), reduction = "tsne", 
        group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.1, 
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL) +
  coord_cartesian(xlim=c(-15,10), ylim=c(-8.5,15))+ NoAxes()

ggsave("Plots/Manuscript/Fig2_Ach_1.pdf", width=8, height = 8)

DimPlot(subset(TrioBrain.integrated_slim_labeled_final_species, ident = "Ach_1", orig.ident != "DsecNoni"), reduction = "tsne", 
        group.by = "orig.ident", split.by= "orig.ident", shuffle = T, raster=F, pt.size = 0.1, 
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL) +
  coord_cartesian(xlim=c(-15,10), ylim=c(-8.5,15))+ NoAxes()

ggsave("Plots/Manuscript/Fig2_Ach_1_species.pdf", width=16, height = 6)

### Poxn_1 ###

DimPlot(subset(TrioBrain.integrated_slim_labeled_final_species, ident = "Poxn_1", orig.ident != "DsecNoni"), reduction = "tsne", 
        group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.5, 
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL) +
  coord_cartesian(xlim=c(-19.2,-16.8), ylim=c(6.8,9.5))+ NoAxes()

ggsave("Plots/Manuscript/Fig2_Poxn_1.pdf", width=8, height = 8)

DimPlot(subset(TrioBrain.integrated_slim_labeled_final_species, ident = "Poxn_1", orig.ident != "DsecNoni"), reduction = "tsne", 
        group.by = "orig.ident", split.by= "orig.ident", shuffle = T, raster=F, pt.size = 0.5, 
        cols=c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) + NoLegend() + labs(title=NULL)  +
  coord_cartesian(xlim=c(-19.2,-16.8), ylim=c(6.8,9.5))+ NoAxes()

ggsave("Plots/Manuscript/Fig2_Poxn_1_species.pdf", width=16, height = 6)

save(df_cell_comp_padj, df_freq_barplot, df_freq_barplot_rep, file="Processed_Data/freq_data.RData")