library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)
library(gtools)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

TrioBrain.integrated_slim_labeled_final <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_labeled_final.rds")
load("Processed_Data/celltype_order.RData")
load("Processed_Data/freq_data.RData")
load("Processed_Data/df_sum_exp_SCT_melrank.RData")
load("Processed_Data/ClusterLabel_cor_permutation.RData")
load("Processed_Data/df_deg_sigs.RData")
load("Processed_Data/Dsec_DEG.RData")
load("Processed_Data/df_sum_exp_reps_SCT.RData")

subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident == "Dsec")

##### Main Figures #####

### Figure 1 ###

TrioBrain.integrated_slim_labeled_final_species <- TrioBrain.integrated_slim_labeled_final

TrioBrain.integrated_slim_labeled_final_species$orig.ident <- gsub("_rep[1-6]", "", TrioBrain.integrated_slim_labeled_final_species$orig.ident)

TrioBrain.integrated_slim_labeled_final_species@meta.data$orig.ident <- factor(TrioBrain.integrated_slim_labeled_final_species@meta.data$orig.ident, 
                                                                 levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))

plot_dmel_tsne <- DimPlot(subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident == "Dmel"), 
                          reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01,
                          cols=c("#E41A1C")) + NoLegend() + labs(title=NULL) + NoAxes()

plot_dsim_tsne <- DimPlot(subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident == "Dsim"), 
                          reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01,
                          cols=c("#4DAF4A")) + NoLegend() + labs(title=NULL) + NoAxes()

plot_dsec_tsne <- DimPlot(subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident == "Dsec"), 
                          reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01,
                          cols=c("#377EB8")) + NoLegend() + labs(title=NULL) + NoAxes()

plot_trio_merged_tsne <- DimPlot(subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident != "DsecNoni"), 
                                 reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01,
                                 cols=c("#E41A1C", "#4DAF4A",  "#377EB8")) + NoLegend() + labs(title=NULL) + NoAxes()

plot_tsne_species <- cowplot::plot_grid(plot_dmel_tsne, plot_dsim_tsne,
                                        plot_dsec_tsne, plot_trio_merged_tsne, ncol=2)

plot_tsne_species

ggsave("Plots/Manuscript/Fig1_tsne_species.png", width=13, height = 12)

## Figure 1d - labeled by class ##

Ach_clusters <- mixedsort(Anno_idents[grep('Ach', Anno_idents)])
Glu_clusters <- mixedsort(Anno_idents[grep('Glu', Anno_idents)])
GABA_clusters <- mixedsort(Anno_idents[grep('GABA', Anno_idents)])
Glial_clusters <- c("AST","CTX","ENS","PRN","SUB")
KC_clusters <-  Anno_idents[grep('-KC', Anno_idents)]
MA_clusters <- c("DOP_1","DOP_2","DOP/SER", "SER", "TY_1", "TY_2", "TY_3", "OCTY", "Tbh")
Clock_clusters <-  mixedsort(Anno_idents[grep('Clock', Anno_idents)])
OPN_clusters <- 'OPN'
Poxn_clusters <-  mixedsort(Anno_idents[grep('Poxn', Anno_idents)])
fru_clusters <-  mixedsort(Anno_idents[grep('fru', Anno_idents)])
NP_clusters <- c( "CCHa2", "Dh31","FMRFa","Mip", "MS",  "RYa", "sNPF", "Tk_1",  "Tk_2", "Tk_3", 
                   "AstC_FMRFa","FMRFa_Dh44", "NPF_AstA", "NPF_Dh31",  "sNPF_Dh44", "Tk_sNPF",
                   "Tk_Mip_RYa", "sNPF_CCHa2_AstC")

Celltype_order_group <- c(Glial_clusters, KC_clusters, MA_clusters, Clock_clusters, OPN_clusters, Poxn_clusters, fru_clusters, NP_clusters, Ach_clusters, Glu_clusters, GABA_clusters)

save(Celltype_order_group, file="Processed_Data/Celltype_order_group.RData")

Ach_cells <- WhichCells(TrioBrain.integrated_slim_labeled_final_species, idents=Ach_clusters)
Glu_cells <- WhichCells(TrioBrain.integrated_slim_labeled_final_species, idents=Glu_clusters)
GABA_cells <- WhichCells(TrioBrain.integrated_slim_labeled_final_species, idents=GABA_clusters)
Glial_cells <- WhichCells(TrioBrain.integrated_slim_labeled_final_species, idents=Glial_clusters)
KC_cells <- WhichCells(TrioBrain.integrated_slim_labeled_final_species, idents=KC_clusters)
MA_cells <- WhichCells(TrioBrain.integrated_slim_labeled_final_species, idents=MA_clusters)
Clock_cells <- WhichCells(TrioBrain.integrated_slim_labeled_final_species, idents=Clock_clusters)
OPN_cells <- WhichCells(TrioBrain.integrated_slim_labeled_final_species, idents=OPN_clusters)
Poxn_cells <- WhichCells(TrioBrain.integrated_slim_labeled_final_species, idents=Poxn_clusters)
fru_cells <- WhichCells(TrioBrain.integrated_slim_labeled_final_species, idents=fru_clusters)
NP_cells <- WhichCells(TrioBrain.integrated_slim_labeled_final_species, idents=NP_clusters)

DimPlot(subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident != "DsecNoni"), 
        cells.highlight= list(Ach=Ach_cells, Glu=Glu_cells, GABA=GABA_cells, Glia=Glial_cells, 
                              KC=KC_cells, MA=MA_cells, Clock=Clock_cells, OPN=OPN_cells,
                              Poxn=Poxn_cells, fru=fru_cells, NP=NP_cells), 
        cols.highlight = c('#800000', '#006400', '#00008B', '#FF8C00', '#483D8B', '#2F4F4F', '#8B4513', '#228B22', '#191970', '#8B008B'),
        reduction = "tsne", group.by = "orig.ident", shuffle = T, raster=F, pt.size = 0.01, sizes.highlight = 0.01) + labs(title=NULL) + NoAxes()

ggsave("Plots/Manuscript/Fig1_tsne_anno.png", width=11, height = 9)

## Figure 1e - markers ##

Final_idents <- as.character(unique(Idents(TrioBrain.integrated_slim_labeled_final)))
exclude <- c(Final_idents[grep('Unanno', Final_idents)], 0:100)

Anno_idents <- setdiff(Final_idents,exclude)

TrioBrain.integrated_slim_labeled_final_anno_only <- subset(TrioBrain.integrated_slim_labeled_final, idents = Anno_idents)

Seurat_for_dotplot <- subset(TrioBrain.integrated_slim_labeled_final_anno_only, orig.ident != "DsecNoni")

Idents(Seurat_for_dotplot, cells = Ach_cells) <- 'Ach(all)'
Idents(Seurat_for_dotplot, cells = Glu_cells) <- 'Glu(all)'
Idents(Seurat_for_dotplot, cells = GABA_cells) <- 'GABA(all)'
Idents(Seurat_for_dotplot, cells = Clock_cells) <- 'Clock(all)'
Idents(Seurat_for_dotplot, cells = Poxn_cells) <- 'Poxn(all)'
Idents(Seurat_for_dotplot, cells = fru_cells) <- 'fru(all)'

Seurat_for_dotplot@active.ident <- factor(Seurat_for_dotplot@active.ident,
                                          levels=c('Ach(all)','Glu(all)','GABA(all)', "αβ-KC_1", "αβ-KC_2","γ-KC","α'β'-KC",
                                                   MA_clusters, 'Clock(all)', 'OPN','Poxn(all)','fru(all)',NP_clusters, Glial_clusters
                                                   ))

DotPlot(Seurat_for_dotplot, 
        dot.scale = 4.5,
        feature=c("pros","Imp", 
                  "VAChT","VGlut", "Gad1", 
                  "CG32204","crb","mamo","Pka-C1", "sNPF",
                  "Vmat", "DAT", "SerT", "Tdc2", "Tbh", 
                  "tim", "Clk", 'otp','acj6',"Poxn", "fru",  
                  "CCHa2", "Dh31","FMRFa","Mip", "Ms",  "RYa", "Tk","AstC","Dh44", "NPF", "AstA",
                  "alrm", "CG40470","axo","Gs2","Tret1-1", "baz")) +
  coord_flip() + labs(x=NULL, y="Clusters") +
  theme(legend.position='bottom',
        axis.title = element_blank(),
        axis.text.x = element_text(color='black', angle=45, hjust=1, vjust=1),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1))

ggsave("Plots/Manuscript/Fig1_cluster_markers.pdf", width=15, height = 8.2)


### Figure 2 ###

freq_summary <- df_freq_barplot %>%
  dplyr::filter(ClusterLabel %in% Anno_idents) %>%
  dplyr::select(CellType=ClusterLabel, species, Mean_Percent=percent_combined) %>%
  tidyr::spread(key=species, value=Mean_Percent) %>%
  dplyr::arrange(-Dmel)

write.table(freq_summary, file="Processed_Data/cell_type_freq.tsv", quote=F, row.names=F)

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

## DotPlot for Tret1-1 ##

DotPlot(subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident != "DsecNoni", idents = Anno_idents), feature="Tret1-1")

Dotplot_Tret_expression_data <- DotPlot(subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident != "DsecNoni", idents = Anno_idents), 
                                        feature='Tret1-1', split.by='orig.ident', cols=trio_col)$data %>%
  dplyr::mutate(species = ifelse(grepl('Dmel', id),'Dmel',
                                 ifelse(grepl('Dsim', id), 'Dsim', 'Dsec')))

Dotplot_Tret_expression_data$id <- gsub("_D(mel|sim|sec)", "", Dotplot_Tret_expression_data$id)

Dotplot_Tret_expression_data <- Dotplot_Tret_expression_data %>%  
  dplyr::rename(cluster=id) %>%
  dplyr::filter(cluster %in% Anno_idents) %>%
  dplyr::select(species, cluster, gene=features.plot, avg.exp, pct.exp) %>%
  na.omit()

Dotplot_Tret_expression_data$species <- factor(Dotplot_Tret_expression_data$species, levels=c('Dmel', 'Dsim', 'Dsec'))
Dotplot_Tret_expression_data$cluster <- factor(Dotplot_Tret_expression_data$cluster, levels=Anno_idents)

Dotplot_Tret_expression_data %>%
  ggplot(.) +
  geom_count() +
  aes(x=cluster, y=species, size=pct.exp, color=avg.exp) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', size=9, angle=45,  hjust=1, vjust=1),
        axis.text.y= element_text(color='black', size=9, face='italic'),
        axis.title.y= element_blank(),
        axis.title.x=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), face='italic', size=10),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
        strip.background = element_blank(),
        legend.title=element_text(color='black', size=9),
        legend.text=element_text(color='black', size=9),
        legend.key.size = unit(0.5, "lines"),
        legend.box = "horizontal") +
  scale_color_gradient(low = "lightgrey", high="blue")+
  scale_size_continuous(range = c(0.1, 4), breaks = c(0,10, 20,40,60,80)) +
  facet_wrap(~gene, ncol=1, scales='free') +
  labs(size='Frequency (%)', color='Expression')

ggsave("Plots/Manuscript/Tret1-1.pdf", width=14, height=2.5)

### Figure 3 ###

clu_exp_cell <- c("OPN",  "Clock_3", "γ-KC", "NPF_AstA")

df_sum_exp_SCT_melrank %>%
  dplyr::filter(cluster %in% clu_exp_cell) %>%
  ggplot(.) +
  geom_point(aes(x=log2(Dmel), y=log2(Dsec)), fill='orange', alpha =0.6, shape=21, stroke=0.3) +
  geom_point(aes(x=log2(Dmel), y=log2(Dsim)), fill='navy', alpha =0.6, shape=21, stroke=0.3) +
  geom_smooth(method='lm', aes(x=log2(Dmel), y=log2(Dsec)), color='orange', fill='orange', alpha =0.2, size=0.5) +
  geom_smooth(method='lm', aes(x=log2(Dmel), y=log2(Dsim)), color='navy', fill='navy', alpha =0.2, size=0.5) +
  geom_text_repel(data=dplyr::filter(df_sum_exp_SCT_melrank, cluster %in% clu_exp_cell, Dmel_rank <=3), 
            aes(label=gene, x=log2(Dmel), y=log2(Dsec)), color='orange', size=3, fontface='italic') +
  geom_text_repel(data=dplyr::filter(df_sum_exp_SCT_melrank, cluster %in% clu_exp_cell, Dmel_rank <=3), 
            aes(label=gene, x=log2(Dmel), y=log2(Dsim)), color='navy', size=3, fontface='italic') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=10, color='black'),
        axis.text = element_text(size=9.5, color='black'),
        strip.text = element_text(size=10, color='black')) +
  labs(x="log2(expression) - Dmel", y= "log2(expression) - Dsim or Dsec") +
  facet_wrap(~cluster, nrow=2, scales='free') 

ggsave("Plots/Manuscript/Fig3a_scatter_example.pdf", width=6, height=6)

non_neuronal_cluster <- c("AST","CTX","ENS","PRN","SUB")

df_per_cor_clean %>% 
  dplyr::filter(cluster %in% Anno_idents) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(melsec=mean(melsec), melsim=mean(melsim)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(group = ifelse(cluster %in% non_neuronal_cluster, "Non-neuronal", "Neuronal")) %>%
  ggplot(.) +
  geom_abline(intercept = 0, slope = 1, linetype=2, alpha = 0.5) +
  #geom_smooth(aes(x=melsec, y=melsim), method='lm') +
  geom_point(aes(x=melsec, y=melsim), shape=21, fill='red') +
  #geom_text_repel(aes(x=melsec, y=melsim, label=cluster), segment.size = 0.1, size = 3) +
  geom_text_repel(nudge_y = -0.02,aes(x=melsec, y=melsim, 
                      label=factor(cluster, levels =c("OPN",  "Clock_3", "γ-KC", "NPF_AstA")))) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position=c(0.85,0.1),
        legend.title=element_blank(),
        axis.title = element_text(size=10, color='black', face='italic'),
        axis.text = element_text(size=9.5, color='black'),
        legend.text = element_text(size=9.5, color='black'),
        legend.background = element_blank()) +
  scale_fill_manual(values=c("red", "blue")) +
  lims(x=c(0.38,0.8), y=c(0.38,0.8)) +
  labs(x="Dmel vs Dsec (ρ)", y="Dmel vs Dsim (ρ)")

ggsave(file="Plots/Manuscript/Fig3b_scatter_ClusterLabel.pdf", width=5.5, height=5.5)

df_order_melsec <- df_per_cor_clean_gather %>%
  dplyr::filter(pair == "melsec") %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(value=mean(cor_rho)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(value)

df_per_cor_clean_gather %>%
  dplyr::filter(cluster %in% c(Ach_clusters, Glu_clusters, GABA_clusters), pair != "secNoni") %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0.7, outlier.size = 0.5, size=0.3, alpha = 0.7) +
  #geom_smooth() +
  aes(x=factor(cluster,levels=df_order_melsec$cluster), y=cor_rho, fill=pair) +
  scale_fill_manual(values=c("orange", "navy",  "brown")) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position=c(0.92,0.15),
        legend.title=element_blank(),
        axis.title.y = element_text(size=10, color='black'),
        axis.text.x = element_text(size=9.5, color='black', angle=45, hjust =1, vjust=1),
        axis.text.y = element_text(size=9.5, color='black'),
        legend.text = element_text(size=9.5, color='black')) +
  labs(y="Spearman's rho (ρ)")

ggsave(file="Plots/Manuscript/Fig3c_boxplot_rho_NT.pdf", width=15, height=4.2)

df_per_cor_clean_gather %>%
  dplyr::filter(cluster %in% Anno_idents & !cluster %in% c(Ach_clusters, Glu_clusters, GABA_clusters), pair != "secNoni") %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0.7, outlier.size = 0.5, size=0.3, alpha = 0.7, linewidth=0.1) +
  #geom_smooth() +
  aes(x=factor(cluster,levels=df_order_melsec$cluster), y=cor_rho, fill=pair) +
  scale_fill_manual(values=c("orange", "navy",  "brown")) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position=c(0.92,0.15),
        legend.title=element_blank(),
        axis.title.y = element_text(size=10, color='black'),
        axis.text.x = element_text(size=9.5, color='black', angle=45, hjust =1, vjust=1),
        axis.text.y = element_text(size=9.5, color='black'),
        legend.text = element_text(size=9.5, color='black'),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1)) +
  labs(y="Spearman's rho (ρ)")

ggsave(file="Plots/Manuscript/Fig3c_boxplot_rho_nonNT.pdf", width=12.5, height=4)


### Figure 4 ###

### Volcano plot ###

## sim vs sec - PRN ##

plot_volcano <- function(sp_pair='simsec', target_cl = 'PRN', n_sig=10, 
                         data1=df_deg_clusters_ref_merge_sigs, data2=df_deg_clusters_ref_merge_min, xlims=c(-2,2), ylims=c(0,70)) {
  
  topdegs_simsec_fc <- df_deg_clusters_ref_merge_sigs %>%
    dplyr::filter(pair == sp_pair, cluster == target_cl, abs(avg_log2FC) >log2(1.5)) %>%
    dplyr::top_n(-p_val_adj, n=n_sig)
  
  simsec_clusters_deg <- df_deg_clusters_ref_merge_min %>%
    dplyr::filter(pair == sp_pair, cluster == target_cl)
  
  plot <- EnhancedVolcano(simsec_clusters_deg,
                  lab = simsec_clusters_deg$gene,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  selectLab = topdegs_simsec_fc$gene,
                  drawConnectors = TRUE,
                  widthConnectors = 0.2,
                  lengthConnectors = unit(0.01, "npc"),
                  xlim =xlims, 
                  ylim = ylims, 
                  gridlines.major = F, gridlines.minor = F,
                  pointSize = 1,
                  labSize = 3,
                  axisLabSize = 11,
                  legendLabSize = 10,
                  legendIconSize = 3,
                  borderWidth = 0.2,
                  vlineWidth = 0.2,
                  hlineWidth = 0.2,
                  cutoffLineWidth = 0.2,
                  FCcutoff = log2(1.5),
                  pCutoff = 0.05) +
    theme(axis.ticks = element_line(linewidth = 0.2),
          axis.line = element_line(linewidth = 0.1))
  
  return(plot)
  
}
plot_volcano(sp_pair='melsec', target_cl = 'PRN', xlims=c(-2.1,2.1), ylims=c(0,80))
ggsave("Plots/Manuscript/Fig4_Dmelsec_volcano_PRN.pdf", width=5, height=6.5)
plot_volcano(sp_pair='melsim', target_cl = 'PRN', xlims=c(-2.1,2.1), ylims=c(0,100))
ggsave("Plots/Manuscript/Fig4_Dmelsim_volcano_PRN.pdf", width=5, height=6.5)
plot_volcano(sp_pair='simsec', target_cl = 'PRN', xlims=c(-1.8,1.8), ylims=c(0,70))
ggsave("Plots/Manuscript/Fig4_Dsimsec_volcano_PRN.pdf", width=5, height=6.5)

plot_volcano(sp_pair='melsec', target_cl = 'Poxn_1', xlims=c(-6,6), ylims=c(0,120))
ggsave("Plots/Manuscript/Fig4_Dmelsec_volcano_Poxn_1.pdf", width=5, height=6.5)
plot_volcano(sp_pair='melsim', target_cl = 'Poxn_1',  xlims=c(-6.2,6.2), ylims=c(0,120))
ggsave("Plots/Manuscript/Fig4_Dmelsim_volcano_Poxn_1.pdf", width=5, height=6.5)
plot_volcano(sp_pair='simsec', target_cl = 'Poxn_1', xlims=c(-1.2,1.2), ylims=c(0,50))
ggsave("Plots/Manuscript/Fig4_Dsimsec_volcano_Poxn_1.pdf", width=5, height=6.5)

## group by n of cluster ##

df_deg_clusters_ref_merge_sigs_nc <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::group_by(pair, gene) %>%
  dplyr::summarise(n_cluster = n()) %>%
  dplyr::ungroup()

df_deg_clusters_ref_merge_sigs_nc %>%
  dplyr::filter( pair != "secNoni") %>%
  ggplot(.) +
  geom_histogram(stat='count', position = 'stack', binwidth = 1) +
  aes(x=n_cluster, fill=pair) +
  theme_bw() +
  theme(axis.title = element_text(size=8, color='black'),
        axis.text = element_text(size=7, color='black'),
        panel.grid=element_blank()) +
  scale_fill_manual(values= c("orange", "navy",  "brown")) +
  labs(x="Differentially expressed clusters", y= "DEG count") +
  scale_x_continuous(expand=c(0.01,0.01))

ggsave(file="Plots/Manuscript/Fig4_histogram_n_cluster_DEGs.pdf", width=3.5, height= 1.8)

## correlation between size and DEG ##

df_deg_clusters_ref_merge_sigs_summary_all <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::filter(pair != "secNoni") %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(mean_DEG_count=n()/3) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-mean_DEG_count)

df_cor_deg_cluster <- df_freq_barplot %>%
  dplyr::filter(species != "DsecNoni") %>%
  dplyr::rename(cluster=ClusterLabel) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(mean_percent=mean(percent_combined)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(., df_deg_clusters_ref_merge_sigs_summary_all, by='cluster') %>%
  na.omit()

df_cor_deg_cluster %>%
  ggplot(.) +
  geom_point(alpha=0.7) +
  geom_text_repel(data=dplyr::filter(df_cor_deg_cluster, mean_percent > 1 | mean_DEG_count > 40), aes(label=cluster),
                  color='blue', size=4, max.overlaps =10) +
  aes(x=mean_percent, y=mean_DEG_count) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.title=element_blank(),
        axis.title = element_text(size=12, color='black'),
        axis.text = element_text(size=11, color='black'),
        legend.background = element_blank()) +
  labs(x='Mean of cell type frequencies (%)', y='Mean of DEG counts')

ggsave("Plots/Manuscript/DEG_freq_comp.pdf", width=6, height=6)

## make plots with 40 to clusters

df_deg_clusters_ref_merge_sigs_summary <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(sum_DEG=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-sum_DEG) %>%
  dplyr::top_n(sum_DEG, n=40)

df_deg_clusters_ref_merge_sigs$cluster <- factor(df_deg_clusters_ref_merge_sigs$cluster, 
                                                 levels=df_deg_clusters_ref_merge_sigs_summary$cluster)

df_deg_clusters_ref_merge_sigs %>%
  dplyr::filter(cluster %in% Anno_idents, pair != "secNoni") %>%
  ggplot(.) +
  geom_histogram(stat='count', position = 'dodge') +
  aes(x=cluster, fill=pair) +
  theme_bw() +
  theme(axis.title.y = element_text(size=9, color='black'),
        axis.text.y = element_text(size=8, color='black'),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1, color='black', size=8),
        axis.title.x = element_blank(),
        panel.grid.major=element_blank()) +
  scale_fill_manual(values= c("orange", "navy",  "brown"))+
  labs(y= "DEG count")

ggsave(file="Plots/Manuscript/Fig4_histogram_cluster_DEGs.pdf", width=12, height= 4)


### overlap summary

df_overlap_deg_sum <- data.frame(overlap = factor(c("melsec only", "melsim only", "simsec only",
                                                    "melsec & melsim", "melsec & simsec", "melsim & simsec", "melsec & melsim & simsec"),
                                                  levels=c("melsec only", "melsim only", "simsec only",
                                                           "melsec & melsim", "melsec & simsec", "melsim & simsec", "melsec & melsim & simsec")),
                                 count = c(162,104,51, 195, 96, 48, 140))

df_overlap_deg_sum %>%
  ggplot(.) +
  geom_col() +
  aes(x=overlap, y=count) +
  theme_bw() +
  theme(axis.title.y = element_text(size=10, color='black'),
        axis.text.y = element_text(size=9, color='black'),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1, color='black', size=9, face='italic'),
        axis.title.x = element_blank(),
        panel.grid.major=element_blank())

ggsave(file="Plots/Manuscript/Fig4_col_cluster_DEGs_species_overlap.pdf", width=3, height= 4)


### Figure 5 ###

df_deg_Dsec_spec_summary1 <- df_deg_Dsec_spec %>%
  na.omit() %>%
  dplyr::filter(fc_test==T & padj_test==T) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n=n())

df_deg_Dsec_spec_summary1 %>%
  #dplyr::filter(!cluster %in% c("FMRFa(2)", "MON", "SUB", "Tbh", "Mip(2)", "Crz")) %>%
  ggplot(.) +
  #geom_point(shape=21) +
  #geom_line() +
  #aes(x=cluster, y=n_DEGs, fill=pair, group=pair) +
  geom_col() +
  aes(x=reorder(cluster,-n) , y=n) +
  theme_bw() +
  theme(axis.title.y = element_text(size=11, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1, color='black', size=9.5),
        axis.title.x = element_blank(),
        panel.grid.major=element_blank()) +
  scale_fill_manual(values= c("orange", "navy",  "brown"))+
  labs(y= "DEG count")

ggsave("Plots/Manuscript/Fig5a.pdf", width=18.5, height=5)


### Dsec_DEG info ###

TrioBrain.integrated_slim_labeled_trio <- subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident != "DsecNoni")
TrioBrain.integrated_slim_labeled_trio$orig.ident <- gsub("_rep[1-6]", "", TrioBrain.integrated_slim_labeled_trio$orig.ident)

trio_col <- c("#E41A1C", "#4DAF4A",  "#377EB8")

glia_meta <- c('Cat', 'dob', 'Sarm', 'CG4991') 
para1 <- c('CG9394', 'CG18135')
para2 <- c('CAH2', 'CAH3')
ecdy <- c('ImpE1', 'Eip63F-1')

## dot plot for 10 genes referred in the MS ##

### Differential expression - paralogs ###

select_genes <- c(glia_meta, para1, para2, ecdy)
n_select_genes <- length(select_genes)

Dotplot_expression_data <- DotPlot(TrioBrain.integrated_slim_labeled_trio, feature=select_genes, split.by='orig.ident', cols=trio_col)$data %>%
  dplyr::mutate(species = ifelse(grepl('Dmel', id),'Dmel',
                                 ifelse(grepl('Dsim', id), 'Dsim', 'Dsec')))

Dotplot_expression_data$id <- gsub("_D(mel|sim|sec)", "", Dotplot_expression_data$id)

Dotplot_expression_data <- Dotplot_expression_data %>%  
  dplyr::rename(cluster=id) %>%
  dplyr::filter(cluster %in% Anno_idents) %>%
  dplyr::select(species, cluster, gene=features.plot, avg.exp, pct.exp) %>%
  na.omit()

load("Processed_Data/Celltype_order_group.RData")

Dotplot_expression_data$species <- factor(Dotplot_expression_data$species, levels=c('Dmel', 'Dsim', 'Dsec'))
Dotplot_expression_data$cluster <- factor(Dotplot_expression_data$cluster, levels=Celltype_order_group)

list_dots <- list()

for (i in 1:(n_select_genes)) {
  
  list_dots[[i]] <- Dotplot_expression_data %>%
    dplyr::filter(gene == select_genes[i]) %>%
    ggplot(.) +
    geom_count() +
    aes(x=cluster, y=species, size=pct.exp, color=avg.exp) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x= element_text(color='black', size=9.5, angle=45,  hjust=1, vjust=1),
          axis.text.y= element_text(color='black', size=9, face='italic'),
          axis.title.y= element_blank(),
          axis.title.x=element_blank(), 
          strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), face='italic', size=12),
          strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=12),
          strip.background = element_blank(),
          legend.title=element_text(color='black', size=10),
          legend.text=element_text(color='black', size=10),
          legend.key.size = unit(1, "lines"),
          legend.key.height = unit(0.5, "lines"),
          legend.box = "horizontal", 
          legend.position='bottom') +
    scale_color_gradient(low = "lightgrey", high="blue")+
    scale_size_continuous(range = c(0.1, 4), breaks = c(0,10,20, 40,80)) +
    facet_wrap(~gene, ncol=1, scales='free') +
    labs(size='Frequency (%)', color='Expression')
  
  
}

plot_Dsec_degs_dots <- cowplot::plot_grid(plotlist=list_dots, ncol=1, align = 'v', axis="btlr")

plot_Dsec_degs_dots

ggsave(file="Plots/Manuscript/Fig5_dots_raw.pdf", width=22, height= 25)


### Figure 6 ###

secNoni_clusters_deg <- df_deg_clusters_ref_merge_min %>%
  dplyr::filter(pair == 'secNoni', gene != "Hr38", p_val_adj != 1) %>%
  dplyr::mutate(avg_log2FC=-avg_log2FC)

EnhancedVolcano(secNoni_clusters_deg,
                lab = secNoni_clusters_deg$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = NULL,
                drawConnectors = TRUE,
                widthConnectors = 0.25,
                lengthConnectors = unit(0.01, "npc"),
                xlim = c(-1.1, 1.1), 
                gridlines.major = F, gridlines.minor = F,
                pointSize = 1,
                labSize = 3,
                axisLabSize = 11,
                legendLabSize = 10,
                legendIconSize = 3,
                borderWidth = 0.6,
                FCcutoff = log2(1.2),
                pCutoff = 0.05)

ggsave("Plots/Manuscript/Fig6_secNoni_volcano12.pdf", width=5, height=6.5)


Plot_PRN <- df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="PRN", gene %in% c("CG5151")) %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.5, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=9),
        axis.text.y= element_text(color='black', size=9),
        axis.title=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8", "purple"))+
  facet_wrap(cluster~gene, scales = 'free', nrow=1) 

Plot_Tk_sNPF <- df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="Tk_sNPF", gene %in% c("Mkp3")) %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.5, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=9),
        axis.text.y= element_text(color='black', size=9),
        axis.title=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8", "purple"))+
  facet_wrap(cluster~gene, scales = 'free', nrow=1) 

Fig5c <-cowplot::plot_grid(Plot_PRN,Plot_Tk_sNPF, ncol=2,  axis='tblr')
Fig5c

ggsave("Plots/Manuscript/Fig6_dotplot_plasticity.pdf", width=6, height =3)

