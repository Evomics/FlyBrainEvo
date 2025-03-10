library(Seurat)
library(tidyverse)
library(cowplot)
library(EnhancedVolcano)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))


load("Processed_Data/df_deg_sigs.RData")
load("Processed_Data/shared_genes.RData")
load("Processed_Data/Celltype_order_group.RData")

TrioBrain.integrated_slim_labeled_final <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_labeled_final.rds")
TrioBrain.integrated_slim_labeled_final_species <- TrioBrain.integrated_slim_labeled_final
TrioBrain.integrated_slim_labeled_final_species$orig.ident <- gsub("_rep[1-6]", "", TrioBrain.integrated_slim_labeled_final_species$orig.ident)
TrioBrain.integrated_slim_labeled_final_species@meta.data$orig.ident <- factor(TrioBrain.integrated_slim_labeled_final_species@meta.data$orig.ident, 
                                                                               levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))

load("Processed_Data/df_sum_exp_SCT.RData")
load("Processed_Data/cluster_specific_DEGs_ref_self.RData")
load("Processed_Data/df_per_expressed_ref.RData")
load("Processed_Data/celltype_order.RData")
load("Processed_Data/celltype_cor_permutation.RData")

## volcano plots ##

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


### noni DEGs - lowered threshold ###

noni_DEGs <- df_deg_clusters_ref_merge_min %>%
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > log2(1.2), 
                pair=="secNoni", gene != "Hr38")

select_genes <- unique(noni_DEGs$gene)
select_clusters <- unique(noni_DEGs$cluster)

sum_exp_reps <- AverageExpression(TrioBrain.integrated_slim_labeled_final, features = shared_genes, slot = 'data', add.ident = 'orig.ident')

df_sum_exp_reps_SCT <- as.data.frame(sum_exp_reps$SCT)  %>%
  tibble::rownames_to_column('gene') %>%
  tidyr::gather(-gene, key='cluster', value='expression') %>%
  dplyr::mutate(species=ifelse(grepl("Dmel", cluster)==T, "Dmel",
                               ifelse(grepl("Dsim", cluster)==T, "Dsim",
                                      ifelse(grepl("DsecNoni", cluster)==T, "DsecNoni","Dsec")))) %>%
  dplyr::mutate(reps = ifelse(grepl("rep1", cluster)==T, 1,
                              ifelse(grepl("rep2", cluster)==T, 2,
                                     ifelse(grepl("rep3", cluster)==T, 3,
                                            ifelse(grepl("rep4", cluster)==T, 4,
                                                   ifelse(grepl("rep5", cluster)==T, 5,6))))))

df_sum_exp_reps_SCT$cluster <- gsub("_Dmel_rep.", "", df_sum_exp_reps_SCT$cluster)
df_sum_exp_reps_SCT$cluster <- gsub("_Dsim_rep.", "", df_sum_exp_reps_SCT$cluster)
df_sum_exp_reps_SCT$cluster <- gsub("_DsecNoni_rep.", "", df_sum_exp_reps_SCT$cluster)
df_sum_exp_reps_SCT$cluster <- gsub("_Dsec_rep.", "", df_sum_exp_reps_SCT$cluster)

unique(df_sum_exp_reps_SCT$cluster)

df_sum_exp_reps_SCT$species <- factor(df_sum_exp_reps_SCT$species, levels = c("Dmel", "Dsim", "Dsec", "DsecNoni"))

save(Anno_idents, df_sum_exp_reps_SCT, file="Processed_Data/df_sum_exp_reps_SCT.RData")
load("Processed_Data/df_sum_exp_reps_SCT.RData")

df_sum_exp_reps_SCT$species <- factor(df_sum_exp_reps_SCT$species, levels=c("Dmel", "Dsim","Dsec","DsecNoni"),
                                      labels = c("Dmel", "Dsim","Dsec","Dsec\n(+Noni)"))


Plot_noni_AST <- df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="AST", gene %in% c("Treh")) %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.5, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=9),
        axis.text.y= element_text(color='black', size=9),
        axis.title=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8", "purple"))+
  facet_wrap(cluster~gene, scales = 'free', nrow=1) 

Plot_ENS <- df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="ENS", gene %in% c("Cipc", "SNF4Agamma")) %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.5, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=9),
        axis.text.y= element_text(color='black', size=9),
        axis.title=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8", "purple"))+
  facet_wrap(cluster~factor(gene, levels= c("Cipc", "SNF4Agamma")), scales = 'free', nrow=1) 

Plot_Poxn_3 <- df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="Poxn_3", gene %in% c("Vha100-1")) %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.5, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=9),
        axis.text.y= element_text(color='black', size=9),
        axis.title=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8", "purple"))+
  facet_wrap(cluster~gene, scales = 'free', nrow=1) 

Plot_Poxn_3

plot_plasticity1 <- cowplot::plot_grid(Plot_noni_AST, Plot_Poxn_3, rel_widths = c(1,1), nrow=1)
cowplot::plot_grid(Plot_ENS, plot_plasticity1, ncol=1, align='hv')

ggsave(glue::glue("Plots/Manuscript/FigS15_dotplot_plasticity.pdf"), width=6, height =6)


Plot_Tk_sNPF <- df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="Tk_sNPF", gene %in% c("Mkp3")) %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.5, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=9),
        axis.text.y= element_text(color='black', size=9),
        axis.title=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8", "purple"))+
  facet_wrap(cluster~gene, scales = 'free', nrow=1) 

Plot_PRN <- df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="PRN", gene %in% c("CG5151")) %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.5, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=9),
        axis.text.y= element_text(color='black', size=9),
        axis.title=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8", "purple"))+
  facet_wrap(cluster~gene, scales = 'free', nrow=1) 

Fig5d1 <- cowplot::plot_grid(Plot_PRN, Plot_noni_AST, rel_widths = c(1,1), nrow=1)
Fig5d3 <- cowplot::plot_grid(Plot_Poxn_3, Plot_Tk_sNPF, rel_widths = c(1,1), nrow=1)

Fig5d <-cowplot::plot_grid(Fig5d1,Plot_ENS, Fig5d3, ncol=1,  axis='tblr', rel_heights = c(1,1,1))
Fig5d

ggsave(glue::glue("Plots/Manuscript/Fig6_dotplot_plasticity.png"), width=8, height =12)



## stat test ##

df_sum_exp_reps_SCT_sNPF <- df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="sNPF(E)", gene == "sNPF")

anova_result <- aov(expression ~ species, data = df_sum_exp_reps_SCT_sNPF)
tukey_result <- TukeyHSD(anova_result)
tukey_result


## frequency ##

load("Processed_Data/freq_data.RData")

df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep_noni <- df_freq_barplot_rep %>%
  dplyr::filter(species %in% c("DsecNoni","Dsec"), ClusterLabel %in% Anno_idents)

df_freq_barplot %>%
  dplyr::filter(species %in% c("DsecNoni","Dsec"), ClusterLabel %in% Anno_idents) %>%
  ggplot(.) +
  #geom_col(aes(x=species, y=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.8)) + 
  geom_errorbar(aes(x=ClusterLabel, y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,log10(0.00000000000001),log10(percent_combined-sem)), 
                    ymax=log10(percent_combined+sem), group=species), width=0.8, position=position_dodge(width=0.8), linewidth=0.2) +
  geom_point(data=df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep_noni,
             size = 2, alpha =0.8,  
             aes(x=ClusterLabel, y=log10(percent), group=species,fill=species), shape=21, position=position_dodge(width=0.8), stroke=0.1) +
  geom_point(aes(x=ClusterLabel, y=log10(percent_combined), group=species,fill=species), shape=3, position=position_dodge(width=0.8), stroke=0.1) +
#  stat_summary(data=df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep_noni, 
#               fun = "mean", colour = "black", size = 2, geom = "point", shape=3,  aes(x=ClusterLabel, y=n/total*100, fill=species), position=position_dodge(width=0.8)) +
  theme_bw() +
  theme(axis.title.y = element_text(size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=9, color='black'),
        axis.text.x=element_text(size=8.5, color='black', angle=45, hjust=1, vjust=1),
        legend.text = element_text(size=11, face='italic', color='black'),
        panel.grid = element_blank(),
        strip.text = element_text(size=11, color='black'),
        legend.position='none') +
  labs(x="Cell type", y="log10(Frequency (% of cells))", fill="")  +
  scale_fill_manual(values= c( "#377EB8","#984EA3"))+
  ylim(-2,2)

ggsave("Plots/Manuscript/DsecNoni_freq.pdf", width=20, height=5)

### significance test ###

df_noni_t_test <- df_TrioBrain.integrated_slim_labeled_final_metadata_summary_rep_noni %>%
  dplyr::filter(species %in% c("DsecNoni","Dsec"), ClusterLabel %in% Anno_idents) %>%
  tidyr::separate(orig.ident, into=c('condition','rep'), sep='_') %>%
  dplyr::select(ClusterLabel, condition, rep, percent) %>%
  tidyr::spread(key='condition', value='percent') %>%
  dplyr::group_by(ClusterLabel) %>%
  dplyr::summarise(p_value = t.test(Dsec, DsecNoni)$p.value) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(padj = p.adjust(p_value, method='fdr'))
  


### dot plot across all cell types ###

## Supplemetnary Fig -  CG5151 expression pattern ##

Dsec_Seurat <- subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident %in% c("Dsec", "DsecNoni"), idents = Celltype_order_group)

Dotplot_expression_data_Dsec <- DotPlot(Dsec_Seurat, feature='CG5151', split.by='orig.ident', cols=c("#377EB8","#984EA3"))$data %>%
  dplyr::mutate(species = ifelse(grepl('Noni', id),'DsecNoni', 'Dsec'))

Dotplot_expression_data_Dsec$id <- gsub("_DsecNoni", "", Dotplot_expression_data_Dsec$id)
Dotplot_expression_data_Dsec$id <- gsub("_Dsec", "", Dotplot_expression_data_Dsec$id)

Dotplot_expression_data_Dsec <- Dotplot_expression_data_Dsec %>%  
  dplyr::rename(cluster=id) %>%
  dplyr::filter(cluster %in% Anno_idents) %>%
  dplyr::select(species, cluster, gene=features.plot, avg.exp, pct.exp) %>%
  na.omit()

Dotplot_expression_data_Dsec$cluster <- factor(Dotplot_expression_data_Dsec$cluster, levels=Celltype_order_group)

Dotplot_expression_data_Dsec %>%
  ggplot(.) +
  geom_count() +
  aes(x=cluster, y=species, size=pct.exp, color=avg.exp) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', size=9, angle=45,  hjust=1, vjust=1),
        axis.text.y= element_text(color='black', size=10, face='italic'),
        axis.title.y= element_blank(),
        axis.title.x=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), face='italic', size=12),
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

ggsave(file="Plots/Manuscript/CG5151_dotplot.pdf", width=22, height= 2)

Dsec_Seurat <- subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident == "Dsec", idents = Celltype_order_group)
DsecNoni_Seurat <- subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident == "DsecNoni", idents = Celltype_order_group)

Dsec_Seurat@active.ident <- factor(Dsec_Seurat@active.ident,levels=Celltype_order_group)
DsecNoni_Seurat@active.ident <- factor(DsecNoni_Seurat@active.ident,levels=Celltype_order_group)

Plot_Dsec_CG5151 <- DotPlot(Dsec_Seurat, feature="CG5151") + scale_radius(limits=c(0,80), range=c(0,6))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=11, color='black'),
        legend.text = element_text(size=10, color='black'),
        legend.title = element_text(size=11, color='black'),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.2, 'in'),
        legend.direction = 'vertical',
        legend.position = 'right') +
  scale_color_gradient2(low=-2)

Plot_DsecNoni_CG5151 <- DotPlot(DsecNoni_Seurat, feature="CG5151") + scale_radius(limits=c(0,70), range=c(0,6))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=11, color='black'),
        legend.text = element_text(size=10, color='black'),
        legend.title = element_text(size=11, color='black'),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.2, 'in'),
        legend.direction = 'vertical',
        legend.position = 'right') +
  scale_fill_continuous(lims=c(-2,2))

Plot_CG5151 <- cowplot::plot_grid(Plot_Dsec_CG5151, Plot_DsecNoni_CG5151, nrow=1)               

Plot_CG5151               

ggsave(file="Plots/Manuscript/CG5151.pdf", width=8, height=8)

## plot HCR

Plot_Dsec_HCR <- DotPlot(subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident == "Dsec", idents = Anno_idents), feature=c("CG5151", "Tret1-1", "ImpE1")) + scale_radius(limits=c(0,80), range=c(0,6))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=11, color='black'),
        legend.text = element_text(size=10, color='black'),
        legend.title = element_text(size=11, color='black'),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.2, 'in'),
        legend.direction = 'vertical',
        legend.position = 'right')

Plot_DsecNoni_HCR <- DotPlot(subset(TrioBrain.integrated_slim_labeled_final_species, orig.ident == "DsecNoni", idents = Anno_idents), feature=c("CG5151", "Tret1-1", "ImpE1")) + scale_radius(limits=c(0,70), range=c(0,6))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=11, color='black'),
        legend.text = element_text(size=10, color='black'),
        legend.title = element_text(size=11, color='black'),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.2, 'in'),
        legend.direction = 'vertical',
        legend.position = 'right')

Plot_HCR <- cowplot::plot_grid(Plot_Dsec_HCR, Plot_DsecNoni_HCR, nrow=1)               

Plot_HCR               

ggsave(file="Plots/HCR.pdf", width=8, height=8)

df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="PRN", gene %in% c("Tret1-1", "CG5151")) %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.5, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=9),
        axis.text.y= element_text(color='black', size=9),
        axis.title=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8", "purple"))+
  facet_wrap(cluster~gene, scales = 'free', nrow=1) 

df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="ENS", gene %in% c("ImpE1")) %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.5, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=9),
        axis.text.y= element_text(color='black', size=9),
        axis.title=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8", "purple"))+
  facet_wrap(cluster~gene, scales = 'free', nrow=1) 

