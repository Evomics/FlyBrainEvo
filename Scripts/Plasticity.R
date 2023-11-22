library(Seurat)
library(tidyverse)
library(cowplot)
library(EnhancedVolcano)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_per_expressed_ref.RData")
load("Processed_Data/shared_genes.RData")
load("Processed_Data/celltype_order.RData")
load("Processed_Data/df_sum_exp_SCT.RData")
load("Processed_Data/celltype_cor_permutation.RData")
load("Processed_Data/df_sum_exp_reps_SCT.RData")
load("Processed_Data/df_deg_sigs.RData")
load("Processed_Data/DEGs/cluster_specific_DEGs_refmelsim.RData")
load("Processed_Data/global_DEGs.RData")

TrioBrain.integrated_slim_labeled <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled.rds")

## frequency ##

df_TrioBrain.integrated_slim_labeled_metadata <- TrioBrain.integrated_slim_labeled@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))

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

df_TrioBrain.integrated_slim_labeled_metadata_summary_rep_noni <- df_TrioBrain.integrated_slim_labeled_metadata_summary_rep %>%
  dplyr::filter(species %in% c("DsecNoni","Dsec"))

df_TrioBrain.integrated_slim_labeled_metadata_summary %>%
  dplyr::filter(species %in% c("DsecNoni","Dsec")) %>%
  ggplot(.) +
  #geom_col(aes(x=species, y=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.8)) + 
  geom_errorbar(aes(x=CellType, y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,log10(0.00000000000001),log10(percent_combined-sem)), 
                    ymax=log10(percent_combined+sem), group=species), width=0.8, position=position_dodge(width=0.8), linewidth=0.1) +
  geom_point(data=df_TrioBrain.integrated_slim_labeled_metadata_summary_rep_noni,
             size = 1.3, alpha =0.8,  
             aes(x=CellType, y=log10(percent), group=species,fill=species), shape=21, position=position_dodge(width=0.8), stroke=0.1) +
  geom_point(aes(x=CellType, y=log10(percent_combined), group=species,fill=species), shape=3, position=position_dodge(width=0.8), stroke=0.1) +
#  stat_summary(data=df_TrioBrain.integrated_slim_labeled_metadata_summary_rep_noni, 
#               fun = "mean", colour = "black", size = 2, geom = "point", shape=3,  aes(x=CellType, y=n/total*100, fill=species), position=position_dodge(width=0.8)) +
  theme_bw() +
  theme(axis.title.y = element_text(size=12, color='black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=11, color='black'),
        axis.text.x=element_text(size=10, color='black', angle=45, hjust=1, vjust=1),
        legend.text = element_text(size=12, face='italic', color='black'),
        panel.grid = element_blank(),
        strip.text = element_text(size=11, color='black'),
        legend.position='none') +
  labs(x="Cell type", y="log10(Frequency (% of cells))", fill="")  +
  scale_fill_manual(values= c( "#377EB8","#984EA3"))+
  ylim(-2,2)

ggsave("Plots/DsecNoni_freq.pdf", width=12, height=5)

### significance test ###

df_noni_t_test <- df_TrioBrain.integrated_slim_labeled_metadata_summary_rep %>%
  dplyr::filter(species %in% c("Dsec","DsecNoni")) %>%
  tidyr::separate(orig.ident, into=c('condition','rep'), sep='_') %>%
  dplyr::select(CellType, condition, rep, percent) %>%
  tidyr::spread(key='condition', value='percent') %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(p_value = t.test(Dsec, DsecNoni)$p.value) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(padj = p.adjust(p_value, method='fdr'))
  



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

noni_DEGs <- df_deg_clusters_ref_merge %>%
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > log2(1.2), 
                pair=="secNoni", gene != "Hr38")

select_genes <- unique(noni_DEGs$gene)
select_clusters <- unique(noni_DEGs$cluster)

df_sum_exp_reps_SCT$species <- factor(df_sum_exp_reps_SCT$species, levels=c("Dmel", "Dsim","Dsec","DsecNoni"),
                                      labels = c("Dmel", "Dsim","Dsec","Dsec\n(+Noni)"))

## stat test ##

df_sum_exp_reps_SCT_sNPF <- df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="sNPF(E)", gene == "sNPF")

anova_result <- aov(expression ~ species, data = df_sum_exp_reps_SCT_sNPF)
tukey_result <- TukeyHSD(anova_result)
tukey_result


Plot_noni_AST <- df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="AST", gene %in% c("CG14989", "CG43740", "Treh")) %>%
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
  facet_wrap(cluster~factor(gene, levels= c("CG14989", "CG43740", "Treh")), scales = 'free', nrow=1) 

Plot_ENS <- df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="ENS", gene %in% c("CG31324", "CG43740")) %>%
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
  facet_wrap(cluster~factor(gene, levels= c("CG31324", "CG43740")), scales = 'free', nrow=1) 

Plot_Proc <- df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="Proc", gene %in% c("Src42A")) %>%
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

Plot_sNPF <- df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="sNPF(E)", gene %in% c("sNPF")) %>%
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

Fig5d1 <- cowplot::plot_grid(Plot_PRN, Plot_noni_AST, rel_widths = c(1,3))
Fig5d2 <- cowplot::plot_grid(Plot_ENS, Plot_sNPF, Plot_Proc, rel_widths = c(2,1,1), nrow=1)

Fig5d <-cowplot::plot_grid(Fig5d1, Fig5d2, ncol=1, align ='hv', axis='tblr')
Fig5d

ggsave(glue::glue("Plots/Manuscript/Fig6_dotplot_plasticity.pdf"), width=8, height =4.5)


### dot plot across all cell types ###

## Supplemetnary Fig -  CG5151 expression pattern ##

Plot_Dsec_CG5151 <- DotPlot(subset(TrioBrain.integrated_slim_labeled_species, orig.ident == "Dsec", idents = anno_cluster), feature="CG5151") + scale_radius(limits=c(0,80), range=c(0,6))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=11, color='black'),
        legend.text = element_text(size=10, color='black'),
        legend.title = element_text(size=11, color='black'),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.2, 'in'),
        legend.direction = 'vertical',
        legend.position = 'right')

Plot_DsecNoni_CG5151 <- DotPlot(subset(TrioBrain.integrated_slim_labeled_species, orig.ident == "DsecNoni", idents = anno_cluster), feature="CG5151") + scale_radius(limits=c(0,70), range=c(0,6))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=11, color='black'),
        legend.text = element_text(size=10, color='black'),
        legend.title = element_text(size=11, color='black'),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.2, 'in'),
        legend.direction = 'vertical',
        legend.position = 'right')

Plot_CG5151 <- cowplot::plot_grid(Plot_Dsec_CG5151, Plot_DsecNoni_CG5151, nrow=1)               

Plot_CG5151               

ggsave(file="Plots/CG5151.pdf", width=8, height=8)
