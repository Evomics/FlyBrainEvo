library(Seurat)
library(tidyverse)
library(ggplot2)

TrioBrain.integrated_slim_labeled <- readRDS(file = "../brain/Processed_Data/TrioBrain_DF.integrated_slim_labeled.rds")

old_sNPFe_cellid <- WhichCells(TrioBrain.integrated_slim_labeled, idents='sNPF(E)')
old_Dh44_cellid <- WhichCells(TrioBrain.integrated_slim_labeled, idents='Dh44')

TrioBrain.integrated_slim_labeled_final <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_labeled_final.rds")

sNPFe_old <- subset(TrioBrain.integrated_slim_labeled_final, cells = old_sNPFe_cellid)

DimPlot(sNPFe_old, label=T)

sNPFe_old[["ClusterLabel"]] <- Idents(object = sNPFe_old)

df_sNPFe_old_metadata <- sNPFe_old@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))

df_sNPFe_old_metadata_summary_rep <- df_sNPFe_old_metadata %>%
  dplyr::group_by(species, ClusterLabel, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

df_sNPFe_old_metadata_summary <- df_sNPFe_old_metadata_summary_rep %>%
  dplyr::group_by(species, ClusterLabel) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

df_sNPFe_old_metadata_summary %>%
  dplyr::filter(species != "DsecNoni") %>%
  ggplot(.) +
  geom_col(aes(x=ClusterLabel, y=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.9)) + 
  geom_errorbar(aes(x=ClusterLabel, y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,0,percent_combined-sem), ymax=percent_combined+sem, group=species), width=0.9, position=position_dodge(width=0.9)) +
  geom_point(data=dplyr::filter(df_sNPFe_old_metadata_summary_rep,species != "DsecNoni"),
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

sNPF_new <- subset(TrioBrain.integrated_slim_labeled_final,idents = "sNPF")

### sNPF(E) -> sNPF / fru_2, sNPF frequency difference looks similar (fru_2 no difference bewteen species)

### Dh44 ###

Dh44_old <- subset(TrioBrain.integrated_slim_labeled_final, cells = old_Dh44_cellid)

DimPlot(Dh44_old, label=T)

FeaturePlot(Dh44_old, feature="Dh44")

Dh44_old[["ClusterLabel"]] <- Idents(object = Dh44_old)

df_Dh44_old_metadata <- Dh44_old@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))

df_Dh44_old_metadata_summary_rep <- df_Dh44_old_metadata %>%
  dplyr::group_by(species, ClusterLabel, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

df_Dh44_old_metadata_summary <- df_Dh44_old_metadata_summary_rep %>%
  dplyr::group_by(species, ClusterLabel) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

df_Dh44_old_metadata_summary %>%
  dplyr::filter(species != "DsecNoni") %>%
  ggplot(.) +
  geom_col(aes(x=ClusterLabel, y=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.9)) + 
  geom_errorbar(aes(x=ClusterLabel, y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,0,percent_combined-sem), ymax=percent_combined+sem, group=species), width=0.9, position=position_dodge(width=0.9)) +
  geom_point(data=dplyr::filter(df_Dh44_old_metadata_summary_rep,species != "DsecNoni"),
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

### Dh44 -> Ach10, frequency difference looks similar





