library(Seurat)
library(tidyverse)
library(ggrepel)
library(EnhancedVolcano)
library(data.table)
library(cowplot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/shared_genes.RData")
load("Processed_Data/ClusterLabel_order.RData")
load("Processed_Data/cluster_specific_DEGs_refmel.RData")
load("Processed_Data/cluster_specific_DEGs_ref_self.RData")
load("Processed_Data/df_per_expressed_ref.RData")

 ######### cell-type specific #########

df_deg_clusters_ref_merge1 <- df_deg_clusters_ref_mel %>%
  dplyr::mutate(ref="MEL")

df_deg_clusters_ref_merge2 <- df_deg_clusters_ref_self %>%
  dplyr::mutate(ref="OWN")

df_deg_clusters_ref_merge <- rbind(df_deg_clusters_ref_merge1, df_deg_clusters_ref_merge2)

### take min log2FC from two refs ##

df_deg_clusters_ref_merge_min <- df_deg_clusters_ref_merge %>%
  dplyr::group_by(pair, gene, cluster) %>%
  dplyr::mutate(avg_log2FC_min=min(abs(avg_log2FC))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(abs(avg_log2FC)==avg_log2FC_min) %>%
  dplyr::mutate(fc_test = ifelse(abs(avg_log2FC) > log2(1.5),T,F), padj_test = ifelse(p_val_adj < 0.05, T,F))

df_deg_clusters_ref_merge_sigs <- df_deg_clusters_ref_merge_min %>%
  dplyr::filter(fc_test==T & padj_test==T)

save(df_deg_clusters_ref_merge_min, df_deg_clusters_ref_merge_sigs, 
     file="Processed_Data/df_deg_sigs.RData")

#load("Processed_Data/df_deg_sigs.RData")

df_deg_clusters_ref_merge_sigs %>%
  dplyr::mutate(direction = ifelse(avg_log2FC>0, "up", "down")) %>%
  dplyr::group_by(cluster, pair, direction) %>%
  dplyr::summarise(n=n())

df_deg_clusters_summary <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::distinct(gene, CellType=cluster, pair, avg_log2FC, p_val_adj) %>%
  dplyr::arrange(gene, pair)

write.csv(df_deg_clusters_summary, file="Processed_Data/DEG_list.csv", quote=F, row.names=F)

### Volcano plot ###

## sim vs sec - PRN ##

topdegs_simsec_fc_PRN <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::filter(pair == 'simsec', cluster == "PRN")

simsec_clusters_deg_PRN <- df_deg_clusters_ref_merge_min %>%
  dplyr::filter(pair == 'simsec', cluster == "PRN")

EnhancedVolcano(simsec_clusters_deg_PRN,
                lab = simsec_clusters_deg_PRN$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = topdegs_simsec_fc_PRN$gene,
                drawConnectors = TRUE,
                widthConnectors = 0.25,
                lengthConnectors = unit(0.01, "npc"),
                xlim = c(-2,2), 
                ylim = c(0,32), 
                gridlines.major = F, gridlines.minor = F,
                pointSize = 1,
                labSize = 3,
                axisLabSize = 11,
                legendLabSize = 10,
                legendIconSize = 3,
                borderWidth = 0.6,
                FCcutoff = log2(1.5),
                pCutoff = 0.05)

ggsave("Plots/Manuscript/Fig4_Dsimsec_volcano_PRN.pdf", width=5, height=6.5)


### stats for DEGs ###

DEG_allpair_ClusterLabel <- unique(dplyr::filter(df_deg_clusters_ref_merge_sigs, pair!="secNoni")$gene) 

df_deg_clusters_ref_merge_sigs_melsec <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::filter(pair=="melsec")
DEG_melsec_ClusterLabel <- unique(df_deg_clusters_ref_merge_sigs_melsec$gene)

df_deg_clusters_ref_merge_sigs_melsim <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::filter(pair=="melsim")
DEG_melsim_ClusterLabel <- unique(df_deg_clusters_ref_merge_sigs_melsim$gene)

df_deg_clusters_ref_merge_sigs_simsec <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::filter(pair=="simsec")
DEG_simsec_ClusterLabel <- unique(df_deg_clusters_ref_merge_sigs_simsec$gene)

#save(df_deg_clusters_ref_merge_sigs_melsec, df_deg_clusters_ref_merge_sigs_melsim, df_deg_clusters_ref_merge_sigs_simsec,
#     file="Processed_Data/DEG_cluster_ref_sigs_genelist.RData")
#load("Processed_Data/DEG_cluster_ref_sigs_genelist.RData")

length(DEG_melsec_ClusterLabel) ## melsec cluster DEGs=487 -> 593
length(DEG_melsim_ClusterLabel) ## melsim cluster DEGs=369 -> 487
length(DEG_simsec_ClusterLabel) ## simsec cluster DEGs=242 -> 335

length(intersect(intersect(DEG_melsec_ClusterLabel, DEG_melsim_ClusterLabel), DEG_simsec_ClusterLabel)) # melsimsec overlap 86 -> 140

length(DEG_melsec_ClusterLabel[!DEG_melsec_ClusterLabel %in% DEG_simsec_ClusterLabel & !DEG_melsec_ClusterLabel %in% DEG_melsim_ClusterLabel]) # melsec only - 152 -> 162
length(DEG_melsim_ClusterLabel[!DEG_melsim_ClusterLabel %in% DEG_simsec_ClusterLabel & !DEG_melsim_ClusterLabel %in% DEG_melsec_ClusterLabel]) # melsim only - 76 -> 104
length(DEG_simsec_ClusterLabel[!DEG_simsec_ClusterLabel %in% DEG_melsim_ClusterLabel & !DEG_simsec_ClusterLabel %in% DEG_melsec_ClusterLabel]) # simsec only - 38 -> 51

Dmel_spe_genes <- DEG_melsec_ClusterLabel[!DEG_melsec_ClusterLabel %in% DEG_simsec_ClusterLabel & DEG_melsec_ClusterLabel %in% DEG_melsim_ClusterLabel]
length(Dmel_spe_genes) # melsec-melsim - 169 -> 195

Dsec_spe_genes <- DEG_melsec_ClusterLabel[DEG_melsec_ClusterLabel %in% DEG_simsec_ClusterLabel & !DEG_melsec_ClusterLabel %in% DEG_melsim_ClusterLabel]
length(Dsec_spe_genes) # melsec-simsec - 80 -> 96

Dsim_spe_genes <- DEG_simsec_ClusterLabel[DEG_simsec_ClusterLabel %in% DEG_melsim_ClusterLabel & !DEG_simsec_ClusterLabel %in% DEG_melsec_ClusterLabel]
length(Dsim_spe_genes) # simsec-melsim - 38 -> 48


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


### cluster DEGs histogram ##

## group by cluster ##

df_deg_clusters_ref_merge_sigs$cluster <- factor(df_deg_clusters_ref_merge_sigs$cluster, 
                                                 levels=ClusterLabel_order$ClusterLabel)

df_deg_clusters_ref_merge_sigs %>%
  #dplyr::filter(!cluster %in% c("FMRFa(2)", "MON", "SUB", "Tbh", "Mip(2)", "Crz", "Cluster37")) %>%
  dplyr::filter(cluster %in% Anno_idents, pair != "secNoni") %>%
  ggplot(.) +
  #geom_point(shape=21) +
  #geom_line() +
  #aes(x=cluster, y=n_DEGs, fill=pair, group=pair) +
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

ggsave(file="Plots/Manuscript/Fig4_histogram_cluster_DEGs.png", width=18, height= 4)

summary_cluster_sigs_n <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::filter(cluster %in% Anno_idents, pair != "secNoni") %>%
  dplyr::group_by(pair, cluster) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::arrange(-n)

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
        panel.grid.major=element_blank()) +
  scale_fill_manual(values= c("orange", "navy",  "brown")) +
  labs(x="Differentially expressed clusters", y= "DEG count") +
  scale_x_continuous(expand=c(0.01,0.01))

ggsave(file="Plots/Manuscript/Fig4_histogram_n_cluster_DEGs.png", width=3.5, height= 1.8)

df_deg_clusters_ref_merge_sigs_nc %>%
  dplyr::filter( pair == "melsec") %>%
  dplyr::group_by(n_cluster) %>%
  dplyr::summarise(n=n())

df_deg_clusters_ref_merge_sigs_nc %>%
  dplyr::filter( pair == "melsec")%>%
  dplyr::distinct(gene) %>%
  dplyr::summarise(n=n())

melsec_deg <- df_deg_clusters_ref_merge_sigs_nc %>%
  dplyr::filter( pair == "melsec") %>%
  pull(gene)

df_per_expressed_ref_melsec_deg <- df_per_expressed_ref %>%
  dplyr::filter(species=="Dmel", gene %in% melsec_deg, ref == "OWN", pct_exp >= 0.05)  %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(fraction=n_distinct(cluster)/107) %>%
  dplyr::ungroup()

mean(df_per_expressed_ref_melsec_deg$fraction) * 107

### Dsec  specific DEGs ###

df_deg_Dsec_spec <- df_deg_clusters_ref_merge_min %>%
  dplyr::filter(gene %in% Dsec_spe_genes, pair=="simsec") %>%
  dplyr::mutate(class="Dsec_speciation")

df_deg_Dsec_spec_raw <- df_deg_clusters_ref_merge_min %>%
  dplyr::filter(gene %in% Dsec_spe_genes, pair=="simsec") %>%
  dplyr::mutate(class="Dsec_speciation")


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

ggsave("Plots/Manuscript/Fig5b.png", width=17, height=8)

df_deg_Dsec_spec_glia <- df_deg_Dsec_spec %>%
  dplyr::filter(fc_test==T & padj_test==T) %>%
  dplyr::filter(cluster %in% c("AST", "PRN", "ENS", "CTX"))

Dsec_spe_glia <- unique(df_deg_Dsec_spec_glia$gene)
length(Dsec_spe_glia) 

Dsec_spe_non_glia <- Dsec_spe_genes[!Dsec_spe_genes %in% Dsec_spe_glia]

df_deg_Dsec_spec_neuronal <- df_deg_Dsec_spec %>%
  dplyr::filter(fc_test==T & padj_test==T) %>%
  dplyr::filter(!gene %in% Dsec_spe_glia) %>%
  dplyr::filter(!cluster %in% c("Fat","Blood"))

unique(df_deg_Dsec_spec_neuronal$gene)
length(unique(df_deg_Dsec_spec_neuronal$gene))

write.table(unique(df_deg_Dsec_spec_neuronal$gene), file="Processed_Data/Dsec_spe_genes_neuronal.csv", quote=F, row.names=F, col.names=F)
write.table(Dsec_spe_glia, file="Processed_Data/Dsec_spe_genes_glia.csv", quote=F, row.names=F, col.names=F)

save(df_deg_Dsec_spec_raw, df_deg_Dsec_spec, Dsec_spe_genes, Dsec_spe_glia, Dsec_spe_non_glia, file="Processed_Data/Dsec_DEG.RData")



