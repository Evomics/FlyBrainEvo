library(Seurat)
library(tidyverse)
library(ggrepel)
library(EnhancedVolcano)
library(data.table)
library(cowplot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/global_DEGs.RData")
load("Processed_Data/shared_genes.RData")
load("Processed_Data/celltype_order.RData")
load("Processed_Data/DEGs/DEG_global_refmel.RData")
load("Processed_Data/DEGs/cluster_specific_DEGs_refmel.RData")
load("Processed_Data/DEGs/cluster_specific_DEGs_refmelsim.RData")
load("Processed_Data/background_vari_features.RData")
load("Processed_Data/df_deg_sigs.RData")
load("Processed_Data/df_sum_exp_reps_SCT.RData")
load("Processed_Data/df_per_Dsec.RData")

### global DEGs ###

## merge ##

df_deg_global_melsec_ref_merge1 <- test_degs_melsec_refmel %>%
  dplyr::mutate(ref="MEL", pair='melsec')
df_deg_global_melsec_ref_merge2 <- test_degs_melsec %>%
  dplyr::mutate(ref="OWN", pair='melsec')

df_deg_global_melsim_ref_merge1 <- test_degs_melsim_refmel %>%
  dplyr::mutate(ref="MEL", pair='melsim')
df_deg_global_melsim_ref_merge2 <- test_degs_melsim %>%
  dplyr::mutate(ref="OWN", pair='melsim')

df_deg_global_simsec_ref_merge1 <- test_degs_simsec_refmel %>%
  dplyr::mutate(ref="MEL", pair='simsec')
df_deg_global_simsec_ref_merge2 <- test_degs_simsec %>%
  dplyr::mutate(ref="OWN", pair='simsec')

df_deg_global_ref_merge <- rbind(df_deg_global_melsec_ref_merge1, df_deg_global_melsec_ref_merge2,
                                 df_deg_global_melsim_ref_merge1, df_deg_global_melsim_ref_merge2,
                                 df_deg_global_simsec_ref_merge1, df_deg_global_simsec_ref_merge2)

save(df_deg_global_ref_merge, file="Processed_Data/df_deg_global_ref_merge.RData")

#load("Processed_Data/df_deg_global_ref_merge.RData")

df_deg_global_ref_merge %>%
  dplyr::filter(pair=='melsec') %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.8) +
  aes(x=avg_log2FC, y=-log10(p_val_adj), fill=ref) +
  theme_bw()

### compare ###

df_deg_global_fc_ref <- df_deg_global_ref_merge %>%
  dplyr::select(pair, gene, avg_log2FC, ref) %>%
  tidyr::spread(key='ref', value='avg_log2FC') %>%
  dplyr::rename(log2FC_MEL=MEL, log2FC_OWN=OWN) %>%
  dplyr::mutate(fc_test=ifelse(abs(log2FC_MEL)>log2(1.5) & abs(log2FC_OWN)>log2(1.5), T, F))

df_deg_global_pval_ref <- df_deg_global_ref_merge %>%
  dplyr::select(pair, gene, p_val_adj, ref) %>%
  tidyr::spread(key='ref', value='p_val_adj') %>%
  dplyr::rename(padj_MEL=MEL, padj_OWN=OWN) %>%
  dplyr::mutate(padj_test=ifelse(padj_MEL < 0.05 & padj_OWN < 0.05, T, F))

df_deg_global_fc_pval_ref <- df_deg_global_fc_ref %>%
  dplyr::left_join(., df_deg_global_pval_ref, by=c('pair', 'gene')) %>%
  dplyr::mutate(test_fc_padj=ifelse(padj_test == T & fc_test == T, T, F))

df_deg_global_fc_pval_ref %>%
  dplyr::filter(pair=='melsec') %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.8) +
  aes(x=log2FC_MEL, y=log2FC_OWN, fill=test_fc_padj) +
  theme_bw() +
  lims(x=c(-1.6,1.6), y=c(-2.6,2.6))

df_deg_global_fc_pval_ref %>%
  dplyr::filter(pair=='melsec') %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.8) +
  geom_text_repel(data=dplyr::filter(df_deg_global_fc_pval_ref, test_fc_padj ==T, pair=='melsec'), aes(label=gene)) +
  aes(x=log2FC_OWN, y=log2FC_MEL-log2FC_OWN,  fill=test_fc_padj) +
  theme_bw() +
  lims(x=c(-2.5, 2.5), y=c(-2.5, 2.5)) +
  ggtitle("mel vs sec")

df_deg_global_fc_pval_ref %>%
  dplyr::filter(pair=='melsim') %>%
  ggplot(.) +
  geom_point(shape=21, alpha=0.8) +
  geom_text_repel(data=dplyr::filter(df_deg_global_fc_pval_ref, test_fc_padj ==T, pair=='melsim'), aes(label=gene)) +
  aes(x=log2FC_OWN, y=log2FC_MEL-log2FC_OWN,  fill=test_fc_padj) +
  theme_bw()  +
  lims(x=c(-2.5, 2.5), y=c(-2.5, 2.5)) +
  ggtitle("mel vs sim")

### take min log2FC from two refs ##

df_deg_global_ref_merge_min <- df_deg_global_ref_merge %>%
  dplyr::group_by(pair, gene) %>%
  dplyr::mutate(avg_log2FC_min=min(abs(avg_log2FC)), ref_count=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(abs(avg_log2FC)==avg_log2FC_min, ref_count >1) %>%
  dplyr::mutate(fc_test = ifelse(abs(avg_log2FC) > log2(1.5),T,F), padj_test = ifelse(p_val_adj < 0.05, T,F))

df_deg_global_ref_merge_sigs <- df_deg_global_ref_merge_min %>%
  dplyr::filter(fc_test==T & padj_test==T)

df_deg_global_ref_merge_sigs %>%
  dplyr::mutate(direction = ifelse(avg_log2FC>0, "up", "down")) %>%
  dplyr::group_by(pair, direction) %>%
  dplyr::summarise(n=n())

df_deg_global_ref_merge_sigs %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(n=n())

## volcano plots ##

## mel vs sim ##

topdegs_melsim_fc <- df_deg_global_ref_merge_sigs %>%
  dplyr::filter(pair == 'melsim')

melsim_global_deg <- df_deg_global_ref_merge_min %>%
  dplyr::filter(pair == 'melsim')

EnhancedVolcano(melsim_global_deg,
                lab = melsim_global_deg$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = topdegs_melsim_fc$gene,
                drawConnectors = TRUE,
                widthConnectors = 0.25,
                lengthConnectors = unit(0.01, "npc"),
                xlim = c(-2.7, 2.7), 
                gridlines.major = T, gridlines.minor = F,
                pointSize = 1,
                labSize = 3,
                axisLabSize = 11,
                legendLabSize = 10,
                legendIconSize = 3,
                FCcutoff = log2(1.5),
                pCutoff = 0.05)

ggsave("Plots/Manuscript/FigS4_Dmelsim_volcano.pdf", width=5, height=6.5)

## mel vs sec ##

topdegs_melsec_fc <- df_deg_global_ref_merge_sigs %>%
  dplyr::filter(pair == 'melsec')

melsec_global_deg <- df_deg_global_ref_merge_min %>%
  dplyr::filter(pair == 'melsec')

EnhancedVolcano(melsec_global_deg,
                lab = melsec_global_deg$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = topdegs_melsec_fc$gene,
                drawConnectors = TRUE,
                widthConnectors = 0.25,
                lengthConnectors = unit(0.01, "npc"),
                xlim = c(-2.7, 2.7), 
                gridlines.major = T, gridlines.minor = F,
                pointSize = 1,
                labSize = 3,
                axisLabSize = 11,
                legendLabSize = 10,
                legendIconSize = 3,
                FCcutoff = log2(1.5),
                pCutoff = 0.05)

ggsave("Plots/Manuscript/FigS4_Dmelsec_volcano.pdf", width=5, height=6.5)


## sim vs sec ##

topdegs_simsec_fc <- df_deg_global_ref_merge_sigs %>%
  dplyr::filter(pair == 'simsec')

simsec_global_deg <- df_deg_global_ref_merge_min %>%
  dplyr::filter(pair == 'simsec')

EnhancedVolcano(simsec_global_deg,
                lab = simsec_global_deg$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = topdegs_simsec_fc$gene,
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
                FCcutoff = log2(1.5),
                pCutoff = 0.05)

ggsave("Plots/Manuscript/Fig4_Dsimsec_volcano.pdf", width=5, height=6.5)

 ######### cell-type specific #########

df_deg_clusters_ref_merge1 <- df_deg_clusters_ref_mel %>%
  dplyr::mutate(ref="MEL")

df_deg_clusters_ref_merge2 <- df_deg_clusters_ref_melsim %>%
  dplyr::mutate(ref="OWN")

df_deg_clusters_ref_merge <- rbind(df_deg_clusters_ref_merge1, df_deg_clusters_ref_merge2) %>%
  dplyr::mutate(cluster=ifelse(cluster=="Cluster9", "MON", cluster))

df_per_expressed_ref <- df_deg_clusters_ref_merge %>%
  dplyr::distinct(gene, cluster, pct.1, pct.2, pair, ref) %>%
  dplyr::mutate(Dmel_refMEL = ifelse(pair =="melsec" & ref=="MEL", pct.1, 
                                     ifelse(pair =="melsim" & ref=="MEL", pct.1, NA)),
                Dsim_refMEL = ifelse(pair =="melsim" & ref=="MEL", pct.2, 
                                     ifelse(pair =="simsec" & ref=="MEL", pct.1, NA)),
                Dsec_refMEL = ifelse(pair =="melsec" & ref=="MEL", pct.2, 
                                     ifelse(pair =="simsec" & ref=="MEL", pct.2, NA)),
                Dmel_refOWN = ifelse(pair =="melsec" & ref=="OWN", pct.1, 
                                     ifelse(pair =="melsim" & ref=="OWN", pct.1, NA)),
                Dsim_refOWN = ifelse(pair =="melsim" & ref=="OWN", pct.2, 
                                     ifelse(pair =="simsec" & ref=="OWN", pct.1, NA)),
                Dsec_refOWN = ifelse(pair =="melsec" & ref=="OWN", pct.2, 
                                     ifelse(pair =="simsec" & ref=="OWN", pct.2, NA)),
                
                ) %>%
  dplyr::distinct(gene, cluster, Dmel_refMEL, Dsim_refMEL, Dsec_refMEL, Dmel_refOWN, Dsim_refOWN, Dsec_refOWN) %>%
  tidyr::gather(-gene, -cluster, key="species_ref", value="pct_exp") %>%
  tidyr::separate(species_ref, into = c('species','ref'), sep='_ref') %>%
  na.omit() %>%
  dplyr::distinct()

#save(df_per_expressed_ref, file="Processed_Data/df_per_expressed_ref.RData")
load("Processed_Data/df_per_expressed_ref.RData")

### take min log2FC from two refs ##

df_deg_clusters_ref_merge_min <- df_deg_clusters_ref_merge %>%
  dplyr::group_by(pair, gene, cluster) %>%
  dplyr::mutate(avg_log2FC_min=min(abs(avg_log2FC))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(abs(avg_log2FC)==avg_log2FC_min) %>%
  dplyr::mutate(fc_test = ifelse(abs(avg_log2FC) > log2(1.5),T,F), padj_test = ifelse(p_val_adj < 0.05, T,F))

df_deg_clusters_ref_merge_sigs <- df_deg_clusters_ref_merge_min %>%
  dplyr::filter(fc_test==T & padj_test==T)

save(df_deg_global_ref_merge_min, df_deg_global_ref_merge_sigs, df_deg_clusters_ref_merge_min, df_deg_clusters_ref_merge_sigs, 
     file="Processed_Data/df_deg_sigs.RData")

#load("Processed_Data/df_deg_sigs.RData")

df_deg_clusters_ref_merge_sigs %>%
  dplyr::mutate(direction = ifelse(avg_log2FC>0, "up", "down")) %>%
  dplyr::group_by(cluster, pair, direction) %>%
  dplyr::summarise(n=n())

### Volcano plot ###

## sim vs sec - RYa ##

topdegs_simsec_fc_RYa <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::filter(pair == 'simsec', cluster == "RYa")

simsec_clusters_deg_RYa <- df_deg_clusters_ref_merge_min %>%
  dplyr::filter(pair == 'simsec', cluster == "RYa")

EnhancedVolcano(simsec_clusters_deg_RYa,
                lab = simsec_clusters_deg_RYa$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = topdegs_simsec_fc_RYa$gene,
                drawConnectors = TRUE,
                widthConnectors = 0.25,
                lengthConnectors = unit(0.01, "npc"),
                xlim = c(-2,2), 
                ylim = c(0,22), 
                gridlines.major = F, gridlines.minor = F,
                pointSize = 1,
                labSize = 3,
                axisLabSize = 11,
                legendLabSize = 10,
                legendIconSize = 3,
                borderWidth = 0.6,
                FCcutoff = log2(1.5),
                pCutoff = 0.05)

ggsave("Plots/Manuscript/Fig4_Dsimsec_volcano_RYa.pdf", width=5, height=6.5)

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

df_deg_global_ref_merge_sigs %>%
  dplyr::mutate(direction = ifelse(avg_log2FC>0, "up", "down")) %>%
  #dplyr::group_by(pair, direction) %>%
  dplyr::group_by(pair) %>%
  dplyr::summarise(n=n())

## global - melsec 14, melsim 10, simsec 4

DEG_allpair_celltype <- unique(dplyr::filter(df_deg_clusters_ref_merge_sigs, pair!="secNoni")$gene) 

df_deg_clusters_ref_merge_sigs_melsec <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::filter(pair=="melsec")
DEG_melsec_celltype <- unique(df_deg_clusters_ref_merge_sigs_melsec$gene)

df_deg_clusters_ref_merge_sigs_melsim <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::filter(pair=="melsim")
DEG_melsim_celltype <- unique(df_deg_clusters_ref_merge_sigs_melsim$gene)

df_deg_clusters_ref_merge_sigs_simsec <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::filter(pair=="simsec")
DEG_simsec_celltype <- unique(df_deg_clusters_ref_merge_sigs_simsec$gene)

save(df_deg_clusters_ref_merge_sigs_melsec, df_deg_clusters_ref_merge_sigs_melsim, df_deg_clusters_ref_merge_sigs_simsec,
     file="Processed_Data/DEG_cluster_ref_sigs_genelist.RData")
load("Processed_Data/DEG_cluster_ref_sigs_genelist.RData")

length(DEG_melsec_celltype) ## melsec cluster DEGs=487
length(DEG_melsim_celltype) ## melsim cluster DEGs=369
length(DEG_simsec_celltype) ## simsec cluster DEGs=242

length(intersect(intersect(DEG_melsec_celltype, DEG_melsim_celltype), DEG_simsec_celltype)) # melsimsec overlap 86

length(DEG_melsec_celltype[!DEG_melsec_celltype %in% DEG_simsec_celltype & !DEG_melsec_celltype %in% DEG_melsim_celltype]) # melsec only - 152
length(DEG_melsim_celltype[!DEG_melsim_celltype %in% DEG_simsec_celltype & !DEG_melsim_celltype %in% DEG_melsec_celltype]) # melsim only - 76
length(DEG_simsec_celltype[!DEG_simsec_celltype %in% DEG_melsim_celltype & !DEG_simsec_celltype %in% DEG_melsec_celltype]) # simsec only - 38

Dmel_spe_genes <- DEG_melsec_celltype[!DEG_melsec_celltype %in% DEG_simsec_celltype & DEG_melsec_celltype %in% DEG_melsim_celltype]
length(Dmel_spe_genes) # melsec-melsim - 169

Dsec_spe_genes <- DEG_melsec_celltype[DEG_melsec_celltype %in% DEG_simsec_celltype & !DEG_melsec_celltype %in% DEG_melsim_celltype]
length(Dsec_spe_genes) # melsec-simsec - 80

Dsim_spe_genes <- DEG_simsec_celltype[DEG_simsec_celltype %in% DEG_melsim_celltype & !DEG_simsec_celltype %in% DEG_melsec_celltype]
length(Dsim_spe_genes) # simsec-melsim - 38


### overlap summary

df_overlap_deg_sum <- data.frame(overlap = factor(c("melsec only", "melsim only", "simsec only",
                                                    "melsec & melsim", "melsec & simsec", "melsim & simsec", "melsec & melsim & simsec"),
                                                  levels=c("melsec only", "melsim only", "simsec only",
                                                           "melsec & melsim", "melsec & simsec", "melsim & simsec", "melsec & melsim & simsec")),
                                 count = c(152,76,38, 169, 80, 38, 86))

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

## Lineage specific gene list ##

df_deg_Dmel_spec <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::filter(gene %in% Dmel_spe_genes, pair %in% c("melsec", "melsim")) %>%
  dplyr::mutate(class="Dmel_speciation")

df_deg_Dsim_spec <- df_deg_clusters_ref_merge_sigs_simsec %>%
  dplyr::filter(gene %in% Dsim_spe_genes) %>%
  dplyr::mutate(class="Dsim_speciation")

df_deg_Dsec_spec <- df_deg_clusters_ref_merge_min %>%
  dplyr::filter(gene %in% Dsec_spe_genes, pair=="simsec") %>%
  dplyr::mutate(class="Dsec_speciation")

df_deg_Dsec_spec_raw <- df_deg_clusters_ref_merge_min %>%
  dplyr::filter(gene %in% Dsec_spe_genes, pair=="simsec") %>%
  dplyr::mutate(class="Dsec_speciation")


## Dsec lineage analysis

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
  aes(x=cluster, y=n) +
  theme_bw() +
  theme(axis.title.y = element_text(size=11, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1, color='black', size=9.5),
        axis.title.x = element_blank(),
        panel.grid.major=element_blank()) +
  scale_fill_manual(values= c("orange", "navy",  "brown"))+
  labs(y= "DEG count")

ggsave("Plots/Manuscript/Fig5b.pdf", width=7.5, height=4)

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

unique(df_deg_Dsec_spec_neuronal$cluster)
length(unique(df_deg_Dsec_spec_neuronal$gene))


save(df_deg_Dsec_spec_raw, df_deg_Dsec_spec, Dsec_spe_genes, Dsec_spe_glia, Dsec_spe_non_glia, file="Processed_Data/Dsec_DEG.RData")

df_deg_Dsec_spec$cluster <- factor(df_deg_Dsec_spec$cluster, levels= anno_cluster)

df_deg_Dsec_spec %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  ggplot(.) +
  geom_tile() +
  aes(x=cluster, y=gene, fill=avg_log2FC) +
  theme_bw() +
  scale_fill_gradient2(high='darkblue',low='darkred', mid='grey90') + 
  labs(x=NULL, y="Clusters") +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(color='black', angle=45, hjust=1, vjust=1))

mt_deg_Dsec_spec_fc <- df_deg_Dsec_spec %>%
  dplyr::select(gene, cluster, avg_log2FC) %>%
  na.omit() %>%
  distinct() %>%
  dplyr::distinct(gene, cluster, abs(avg_log2FC), .keep_all=T) %>%
  dplyr::select(-'abs(avg_log2FC)') %>%
  tidyr::spread(key='gene', value='avg_log2FC') %>%
  dplyr::select(-cluster) %>%
  dplyr::mutate(across(where(anyNA), ~ replace_na(., 0))) %>%
  as.matrix() 

rownames(mt_deg_Dsec_spec_fc) <- anno_cluster

paletteLength <- 50

myBreaks <- c(seq(min(mt_deg_Dsec_spec_fc), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(mt_deg_Dsec_spec_fc)/paletteLength, max(mt_deg_Dsec_spec_fc), length.out=floor(paletteLength/2)))

pheatmap::pheatmap(mt_deg_Dsec_spec_fc, angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                   annotation_names_col = T, show_colnames = T, cluster_rows = T,
                   treeheight_col =20, treeheight_row = 20, border_color = "black", color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
                   clustering_method="average", breaks=myBreaks, na_col = "grey90")

mt_deg_Dsec_spec_fc_scaled <- scale(mt_deg_Dsec_spec_fc, center=F)

Fig5a <- pheatmap::pheatmap(mt_deg_Dsec_spec_fc_scaled, angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                                       annotation_names_col = T, show_colnames = T, cluster_rows = T,
                                       treeheight_col =20, treeheight_row = 20, border_color = "black", color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
                   clustering_method="average", breaks=myBreaks,na_col = "grey90")

ggsave(Fig5a, file="Plots/Manuscript/Fig5a.pdf", width=17, height=9)


# pval

mt_deg_Dsec_spec_pval <- df_deg_Dsec_spec %>%
  dplyr::select(gene, cluster, p_val_adj) %>%
  dplyr::mutate(p_val_adj=-log10(p_val_adj)) %>%
  na.omit() %>%
  distinct() %>%
  dplyr::distinct(gene, cluster, p_val_adj, .keep_all=T) %>%
  dplyr::group_by(gene, cluster) %>%
  dplyr::summarise(p_val_adj=max(p_val_adj)) %>%
  dplyr::ungroup() %>%
  tidyr::spread(key='gene', value='p_val_adj') %>%
  dplyr::select(-cluster) %>%
  dplyr::mutate(across(where(anyNA), ~ replace_na(., 0))) %>%
  as.matrix() 

rownames(mt_deg_Dsec_spec_pval) <- anno_cluster

pheatmap::pheatmap(mt_deg_Dsec_spec_pval, angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                   annotation_names_col = T, show_colnames = T, cluster_rows = T,
                   treeheight_col =20, treeheight_row = 20, border_color = "black", color = colorRampPalette(c('#2471A3','white','#C0392B'))(0),
                   clustering_method="average", na_col = "grey90")

mt_deg_Dsec_spec_pval_scaled <- scale(mt_deg_Dsec_spec_pval)

pheatmap::pheatmap(mt_deg_Dsec_spec_pval_scaled, angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                   annotation_names_col = T, show_colnames = T, cluster_rows = T,
                   treeheight_col =20, treeheight_row = 20, border_color = "black", color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
                   clustering_method="average",na_col = "grey90")




df_deg_speciation <- rbind(df_deg_Dmel_spec, df_deg_Dsim_spec, df_deg_Dsec_spec) %>%
  dplyr::group_by(class) %>%
  dplyr::mutate(n_total=n()) %>%
  dplyr::group_by(class, cluster) %>%
  dplyr::summarise(n_gene=n(), freq=n_gene/n_total*100) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(class=factor(class, levels = c("Dmel_speciation", "Dsim_speciation", "Dsec_speciation"), labels=c("Dmel","Dsim","Dsec")))

df_deg_speciation %>%
  dplyr::group_by(class,cluster) %>%
  dplyr::summarise(n=n()) %>%
  View()

df_deg_speciation$cluster <- factor(df_deg_speciation$cluster, 
                                    levels=CellType_order$CellType)

df_deg_speciation %>%
  #dplyr::filter(!cluster %in% c("FMRFa(2)", "MON", "SUB", "Tbh", "Mip(2)", "Crz", "Cluster37")) %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  ggplot(.) +
  #geom_point(shape=21) +
  #geom_line() +
  #aes(x=cluster, y=n_DEGs, fill=pair, group=pair) +
  geom_col(stat='count', position = 'dodge') +
  aes(x=cluster, y=freq, fill=class) +
  theme_bw() +
  theme(axis.title.y = element_text(size=9, color='black'),
        axis.text.y = element_text(size=8, color='black'),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1, color='black', size=8),
        axis.title.x = element_blank(),
        panel.grid.major=element_blank()) +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8")) +
  labs(y= "DEG frequency (%)")


length(unique(df_deg_Dsec_spec$gene))

df_deg_Dsec_speciation <- df_deg_clusters_ref_merge_sigs_simsec %>%
  dplyr::left_join(., df_deg_clusters_ref_merge_sigs_melsec, by=c('gene','cluster')) %>%
  na.omit() %>%
  dplyr::left_join(., df_deg_clusters_ref_merge_sigs_melsim, by=c('gene','cluster')) %>%
  dplyr::filter(is.na(p_val)) %>%
  dplyr::select(gene, cluster, pair.x, p_val.x, avg_log2FC.x, pair.y, p_val.y, avg_log2FC.y) %>%
  dplyr::mutate(class="Dsec_speciation")

Dsec_spe_genes <- unique(df_deg_Dsec_speciation$gene)

length(Dsec_spe_genes) # 101 specific genes

df_deg_Dsim_speciation <- df_deg_clusters_ref_merge_sigs_simsec %>%
  dplyr::left_join(., df_deg_clusters_ref_merge_sigs_melsim, by=c('gene','cluster')) %>%
  na.omit() %>%
  dplyr::left_join(., df_deg_clusters_ref_merge_sigs_melsec, by=c('gene','cluster')) %>%
  dplyr::filter(is.na(p_val)) %>%
  dplyr::select(gene, cluster, pair.x, p_val.x, avg_log2FC.x, pair.y, p_val.y, avg_log2FC.y) %>%
  dplyr::mutate(class="Dsim_speciation")

length(Dsim_spe_genes) # 64 specific genes

df_deg_Dmel_speciation <- df_deg_clusters_ref_merge_sigs_melsec %>%
  dplyr::left_join(., df_deg_clusters_ref_merge_sigs_melsim, by=c('gene','cluster')) %>%
  na.omit() %>%
  dplyr::left_join(., df_deg_clusters_ref_merge_sigs_simsec, by=c('gene','cluster')) %>%
  dplyr::filter(is.na(p_val)) %>%
  dplyr::select(gene, cluster, pair.x, p_val.x, avg_log2FC.x, pair.y, p_val.y, avg_log2FC.y) %>%
  dplyr::mutate(class="Dmel_speciation")


### cluster DEGs histogram ##

## group by cluster ##

df_deg_clusters_ref_merge_sigs$cluster <- factor(df_deg_clusters_ref_merge_sigs$cluster, 
                                               levels=CellType_order$CellType)

df_deg_clusters_ref_merge_sigs %>%
  dplyr::filter(!cluster %in% c("FMRFa(2)", "MON", "SUB", "Tbh", "Mip(2)", "Crz", "Cluster37"), pair != "secNoni") %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
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

ggsave(file="Plots/Manuscript/Fig4_histogram_cluster_DEGs.pdf", width=9, height= 2.5)

summary_cluster_sigs_n <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::filter(!cluster %in% c("FMRFa(2)", "MON", "SUB", "Tbh", "Mip(2)", "Crz", "Cluster37"), pair != "secNoni") %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
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

ggsave(file="Plots/Manuscript/Fig4_histogram_n_cluster_DEGs.pdf", width=3.5, height= 1.8)

df_deg_clusters_ref_merge_sigs_nc %>%
  dplyr::filter( pair == "melsec") %>%
  dplyr::group_by(n_cluster) %>%
  dplyr::summarise(n=n())

df_deg_clusters_ref_merge_sigs_nc %>%
  dplyr::filter( pair == "melsec")%>%
  dplyr::distinct(gene) %>%
  dplyr::summarise(n=n())

### Differential expression - paralogs ###

select_genes <- c("ImpE1", "CG2256", "CG5149", "dob", "fiz", "MtnA", "Strn-Mlck")
trio <- c("Dmel","Dsim","Dsec")

n_select_genes <- length(select_genes)

df_per_Dsec_NA_filled <- data.frame(cluster=rep(anno_cluster, each=3*n_select_genes), gene=rep(select_genes, each=3), species=rep(trio), pct_exp2=NA) %>%
  dplyr::left_join(., dplyr::filter(df_per_Dsec, gene %in% select_genes, ref=="OWN"), by=c('gene','cluster','species')) %>%
  dplyr::mutate(pct_exp = ifelse(is.na(pct_exp), 0, pct_exp))
  
df_per_Dsec_NA_filled$cluster <- factor(df_per_Dsec_NA_filled$cluster, levels= anno_cluster)
df_per_Dsec_NA_filled$species <- factor(df_per_Dsec_NA_filled$species, levels= trio)

list_dots <- list()

for (i in 1:(n_select_genes)) {
  
  if (i==1) {
    
    list_dots[[i]] <- df_per_Dsec_NA_filled %>%
      dplyr::filter(gene == select_genes[i]) %>%
      ggplot(.) +
      geom_count() +
      aes(x=cluster, y=species, size=pct_exp*100, color=species) +
      theme_bw() +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.x= element_blank(),
            axis.text.y= element_text(color='black', size=9, face='italic'),
            axis.title.y= element_blank(),
            axis.title.x=element_blank(), 
            strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
            strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
            legend.title=element_text(color='black', size=8),
            legend.text=element_text(color='black', size=8),
            legend.key.size = unit(0.5, "lines")) +
      scale_color_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8"),guide = "none")+
      scale_size_continuous(range = c(0.1, 4), breaks = c(0,10, 20,40,60)) +
      facet_wrap(~gene, ncol=1, scales='free') +
      labs(size='Expression (%)')
    
  }
  
  else {
    
    if (i!=n_select_genes) {
      list_dots[[i]] <- df_per_Dsec_NA_filled %>%
        dplyr::filter(gene == select_genes[i]) %>%
        ggplot(.) +
        geom_count() +
        aes(x=cluster, y=species, size=pct_exp*100, color=species) +
        theme_bw() +
        theme(panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.x= element_blank(),
              axis.text.y= element_text(color='black', size=9, face='italic'),
              axis.title.y= element_blank(),
              axis.title.x=element_blank(), 
              strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
              strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
              legend.title=element_blank(),
              legend.text=element_text(color='black', size=8),
              legend.key.size = unit(0.5, "lines")) +
        scale_color_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8"),guide = "none")+
        scale_size_continuous(range = c(0.1, 4), breaks = c(0,10, 20,40,60)) +
        facet_wrap(~gene, ncol=1, scales='free') +
        labs(size='Expression (%)')
    }
    
    else {
      
      list_dots[[i]] <-df_per_Dsec_NA_filled %>%
        dplyr::filter(gene == select_genes[n_select_genes]) %>%
        ggplot(.) +
        geom_count() +
        aes(x=cluster, y=species, size=pct_exp*100, color=species) +
        theme_bw() +
        theme(panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.x= element_text(color='black', size=9, angle=45,  hjust=1, vjust=1),
              axis.text.y= element_text(color='black', size=9, face='italic'),
              axis.title.y= element_blank(),
              axis.title.x=element_blank(), 
              strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
              strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
              legend.title=element_blank(),
              legend.text=element_text(color='black', size=8),
              legend.key.size = unit(0.5, "lines")) +
        scale_color_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8"),guide = "none")+
        scale_size_continuous(range = c(0.1, 4), breaks = c(0,10, 20,40,60)) +
        facet_wrap(~gene, ncol=1, scales='free')  +
        labs(size='Expression (%)') 
      
    }
    
  }
  
  
}


plot_Dsec_degs_dots <- cowplot::plot_grid(plotlist=list_dots, ncol=1, rel_heights = c(1,1,1,1,1,1,1.4), align = 'v', axis="btlr")

ggsave(file="Plots/Manuscript/Fig5_dots.pdf", width=10, height= 7)

DotPlot(TrioBrain.integrated_slim_labeled_species, feature=c("CG9394", "CG18135"), split.by='orig.ident')

df_sum_exp_reps_SCT %>%
  dplyr::filter(gene %in% c("CG9394", "CG18135"), cluster=="AST")
  
df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="AST", gene %in% c("CG9394", "CG18135"), species != "DsecNoni") %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.1, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=9),
        axis.text.y= element_text(color='black', size=9),
        axis.title.y= element_text(color='black', size=9),
        axis.title.x=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8"))+
  facet_wrap(cluster~factor(gene, levels=c("CG9394", "CG18135")), scales = 'free', ncol=2) 

ggsave(glue::glue("Plots/Manuscript/Fig5_dotplot_AST.pdf"), width=3.5, height =3)

df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="PRN", gene %in% c("CAH2", "CAH3"), species != "DsecNoni") %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.1, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=9),
        axis.text.y= element_text(color='black', size=9),
        axis.title.y= element_text(color='black', size=9),
        axis.title.x=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8"))+
  facet_wrap(cluster~factor(gene, levels=c("CAH2", "CAH3")), scales = 'free', ncol=2) 

ggsave(glue::glue("Plots/Manuscript/Fig5_dotplot_PRN.pdf"), width=3.5, height =3)


df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="αβ-KC", gene %in% c("Hr38")) %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.1, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=9),
        axis.text.y= element_text(color='black', size=9),
        axis.title.y= element_text(color='black', size=9),
        axis.title.x=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8", "purple"))+
  facet_wrap(cluster~factor(gene, levels=c("Hr38")), scales = 'free', ncol=2) 

df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster =="PRN", gene %in% c("CG5151")) %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.1, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=9),
        axis.text.y= element_text(color='black', size=9),
        axis.title.y= element_text(color='black', size=9),
        axis.title.x=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=9),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8", "purple"))+
  facet_wrap(cluster~factor(gene, levels=c("CG5151")), scales = 'free', ncol=2) 

## MF ##

bg_go_MF_select <- bg_go_MF@result %>%
  dplyr::select(ID, Description, GeneRatio_bg=GeneRatio,BgRatio, p.adjust.bg = p.adjust) %>%
  dplyr::mutate(group="bg")

ego_MF_melsec_select <- deg_melsec_go_MF@result %>%
  dplyr::select(ID, Description, GeneRatio_pair=GeneRatio,BgRatio, p.adjust.pair = p.adjust) %>%
  dplyr::mutate(group="melsec")

ego_MF_melsim_select <- deg_melsim_go_MF@result %>%
  dplyr::select(ID, Description, GeneRatio_pair=GeneRatio,BgRatio, p.adjust.pair = p.adjust) %>%
  dplyr::mutate(group="melsim")

ego_MF_simsec_select <- deg_simsec_go_MF@result %>%
  dplyr::select(ID, Description, GeneRatio_pair=GeneRatio,BgRatio, p.adjust.pair = p.adjust) %>%
  dplyr::mutate(group="simsec")

ego_MF_pair_select_merge <- rbind(ego_MF_melsec_select, ego_MF_melsim_select, ego_MF_simsec_select)

go_MF_comp_pair <- ego_MF_pair_select_merge %>%
  dplyr::left_join(., bg_go_MF_select, by=c('ID','Description','BgRatio')) %>%
  tidyr::separate(col=GeneRatio_pair, into=c("A1","A2"),sep="/") %>%
  dplyr::mutate(GeneRatio_pair=as.numeric(A1)/as.numeric(A2))%>%
  tidyr::separate(col=GeneRatio_bg, into=c("B1","B2"),sep="/") %>%
  dplyr::mutate(GeneRatio_bg=as.numeric(B1)/as.numeric(B2)) %>%
  dplyr::mutate(ratio_diff = GeneRatio_pair/GeneRatio_bg)

go_MF_comp_pair_filtered <- go_MF_comp_pair %>%
  dplyr::filter(p.adjust.pair < 0.05 & as.numeric(A1) >= 10) %>%
  dplyr::group_by(group.x) %>%
  dplyr::top_n(n=10, wt=ratio_diff)

df_GO_enrich_MF_pair <- rbind(dplyr::mutate(deg_melsec_go_MF@result, group="melsec"), 
                              dplyr::mutate(deg_melsim_go_MF@result, group="melsim"), 
                              dplyr::mutate(deg_simsec_go_MF@result, group="simsec"), 
                              dplyr::mutate(bg_go_MF@result, group="bg")) %>%
  tidyr::separate(GeneRatio, c("N_div_class", "N_div_all")) %>%
  tidyr::separate(BgRatio, c("Total_class", "Total_genome")) %>%
  dplyr::mutate(N_div_class = as.numeric(N_div_class), N_div_all = as.numeric(N_div_all), Total_class = as.numeric(Total_class), Total_genome = as.numeric(Total_genome)) %>%
  dplyr::mutate(GeneRatio = N_div_class/N_div_all, BgRatio = Total_class/Total_genome) %>%
  dplyr::mutate(EnrichRatio = GeneRatio/BgRatio) %>%
  dplyr::arrange(-EnrichRatio) %>%
  dplyr::filter(ID %in% go_MF_comp_pair_filtered$ID)

GO_list_MF_pair <- unique(df_GO_enrich_MF_pair$Description)

GO_list_MF_pair_plotpoint <- data.frame(Description=GO_list_MF_pair, plotpoint=length(GO_list_MF_pair):1)

df_GO_enrich_MF_pair_sum <- df_GO_enrich_MF_pair %>%
  dplyr::filter(Description %in% GO_list_MF_pair) %>%
  dplyr::left_join(., GO_list_MF_pair_plotpoint, by = "Description") %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(class_gene_total = max(Total_class)) %>%
  dplyr::ungroup()

df_class_total_MF_pair <- df_GO_enrich_MF_pair_sum %>%
  dplyr::distinct(Description, class_gene_total, plotpoint) %>%
  dplyr::arrange(-plotpoint)

plot_GO_MF_pair <- df_GO_enrich_MF_pair_sum %>%
  ggplot(.) +
  geom_vline(xintercept = 1, color='blue', size=0.4) +
  geom_point(alpha=0.9, size=2) +
  aes(x=EnrichRatio, y=plotpoint, shape=group, fill=-log10(p.adjust), alpha=0.9) +
  theme_bw() +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10), expand = c(0.2,0.2)) +
  scale_y_continuous(breaks = length(GO_list_MF_pair):1, labels = GO_list_MF_pair, name = "", expand = c(0.05,0.05)) +
  #scale_y_continuous(breaks = length(GO_list_MF_pair):1, labels = GO_list_MF_pair, name = "",pair.axis = dup_axis(labels = df_class_total_MF$class_gene_total, name = "Total gene counts")) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  scale_shape_manual(values=c(21, 22,23,24)) +
  theme(axis.text = element_text(size=9, color='black'), axis.title = element_text(size=10, color='black'), 
        legend.title = element_text(size=9, color='black'), legend.text = element_text(size=8, color='black'), 
        legend.position = 'bottom', panel.grid = element_blank(),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.05, 'in')) +
  guides(shape= guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3)), 
         size= guide_legend(nrow=1, order = 2, title.position = "top"), 
         fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2)) +
  labs(x="Fold enrichment", shape = "Comparison", fill = "-log10(Corrected P-value)") 

plot_GO_MF_pair

ggsave("Plots/Manuscript/Fig4_pair_GO_MF.pdf", width=7, height=2.2)