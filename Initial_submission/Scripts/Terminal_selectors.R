library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_cc_cor_gather_pct_join.RData")
load("Processed_Data/df_sum_exp_SCT.RData")
load("Processed_Data/celltype_order.RData")
TrioBrain.integrated_slim_labeled <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled.rds")

#### TFs ###

Dmel_gene_name <- read.table(file="Dmel/Dmel_ID_genename.tsv") %>% dplyr::distinct(V1,V2) %>%
  dplyr::rename(FBID=V1, gene=V2)

df_TF <- data.table::fread("genelist/list_TF_GO0000981.txt", header = F) %>%
  dplyr::mutate(class="TF") %>%
  dplyr::select(FBID=V1, class) %>%
  dplyr::left_join(., Dmel_gene_name, by="FBID") %>% ### few genes without matching gene_name
  na.omit()

TF_genes <- df_TF$gene

## TF genes EXP summaary ##

df_sum_exp_SCT_TF <- df_sum_exp_SCT %>%
  dplyr::filter(gene %in% df_TF$gene, species != "DsecNoni") 

df_sum_exp_SCT_TF_top <- df_sum_exp_SCT_TF %>%
  dplyr::group_by(cluster, species) %>%
  dplyr::top_n(n = 20, wt = expression) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(cluster, gene) %>%
  dplyr::filter(n()==3) %>%
  dplyr::ungroup() 

Ter_sel <- unique(df_sum_exp_SCT_TF_top$gene)

df_sum_exp_SCT_TF_summary <- df_sum_exp_SCT_TF %>%
  dplyr::filter(expression>0) %>%
  dplyr::group_by(cluster, species) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() 

df_sum_exp_SCT_TF_summary %>%
  dplyr::summarise(mean=mean(n))

### normalized expression level

df_scale_ts <- df_sum_exp_SCT %>%
  dplyr::filter(species == "Dmel", gene %in% Ter_sel, cluster %in% anno_cluster) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(max_expression=max(expression)) %>%
  dplyr::ungroup() %>%
  #dplyr::filter(max_expression > 0.3) %>%
  dplyr::mutate(norm_exp = expression/max_expression) 

df_scale_ts %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  ggplot(.) +
  geom_tile(color='black', size=0.1) +
  aes(y=factor(gene, levels=rev(Ter_sel)), x=factor(cluster, levels=anno_cluster), fill=norm_exp) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y= element_text(color='black', face='italic', size=9),
        axis.text.x = element_text(color='black', size=9, angle=45, hjust = 1, vjust = 1),
        legend.title = element_text(size=9.5, color='black'),
        legend.text = element_text(size=9, color = 'black'),
        axis.title=element_blank(),
        legend.position='bottom',
        strip.text.x = element_text(margin = margin(0.02,0,0.02,0, "in"))) +
  scale_fill_gradient(low="white", high="darkblue", breaks = c(0, 0.2,0.4,0.6,0.8,1)) +
  #scale_fill_manual(values=c("white","darkblue")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="expression")

#ggsave(file="Plots/Manuscript/Fig2_heatmap_ts_gene_cluster_exp.pdf", width=10, height= 8)

df_scale_ts$cluster <- factor(df_scale_ts$cluster, levels= anno_cluster)

mt_scale_ts <- df_scale_ts %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  dplyr::select(gene, cluster, norm_exp) %>%
  tidyr::spread(key='gene', value='norm_exp') %>%
  dplyr::select(-cluster) %>%
  as.matrix()

rownames(mt_scale_ts) <- anno_cluster

plot_ts <- pheatmap::pheatmap(mt_scale_ts, angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                   annotation_names_col = T, show_colnames = T, cluster_rows = T,
                   treeheight_col =20, treeheight_row = 20, border_color = "black", na_col = "grey90")

ggsave(plot_ts, file="Plots/Manuscript/working/230403_Extended_Data_Fig2d_ver4.pdf", width=16, height= 8)


## Example - Sp1 and Rx ###

Plot_Rx_Sp1 <- df_sum_exp_SCT_TF %>%
  dplyr::filter(cluster %in% anno_cluster, gene %in% c("Sp1","Rx")) %>%
  dplyr::group_by(species, cluster, gene) %>%
  dplyr::mutate(mean_exp = mean(expression)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(gene, species) %>%
  dplyr::mutate(max_expression=max(mean_exp)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec"))) %>%
  #dplyr::filter(max_expression > 0.3) %>%
  dplyr::mutate(norm_exp = mean_exp/max_expression) %>%
  ggplot(.) +
  geom_tile(color='black', size=0.1) +
  aes(y=species, x=factor(cluster, levels=anno_cluster), fill=norm_exp) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y= element_text(color='black', face='italic', size=9),
        axis.text.x = element_text(color='black', size=9, angle=45, hjust = 1, vjust = 1),
        legend.title = element_text(size=9.5, color='black'),
        legend.text = element_text(size=9, color = 'black'),
        axis.title=element_blank(),
        legend.position='bottom',
        strip.text.x = element_text(margin = margin(0.02,0,0.02,0, "in"))) +
  scale_fill_gradient(low="white", high="darkblue", breaks = c(0, 0.2,0.4,0.6,0.8,1)) +
  #scale_fill_manual(values=c("white","darkblue")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="expression") +
  facet_wrap(~gene, ncol=1)

ggsave(Plot_Rx_Sp1, file="Plots/Plot_Rx_Sp1.pdf", width=10, height= 2.5)

## range of expression ##

load("Processed_Data/df_per_expressed_ref.RData")

## gene centric ##

df_per_expressed_ref_TF1 <- df_per_expressed_ref %>%
  dplyr::filter(gene %in% TF_genes, pct_exp > 0.1) %>%
  dplyr::filter(ref == "OWN") %>%
  dplyr::group_by(species, gene) %>%
  dplyr::summarise(n_celltype=n()) %>%
  dplyr::ungroup()

df_per_expressed_ref_TF1 %>%
  ggplot(.) +
  geom_histogram() +
  aes(x=n_celltype, fill=species) +
  theme_bw() +
  facet_grid(~species)

## cell type centric ##

df_per_expressed_ref_TF2 <- df_per_expressed_ref %>%
  dplyr::filter(gene %in% TF_genes, pct_exp > 0.1) %>%
  dplyr::filter(ref == "OWN") %>%
  dplyr::group_by(species, cluster) %>%
  dplyr::summarise(n_gene=n()) %>%
  dplyr::ungroup()

df_per_expressed_ref_TF2 %>%
  ggplot(.) +
  geom_histogram(binwidth=5) +
  aes(x=n_gene, fill=species) +
  theme_bw() +
  facet_grid(~species)

df_per_expressed_ref_TF1_specific <- df_per_expressed_ref %>%
  dplyr::filter(gene %in% TF_genes, pct_exp > 0.2) %>%
  dplyr::filter(ref == "OWN") %>%
  dplyr::group_by(species, gene) %>%
  dplyr::mutate(n_celltype=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_celltype <= 64*0.2) %>%
  dplyr::group_by(gene, cluster) %>%
  dplyr::summarise(n_species=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_species ==3)


### Dotplots ###

TrioBrain.integrated_slim_labeled_filter <- subset(TrioBrain.integrated_slim_labeled, idents = anno_cluster)

TrioBrain.integrated_slim_labeled_filter@active.ident <- factor(TrioBrain.integrated_slim_labeled_filter@active.ident,
                                                                levels= anno_cluster)

Ter_cans <- unique(df_per_expressed_ref_TF1_specific$gene)

DotPlot(TrioBrain.integrated_slim_labeled_filter, 
        feature=Ter_cans, scale = T,
        scale.min=10) +
  coord_flip() + labs(x=NULL, y="Cell types") +
  theme(legend.position='bottom',
        axis.title = element_blank(),
        axis.text.x = element_text(color='black', angle=45, hjust=1, vjust=1)) +
  scale_size_continuous(trans='log2',
                        range = c(0, 3.5), breaks = c(25,50,100))



df_per_expressed_ref_TF2 %>%
  ggplot(.) +
  geom_histogram(binwidth=5) +
  aes(x=n_gene, fill=species) +
  theme_bw() +
  facet_grid(~species)


## conserved expression (r>0.6) ##

df_TF_cons_E <- df_cc_cor_gather_pct_join %>%
  dplyr::filter(value=='r', pair!="secNoni") %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(min(cor) > 0.5) %>%
  dplyr::filter(gene %in% TF_genes)

Ter_sel <- unique(df_TF_cons_E$gene)

cluster_list <- CellType_order$CellType




