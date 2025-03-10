library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_cc_cor_gather_pct_join.RData")
load("Processed_Data/df_sum_exp_SCT.RData")
load("Processed_Data/celltype_order.RData")

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

length(unique(df_sum_exp_SCT_TF$gene))

df_sum_exp_SCT_TF_top <- df_sum_exp_SCT_TF %>%
  dplyr::group_by(cluster, species) %>%
  dplyr::top_n(n = 20, wt = expression) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(cluster, gene) %>%
  dplyr::filter(n()==3) %>%
  dplyr::ungroup() 

Ter_sel <- unique(df_sum_exp_SCT_TF_top$gene)

length(Ter_sel)

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

ggsave(plot_ts, file="Plots/Manuscript/Terminal_selector.png", width=15, height= 8)


## conserved expression (r>0.6) ##

df_TF_cons_E <- df_cc_cor_gather_pct_join %>%
  dplyr::filter(value=='r', pair!="secNoni") %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(min(cor) > 0.5) %>%
  dplyr::filter(gene %in% TF_genes)

Ter_sel <- unique(df_TF_cons_E$gene)

cluster_list <- CellType_order$CellType




