library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)
library(pheatmap)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/celltype_order.RData")
load("Processed_Data/df_per_expressed_ref.RData")

TrioBrain.integrated_slim_labeled_final <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_labeled_final.rds")
TrioBrain.integrated_slim_labeled_final_species <- TrioBrain.integrated_slim_labeled_final
TrioBrain.integrated_slim_labeled_final_species$orig.ident <- gsub("_rep[1-6]", "", TrioBrain.integrated_slim_labeled_final_species$orig.ident)

df_per_expressed_freq_con <- df_per_expressed_ref %>%
  dplyr::filter(cluster %in% Anno_idents) %>%
  dplyr::group_by(gene, cluster, species) %>%
  dplyr::summarise(pct_exp_max=max(pct_exp)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(pct_exp_max >= 0.3) %>%
  dplyr::group_by(gene, cluster) %>%
  dplyr::mutate(n_species=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_species ==3) %>%
  dplyr::group_by(species, gene) %>%
  dplyr::mutate(n_celltype=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_celltype <= 107*0.2) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(species, cluster) %>%
  dplyr::mutate(n_gene=n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(-n_species) %>%
  dplyr::arrange(cluster, gene, species) %>%
  dplyr::rename(percent_expressed=pct_exp_max, n_high_expression_celltype=n_celltype, n_gene_per_cluster=n_gene)

write.csv(df_per_expressed_freq_con, file="Processed_Data/conserved_specific_genes.csv", row.names = F, quote= F)
  
length(unique(df_per_expressed_freq_con$gene))

df_per_expressed_freq_con %>%
  dplyr::distinct(cluster, n_gene_per_cluster) %>%
  View()

#### TFs ###

Dmel_gene_name <- read.table(file="genelist/Dmel_ID_genename.tsv") %>% dplyr::distinct(V1,V2) %>%
  dplyr::rename(FBID=V1, gene=V2)

df_TF <- data.table::fread("genelist/list_TF_GO0000981.txt", header = F) %>%
  dplyr::mutate(class="TF") %>%
  dplyr::select(FBID=V1, class) %>%
  dplyr::left_join(., Dmel_gene_name, by="FBID") %>% ### few genes without matching gene_name
  na.omit()

TF_genes <- df_TF$gene

df_per_expressed_freq_con_TF <- df_per_expressed_freq_con %>%
  dplyr::filter(gene %in% TF_genes)

length(unique(df_per_expressed_freq_con_TF$gene))

write.csv(df_per_expressed_freq_con_TF, file="Processed_Data/conserved_specific_TFs.csv", row.names = F)

df_per_expressed_freq_con_TF %>%
  dplyr::group_by(species, cluster) %>%
  dplyr::mutate(n_gene_TF=n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(cluster, n_gene_TF) %>% 
    View()

### Dotplots ###

TrioBrain.integrated_slim_labeled_Dmel <- subset(TrioBrain.integrated_slim_labeled_final_species, idents = Anno_idents, orig.ident!="Dmel")

TrioBrain.integrated_slim_labeled_Dmel@active.ident <- factor(TrioBrain.integrated_slim_labeled_Dmel@active.ident,
                                                                        levels= Anno_idents)
ConGenes_TF <- unique(df_per_expressed_freq_con_TF$gene)

DotPlot_ConGenes_TF_Dmel <- DotPlot(TrioBrain.integrated_slim_labeled_Dmel, 
        feature=ConGenes_TF, scale = T,
        scale.min=0) +
  coord_flip() + labs(x=NULL, y="Cell types") +
  theme(legend.position='bottom',
        axis.title = element_blank(),
        axis.text.x = element_text(color='black', angle=45, hjust=1, vjust=1)) +
  scale_size_continuous(trans='identity',
                        range = c(0, 4), breaks = c(25,50,100))

df_DotPlot_ConGenes_TF_Dmel_pct <- DotPlot_ConGenes_TF_Dmel$data %>%
  dplyr::select(gene=features.plot, cluster=id, pct_exp=pct.exp) %>%
  tidyr::spread(key='gene', value='pct_exp')

df_DotPlot_ConGenes_TF_Dmel_pct$cluster <- factor(df_DotPlot_ConGenes_TF_Dmel_pct$cluster, levels= Anno_idents)

mt_DotPlot_ConGenes_TF_Dmel_pct <- df_DotPlot_ConGenes_TF_Dmel_pct %>%
  dplyr::select(-cluster) %>%
  as.matrix()

rownames(mt_DotPlot_ConGenes_TF_Dmel_pct) <- Anno_idents

plot_TFs <- pheatmap::pheatmap(mt_DotPlot_ConGenes_TF_Dmel_pct, angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                              annotation_names_col = T, show_colnames = T, cluster_rows = T,
                              treeheight_col =20, treeheight_row = 20, border_color = "black", na_col = "grey90")

ggsave(plot_TFs, file="Plots/Manuscript/Terminal_selectors.pdf", width=20, height= 15)




