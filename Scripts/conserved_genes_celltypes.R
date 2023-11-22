library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)
library(pheatmap)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/celltype_order.RData")
load("Processed_Data/df_per_expressed_ref.RData")

TrioBrain.integrated_slim_labeled_species <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled_species.rds")

df_per_expressed_freq_con <- df_per_expressed_ref %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
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
  dplyr::filter(n_celltype <= 48*0.2) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(species, cluster) %>%
  dplyr::mutate(n_gene=n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(-n_species) %>%
  dplyr::arrange(cluster, gene, species) %>%
  dplyr::rename(percent_expressed=pct_exp_max, n_high_expression_celltype=n_celltype, n_gene_per_cluster=n_gene)

write.csv(df_per_expressed_freq_con, file="conserved_specific_genes.csv", row.names = F)
  
length(unique(df_per_expressed_freq_con$gene))

df_per_expressed_freq_con %>%
  dplyr::distinct(cluster, n_gene) %>%
  View()

#### TFs ###

Dmel_gene_name <- read.table(file="Dmel/Dmel_ID_genename.tsv") %>% dplyr::distinct(V1,V2) %>%
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

write.csv(df_per_expressed_freq_con_TF, file="conserved_specific_TFs.csv", row.names = F)

df_per_expressed_freq_con_TF %>%
  dplyr::group_by(species, cluster) %>%
  dplyr::mutate(n_gene_TF=n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(cluster, n_gene_TF) %>% 
    View()

### Dotplots ###

TrioBrain.integrated_slim_labeled_Dmel <- subset(TrioBrain.integrated_slim_labeled_species, idents = anno_cluster, orig.ident!="Dmel")

TrioBrain.integrated_slim_labeled_Dmel@active.ident <- factor(TrioBrain.integrated_slim_labeled_Dmel@active.ident,
                                                                        levels= anno_cluster)
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

DotPlot_ConGenes_TF_Dmel

ggsave(file="Plots/con_genes_TF_dotplot_Dmel.png", width=12, height= 30)

df_DotPlot_ConGenes_TF_Dmel_pct <- DotPlot_ConGenes_TF_Dmel$data %>%
  dplyr::select(gene=features.plot, cluster=id, pct_exp=pct.exp) %>%
  tidyr::spread(key='gene', value='pct_exp')

df_DotPlot_ConGenes_TF_Dmel_pct$cluster <- factor(df_DotPlot_ConGenes_TF_Dmel_pct$cluster, levels= anno_cluster)

mt_DotPlot_ConGenes_TF_Dmel_pct <- df_DotPlot_ConGenes_TF_Dmel_pct %>%
  dplyr::select(-cluster) %>%
  as.matrix()

rownames(mt_DotPlot_ConGenes_TF_Dmel_pct) <- anno_cluster

plot_TFs <- pheatmap::pheatmap(mt_DotPlot_ConGenes_TF_Dmel_pct, angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                              annotation_names_col = T, show_colnames = T, cluster_rows = T,
                              treeheight_col =20, treeheight_row = 20, border_color = "black", na_col = "grey90")

ggsave(plot_TFs, file="Terminal_selectors.pdf", width=16, height= 8)




### Archive ###

df_per_expressed_ref_specific <- df_per_expressed_ref %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  dplyr::filter(pct_exp > 0.5) %>%
  dplyr::filter(ref == "OWN") %>%
  dplyr::group_by(species, gene) %>%
  dplyr::mutate(n_celltype=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_celltype <= 48*0.1) %>%
  dplyr::group_by(gene, cluster) %>%
  dplyr::mutate(n_species=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_species ==3)

df_per_expressed_ref_specific_top20 <- df_per_expressed_ref_specific %>%
  dplyr::group_by(gene, cluster) %>%
  dplyr::summarise(mean_pct_exp = mean(pct_exp)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = mean_pct_exp)

df_per_expressed_ref_specific_summary <- df_per_expressed_ref_specific %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n_congene = n()) %>%
  dplyr::ungroup()

### Dotplots ###

TrioBrain.integrated_slim_labeled_species_filter <- subset(TrioBrain.integrated_slim_labeled_species, idents = anno_cluster)

TrioBrain.integrated_slim_labeled_species_filter@active.ident <- factor(TrioBrain.integrated_slim_labeled_species_filter@active.ident,
                                                                levels= anno_cluster)

ConGenes <- unique(df_per_expressed_ref_specific_top20$gene)

DotPlot(TrioBrain.integrated_slim_labeled_species_filter, 
        feature=ConGenes, scale = T,
        scale.min=10) +
  coord_flip() + labs(x=NULL, y="Cell types") +
  theme(legend.position='bottom',
        axis.title = element_blank(),
        axis.text.x = element_text(color='black', angle=45, hjust=1, vjust=1)) +
  scale_size_continuous(trans='log2',
                        range = c(0, 3.5), breaks = c(25,50,100))
ggsave(file="Plots/con_genes_dotplot.png", width=12, height= 30)

### example ##

TrioBrain.integrated_slim_labeled_species_KC <- subset(TrioBrain.integrated_slim_labeled_species, idents = c("αβ-KC", "γ-KC","α'β'-KC"), orig.ident != "DsecNoni")

ConGenes_KC <- unique(df_per_expressed_ref_specific_KC$gene)

col3 <- sample(colors(), 3)

KC_dotplot <- DotPlot(TrioBrain.integrated_slim_labeled_species_KC, 
        feature=ConGenes_KC, split.by='orig.ident', cols=col3) +
  theme(legend.position='bottom',
        axis.title = element_blank(),
        axis.text.x = element_text(color='black', angle=90, hjust=1, vjust=0.5))

KC_dotplot

ggsave(file="Plots/KC_dotplot.png", width=10, height= 7)

KC_dotplot_data <- KC_dotplot$data %>%
  tidyr::separate(id, sep='_', into=c('cluster','species'))

df_per_expressed_ref_specific_KC <- df_per_expressed_ref_specific %>%
  dplyr::filter(grepl("KC", cluster))



df_per_expressed_ref %>%
  dplyr::filter(grepl("KC", cluster), gene %in% ConGenes_KC, ref == "OWN") %>%
  ggplot(.) +
  geom_tile(color='black', size=0.1) +
  aes(x=gene, y=species, fill=pct_exp) +
  facet_wrap(~cluster, ncol=1) +
  theme_bw()



crb_rKC_dotplot <- DotPlot(TrioBrain.integrated_slim_labeled_species_KC, idents = "γ-KC", feature=c('crb', 'cmpy'), group.by='orig.ident')
crb_rKC_dotplot_data <- crb_rKC_dotplot$data
percent_expressed_crb_rKC<- crb_rKC_dotplot_data$pct.exp ## Dmel 2.420 Dsim 3.297 Dsec 2.483 DsecNoni 2.893


DotPlot(TrioBrain.integrated_slim_labeled_species_KC, 
        feature=ConGenes_KC, scale = T, group.by='orig.ident',
        scale.min=10) +
  coord_flip() + labs(x=NULL, y="Cell types") +
  theme(legend.position='bottom',
        axis.title = element_blank(),
        axis.text.x = element_text(color='black', angle=45, hjust=1, vjust=1)) +
  scale_size_continuous(trans='log2',
                        range = c(0, 3.5), breaks = c(0,10, 25,50,75,100))
ggsave(file="Plots/con_genes_dotplot.png", width=12, height= 30)

