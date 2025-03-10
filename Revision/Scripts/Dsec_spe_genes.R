library(Seurat)
library(tidyverse)
library(cowplot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/celltype_order.RData")
load("Processed_Data/Dsec_DEG.RData")
load("Processed_Data/df_per_expressed_ref.RData")
load("Processed_Data/df_per_Dsec.RData")

TrioBrain.integrated_slim_labeled_final <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_labeled_final.rds")
TrioBrain.integrated_slim_labeled_final_species <- TrioBrain.integrated_slim_labeled_final
TrioBrain.integrated_slim_labeled_final_species$orig.ident <- gsub("_rep[1-6]", "", TrioBrain.integrated_slim_labeled_final_species$orig.ident)
TrioBrain.integrated_slim_labeled_final_species@meta.data$orig.ident <- factor(TrioBrain.integrated_slim_labeled_final_species@meta.data$orig.ident, 
                                                                               levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))

df_deg_Dsec_spec_raw_summary <- df_deg_Dsec_spec_raw %>%
  dplyr::filter(padj_test == T, fc_test == T)  %>%
  dplyr::select(gene, cluster, p_val_adj, avg_log2FC, pct_exp_Dsim=pct.1, pct_exp_Dsec=pct.2) %>%
  dplyr::arrange(gene, p_val_adj)

#write.csv(df_deg_Dsec_spec_raw_summary, file="Processed_Data/Dsec_spe_gene_info.csv", row.names = F)

length(unique(df_deg_Dsec_spec_raw_summary$gene)) 

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



Dotplot_expression_data$species <- factor(Dotplot_expression_data$species, levels=c('Dmel', 'Dsim', 'Dsec'))
Dotplot_expression_data$cluster <- factor(Dotplot_expression_data$cluster, levels=Anno_idents)

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
            axis.text.x= element_text(color='black', size=9, angle=45,  hjust=1, vjust=1),
            axis.text.y= element_text(color='black', size=9, face='italic'),
            axis.title.y= element_blank(),
            axis.title.x=element_blank(), 
            strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), face='italic', size=10),
            strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=9),
            strip.background = element_blank(),
            legend.title=element_text(color='black', size=8),
            legend.text=element_text(color='black', size=8),
            legend.key.size = unit(0.5, "lines"),
            legend.box = "horizontal") +
      scale_color_gradient(low = "lightgrey", high="blue")+
      scale_size_continuous(range = c(0.1, 4), breaks = c(0,10, 20,40,60,80)) +
      facet_wrap(~gene, ncol=1, scales='free') +
      labs(size='Frequency (%)', color='Expression')
    
    
}

plot_Dsec_degs_dots <- cowplot::plot_grid(plotlist=list_dots, ncol=1, align = 'v', axis="btlr")

plot_Dsec_degs_dots

ggsave(file="Plots/Manuscript/Fig5_dots_raw.pdf", width=22, height= 20)
