library(Seurat)
library(tidyverse)
library(cowplot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/celltype_order.RData")
load("Processed_Data/Dsec_DEG.RData")
load("Processed_Data/df_per_expressed_ref.RData")
load("Processed_Data/df_per_Dsec.RData")

df_deg_Dsec_spec_raw_summary <- df_deg_Dsec_spec_raw %>%
  dplyr::filter(padj_test == T, fc_test == T)  %>%
  dplyr::select(gene, cluster, p_val_adj, avg_log2FC, pct_exp_Dsim=pct.1, pct_exp_Dsec=pct.2) %>%
  dplyr::arrange(gene, p_val_adj)

write.csv(df_deg_Dsec_spec_raw_summary, file="Dsec_spe_gene_info.csv", row.names = F)

length(unique(df_deg_Dsec_spec_raw_summary$gene)) 

### Dsec_DEG info ###

TrioBrain.integrated_slim_labeled_species <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled_species.rds")

TrioBrain.integrated_slim_labeled_trio <- subset(TrioBrain.integrated_slim_labeled_species, orig.ident != "DsecNoni")

trio_col <- c("#E41A1C", "#4DAF4A",  "#377EB8")

amineR <- c('5-HT1A', 'Octbeta2R', 'Octbeta3R')
neuro_etc <- c('amon', 'Pde1c', 'Glut1')
glia_meta <- c('Cat', 'dob', 'CG3961', 'Sarm', 'sbm', 'CG4991') 
para1 <- c('CG9394', 'CG18135')
para2 <- c('CAH2', 'CAH3')
ecdy <- c('tai', 'ImpE1', 'CG10513', 'fiz')

## dot plot for 14 genes referred in the MS ##

### Differential expression - paralogs ###

select_genes <- c(glia_meta, para1, para2, ecdy)

Dotplot_expression_data <- DotPlot(TrioBrain.integrated_slim_labeled_trio, feature=select_genes, split.by='orig.ident', cols=trio_col)$data %>%
  tidyr::separate(id, into=c('cluster','species'), sep='_') %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  dplyr::select(species, cluster, gene=features.plot, avg.exp, pct.exp)

Dotplot_expression_data$species <- factor(Dotplot_expression_data$species, levels=c('Dmel', 'Dsim', 'Dsec'))
Dotplot_expression_data$cluster <- factor(Dotplot_expression_data$cluster, levels=anno_cluster)

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

plot_Dsec_degs_dots <- cowplot::plot_grid(plotlist=list_dots, ncol=2, align = 'v', axis="btlr")

plot_Dsec_degs_dots

ggsave(file="Plots/Fig5_dots_raw.pdf", width=22, height= 11)
  
trio <- c("Dmel","Dsim","Dsec")

n_select_genes <- length(select_genes)

df_per_Dsec_NA_filled <- data.frame(cluster=rep(anno_cluster, each=3*n_select_genes), gene=rep(select_genes, each=3), species=rep(trio), pct_exp2=NA) %>%
  dplyr::left_join(., dplyr::filter(df_per_Dsec, gene %in% select_genes, ref=="OWN"), by=c('gene','cluster','species')) %>%
  dplyr::mutate(pct_exp = ifelse(is.na(pct_exp), 0, pct_exp))

df_per_Dsec_NA_filled$cluster <- factor(df_per_Dsec_NA_filled$cluster, levels= anno_cluster)
df_per_Dsec_NA_filled$species <- factor(df_per_Dsec_NA_filled$species, levels= trio)






DotPlot(TrioBrain.integrated_slim_labeled_trio, feature=glia_meta)

DotPlot(TrioBrain.integrated_slim_labeled_trio, feature=c('Cat','dob','Sarm'))

### Glial Dsec DEGs ###

DotPlot(subset(TrioBrain.integrated_slim_labeled_trio, idents = c('CTX')), feature=c('Cat','dob','Sarm'), group.by='orig.ident') + 
  theme(title = element_text(size=11, color='black'),
        axis.text=element_text(size=11, face='italic', color='black'),
        axis.title=element_blank(),
        legend.title=element_text(size=10, color='black'),
        legend.text=element_text(size=9, color='black'),
        legend.key.width = unit(0.2, 'in'),
        legend.key.height = unit(0.08, 'in'),
        legend.direction = 'horizontal',
        legend.position = 'bottom',
        axis.line=element_line(size=0.2),
        axis.ticks = element_line(size=0.2)) +
  ggtitle("CTX")

ggsave(file="Plots/CTX_meta.pdf", width=6, height=2)

DotPlot(subset(TrioBrain.integrated_slim_labeled_trio, idents = c('PRN')), feature=c('sbm','CG3961','CG4991', 'CAH2', 'CAH3'), group.by='orig.ident') + 
  theme(title = element_text(size=11, color='black'),
        axis.text=element_text(size=11, face='italic', color='black'),
        axis.title=element_blank(),
        legend.title=element_text(size=10, color='black'),
        legend.text=element_text(size=9, color='black'),
        legend.key.width = unit(0.2, 'in'),
        legend.key.height = unit(0.08, 'in'),
        legend.direction = 'horizontal',
        legend.position = 'bottom',
        axis.line=element_line(size=0.2),
        axis.ticks = element_line(size=0.2)) +
  ggtitle("PRN")

ggsave(file="Plots/PRN_meta.pdf", width=6, height=2)

DotPlot(subset(TrioBrain.integrated_slim_labeled_trio, idents = c('AST')), feature=c('CG9394', 'CG18135', "tai", "CG10513"), group.by='orig.ident') + 
  theme(title = element_text(size=11, color='black'),
        axis.text=element_text(size=11, face='italic', color='black'),
        axis.title=element_blank(),
        legend.title=element_text(size=10, color='black'),
        legend.text=element_text(size=9, color='black'),
        legend.key.width = unit(0.2, 'in'),
        legend.key.height = unit(0.08, 'in'),
        legend.direction = 'horizontal',
        legend.position = 'bottom',
        axis.line=element_line(size=0.2),
        axis.ticks = element_line(size=0.2)) +
  ggtitle("AST")

ggsave(file="Plots/AST_meta'.pdf", width=6, height=2)

### Ecdysone signaling ###

DotPlot(TrioBrain.integrated_slim_labeled_trio, feature=c('tai','ImpE1','CG10513'))

DotPlot(subset(TrioBrain.integrated_slim_labeled_trio, idents= c("AST", "CTX", "ENS", "PRN", "SUB", "DOP", "GABA")),
               feature=c('tai','ImpE1','CG10513'), split.by='orig.ident', cols=trio_col)

Dotplot_tai_data <- DotPlot(TrioBrain.integrated_slim_labeled_trio, feature=c('tai','ImpE1','CG10513'), split.by='orig.ident', cols=trio_col)$data %>%
  tidyr::separate(id, into=c('cluster','species'), sep='_') %>%
  dplyr::filter(cluster %in% anno_cluster)

Dotplot_tai_data$species <- factor(Dotplot_tai_data$species, levels=c('Dmel', 'Dsim', 'Dsec'))
Dotplot_tai_data$cluster <- factor(Dotplot_tai_data$cluster, levels=anno_cluster)

Dotplot_tai_data %>%
  dplyr::filter(features.plot == 'tai') %>%
  ggplot(.) +
  geom_point(shape=21) +
  aes(x=cluster, y=species, size=pct.exp, fill=avg.exp) +
  theme_bw() +
  facet_wrap(~features.plot, ncol=1, scales = 'free')
  
## Supplemetnary Fig -  ImpE1 expression pattern ##

Plot_Dmel_ImpE1 <- DotPlot(subset(TrioBrain.integrated_slim_labeled_trio, orig.ident == "Dmel", idents = anno_cluster), feature="ImpE1") + scale_radius(limits=c(0,50), range=c(0,6)) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size=11, color='black'),
        legend.text = element_text(size=10, color='black'),
        legend.title = element_text(size=11, color='black'),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.2, 'in'),
        legend.direction = 'vertical',
        legend.position = 'right')
Plot_Dsim_ImpE1 <- DotPlot(subset(TrioBrain.integrated_slim_labeled_trio, orig.ident == "Dsim", idents = anno_cluster), feature="ImpE1") + scale_radius(limits=c(0,50), range=c(0,6))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=11, color='black'),
        legend.text = element_text(size=10, color='black'),
        legend.title = element_text(size=11, color='black'),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.2, 'in'),
        legend.direction = 'vertical',
        legend.position = 'right')

Plot_Dsec_ImpE1 <- DotPlot(subset(TrioBrain.integrated_slim_labeled_trio, orig.ident == "Dsec", idents = anno_cluster), feature="ImpE1") + scale_radius(limits=c(0,50), range=c(0,6))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=11, color='black'),
        legend.text = element_text(size=10, color='black'),
        legend.title = element_text(size=11, color='black'),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.2, 'in'),
        legend.direction = 'vertical',
        legend.position = 'right')

Plot_ImpE1 <- cowplot::plot_grid(Plot_Dmel_ImpE1, Plot_Dsim_ImpE1, Plot_Dsec_ImpE1, nrow=1)               

Plot_ImpE1               

ggsave(file="Plots/ImpE1.pdf", width=12, height=8)

## Supplemetnary Fig -  fiz expression pattern ##

Plot_Dmel_fiz <- DotPlot(subset(TrioBrain.integrated_slim_labeled_trio, orig.ident == "Dmel", idents = anno_cluster), feature="fiz") + scale_radius(limits=c(0,70), range=c(0,6)) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size=11, color='black'),
        legend.text = element_text(size=10, color='black'),
        legend.title = element_text(size=11, color='black'),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.2, 'in'),
        legend.direction = 'vertical',
        legend.position = 'right')
Plot_Dsim_fiz <- DotPlot(subset(TrioBrain.integrated_slim_labeled_trio, orig.ident == "Dsim", idents = anno_cluster), feature="fiz") + scale_radius(limits=c(0,70), range=c(0,6))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=11, color='black'),
        legend.text = element_text(size=10, color='black'),
        legend.title = element_text(size=11, color='black'),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.2, 'in'),
        legend.direction = 'vertical',
        legend.position = 'right')

Plot_Dsec_fiz <- DotPlot(subset(TrioBrain.integrated_slim_labeled_trio, orig.ident == "Dsec", idents = anno_cluster), feature="fiz") + scale_radius(limits=c(0,70), range=c(0,6))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=11, color='black'),
        legend.text = element_text(size=10, color='black'),
        legend.title = element_text(size=11, color='black'),
        legend.key.width = unit(0.1, 'in'),
        legend.key.height = unit(0.2, 'in'),
        legend.direction = 'vertical',
        legend.position = 'right')

Plot_fiz <- cowplot::plot_grid(Plot_Dmel_fiz, Plot_Dsim_fiz, Plot_Dsec_fiz, nrow=1)               

Plot_fiz               

ggsave(file="Plots/fiz.pdf", width=12, height=8)

