library(tidyverse)
library(Seurat)
library(ggrepel)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/Dsec_DEG.RData")
load("Processed_Data/df_sum_exp_reps_SCT.RData")
load("Processed_Data/df_sum_exp_SCT.RData")
load("Processed_Data/df_per_expressed_ref.RData")

cluster_list <- CellType_order$CellType

for (gen in Dsec_spe_genes) {
  
  df_sum_exp_SCT %>%
    dplyr::filter(gene == gen, species != "DsecNoni", cluster %in% anno_cluster) %>%
    ggplot(.) +
    geom_tile(color='black') +
    aes(x=factor(cluster, levels=anno_cluster), y=factor(species, levels = c("Dsec", "Dsim", "Dmel")), fill=expression) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.y= element_text(color='black', face='italic', size=9),
          axis.text.x = element_text(color='black', size=9, angle=90, hjust = 1, vjust = 0.5),
          legend.title = element_text(size=9.5, color='black'),
          legend.text = element_text(size=9, color = 'black'),
          axis.title=element_blank(),
          legend.position='bottom',
          strip.text.x = element_text(margin = margin(0.02,0,0.02,0, "in"))) +
    scale_fill_gradient(low="white", high="blue") +
    facet_wrap(~gene, ncol=1, scales='free') +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))
  
  ggsave(glue::glue("Plots/DsecDEG/heatplot/{gen}_annotated.png"), width=5, height = 2)
  
  df_sum_exp_SCT %>%
    dplyr::filter(gene == gen, species != "DsecNoni") %>%
    ggplot(.) +
    geom_tile(color='black') +
    aes(x=cluster, y=factor(species, levels = c("Dsec", "Dsim", "Dmel")), fill=expression) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.y= element_text(color='black', face='italic', size=9),
          axis.text.x = element_text(color='black', size=9, angle=90, hjust = 1, vjust = 0.5),
          legend.title = element_text(size=9.5, color='black'),
          legend.text = element_text(size=9, color = 'black'),
          axis.title=element_blank(),
          legend.position='bottom',
          strip.text.x = element_text(margin = margin(0.02,0,0.02,0, "in"))) +
    scale_fill_gradient(low="white", high="blue") +
    facet_wrap(~gene, ncol=1, scales='free') +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))
  
  ggsave(glue::glue("Plots/DsecDEG/heatplot/{gen}_allcluster.png"), width=10, height = 2)
  
  
}

## scatter plots ##

df_DsecDEG_cor <- NULL

for (gen in Dsec_spe_genes) {
  
  df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
    dplyr::filter(gene==gen) %>%
    tidyr::spread(key='species', value='expression')
  
  DsecDEG_cor_melsec <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec)
  DsecDEG_cor_melsim <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsim)
  DsecDEG_cor_simsec <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$Dsim)
  DsecDEG_cor_secNoni <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$DsecNoni)
  
  df_DsecDEG_cor_temp <- data.frame(gene=gen, 
                               melsec=DsecDEG_cor_melsec,
                               melsim=DsecDEG_cor_melsim, 
                               simsec=DsecDEG_cor_simsec, 
                               secNoni=DsecDEG_cor_secNoni)
  
  df_DsecDEG_cor <- rbind(df_DsecDEG_cor, df_DsecDEG_cor_temp)
  
  df_sum_exp_SCT_select %>%
    ggplot(.) +
    geom_smooth(size= 0.5, method="glm") +
    geom_point(alpha = 0.8, size=1) +
    geom_label_repel(aes(label=cluster), label.size=0.1, segment.size = 0.1, box.padding = 0.1,
                     label.padding = 0.1, size = 3) +
    aes(x=Dmel, y=Dsec) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color='black', size=10),
          axis.title= element_text(color='black', size=11, face='italic')) +
    ggtitle(paste(gen, "Dmel vs Dsec | r =", round(cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec),3)))
  
  ggsave(glue::glue("Plots/DsecDEG/scatterplot/{gen}_allcluster_scatter_melsec.pdf"), width=5, height =5)
  
  
  df_sum_exp_SCT_select %>%
    ggplot(.) +
    geom_smooth(size= 0.5, method="glm") +
    geom_point(alpha = 0.8, size=1) +
    geom_label_repel(aes(label=cluster), label.size=0.1, segment.size = 0.1, box.padding = 0.1,
                     label.padding = 0.1, size = 3) +
    aes(x=Dmel, y=Dsim) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color='black', size=10),
          axis.title= element_text(color='black', size=11, face='italic')) +
    ggtitle(paste(gen, "Dmel vs Dsim | r =", round(cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsim),3)))
  
  ggsave(glue::glue("Plots/DsecDEG/scatterplot/{gen}_allcluster_scatter_melsim.pdf"), width=5, height = 5)
  
  df_sum_exp_SCT_select %>%
    ggplot(.) +
    geom_smooth(size= 0.5, method="glm") +
    geom_point(alpha = 0.8, size=1) +
    geom_label_repel(aes(label=cluster), label.size=0.1, segment.size = 0.1, box.padding = 0.1,
                     label.padding = 0.1, size = 3) +
    aes(x=Dsim, y=Dsec) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color='black', size=10),
          axis.title= element_text(color='black', size=11, face='italic')) +
    ggtitle(paste(gen, "Dsim vs Dsec | r =", round(cor(df_sum_exp_SCT_select$Dsim, df_sum_exp_SCT_select$Dsec),3)))
  
  ggsave(glue::glue("Plots/DsecDEG/scatterplot/{gen}_allcluster_scatter_simsec.pdf"), width=5, height =5)
  
  df_sum_exp_SCT_select %>%
    ggplot(.) +
    geom_smooth(size= 0.5, method="glm") +
    geom_point(alpha = 0.8, size=1) +
    geom_label_repel(aes(label=cluster), label.size=0.1, segment.size = 0.1, box.padding = 0.1,
                     label.padding = 0.1, size = 3) +
    aes(x=Dsec, y=DsecNoni) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color='black', size=10),
          axis.title= element_text(color='black', size=11, face='italic')) +
    ggtitle(paste(gen, "Dsec vs DsecNoni | r =", round(cor(df_sum_exp_SCT_select$DsecNoni, df_sum_exp_SCT_select$Dsec),3)))
  
  ggsave(glue::glue("Plots/DsecDEG/scatterplot/{gen}_allcluster_scatter_secNoni.pdf"), width=5, height = 5)
  
}

df_DsecDEG_cor_gather <- df_DsecDEG_cor %>%
  tidyr::gather(-gene, key='pair', value='cor')



df_sum_exp_SCT_select %>%
  dplyr::filter(gene %in% Dsec_spe_genes) %>%
  ggplot(.) +
  geom_smooth(size= 0.3, method="glm") +
  geom_point(alpha = 0.6, size=1.2) +
  geom_text_repel(data=dplyr::filter(df_sum_exp_SCT_select, 
                                     cluster %in% anno_cluster),
                  #                                    cluster %in% select_cluster),
                  aes(label=cluster), label.size=0.1, segment.size = 0.1, box.padding = 0.1,nudge_x = 0.1,
                  label.padding = 0.1, size = 3) +
  aes(x=Dmel, y=Dsec) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title= element_text(color='black', size=11, face='italic')) +
  ggtitle(paste(gen, "| rho =", round(cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method='spearman'),3)))

## Dynamic expression heatmap by pct exp ##

df_per_Dsec <- df_per_expressed_ref %>%
  dplyr::filter(gene %in% Dsec_spe_genes)

df_per_Dsec$cluster <- factor(df_per_Dsec$cluster, 
                                    levels=CellType_order$CellType)

df_per_Dsec$species <- factor(df_per_Dsec$species, 
                              levels=c("Dmel","Dsim","Dsec"))

df_per_Dsec$cluster <- factor(df_per_Dsec$cluster, 
                              levels=anno_cluster)

df_per_Dsec %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  dplyr::filter(gene == "Strn-Mlck") %>%
  ggplot(.) +
  geom_tile() +
  aes(x=species, y=cluster, fill=pct_exp*100) +
  scale_fill_gradient(high='darkblue',low='white') + 
  theme_bw() +
  facet_grid(~ref)


### gene vs gene correlation ###

df_sum_exp_SCT_Dsec_spe <- df_sum_exp_SCT %>%
  dplyr::filter(gene %in% Dsec_spe_genes, species == "Dsec") %>%
  dplyr::select(gene, expression, cluster) %>%
  tidyr::spread(key='gene', value='expression') %>%
  dplyr::select(-cluster)

data.cor <- cor(df_sum_exp_SCT_Dsec_spe)

min(data.cor)

for (ge in Dsec_spe_genes) {
  
  df1 <- df_sum_exp_SCT_Dsec_spe %>%
    dplyr::filter(gene == ge) %>%
    dplyr::select(expression)
  
  cor_re <- cor(df1$expression, )

  
  df_sum_exp_SCT_Dsec_spe_spread <- df_sum_exp_SCT_Dsec_spe %>%
    tidyr::spread(key='gene', value='expression')  
    
  
}


Dsec_spe_genes

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_smooth(size= 0.5, method="glm") +
  geom_point(alpha = 0.8, size=1) +
  geom_label_repel(aes(label=cluster), label.size=0.1, segment.size = 0.1, box.padding = 0.1,
                   label.padding = 0.1, size = 3) +
  aes(x=Dmel, y=Dsec) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title= element_text(color='black', size=11, face='italic')) +
  ggtitle(paste(gen, "Dmel vs Dsec | r =", round(cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec),3)))

