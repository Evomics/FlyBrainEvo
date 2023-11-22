library(Seurat)
library(tidyverse)
library(cowplot)
library(combinat)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_per_expressed_ref.RData")
load("Processed_Data/shared_genes.RData")
load("Processed_Data/celltype_order.RData")
TrioBrain.integrated_slim_labeled_species <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled_species.rds")

#df_per_expressed_ref <- df_per_expressed_ref %>%
#  dplyr::mutate(cluster=ifelse(cluster=="Cluster9", "MON", cluster))
#save(df_per_expressed_ref, file="Processed_Data/df_per_expressed_ref.RData")

### averge expression per cluster, species ###

sum_exp <- AverageExpression(TrioBrain.integrated_slim_labeled_species, features = shared_genes, slot = 'data', add.ident = 'orig.ident')

df_sum_exp_SCT <- as.data.frame(sum_exp$SCT)  %>%
  tibble::rownames_to_column('gene') %>%
  tidyr::gather(-gene, key='cluster', value='expression') %>%
  dplyr::mutate(species=ifelse(grepl("Dmel", cluster)==T, "Dmel",
                               ifelse(grepl("Dsim", cluster)==T, "Dsim",
                                      ifelse(grepl("DsecNoni", cluster)==T, "DsecNoni","Dsec"))))

df_sum_exp_SCT$cluster <- gsub("_Dmel", "", df_sum_exp_SCT$cluster)
df_sum_exp_SCT$cluster <- gsub("_Dsim", "", df_sum_exp_SCT$cluster)
df_sum_exp_SCT$cluster <- gsub("_DsecNoni", "", df_sum_exp_SCT$cluster)
df_sum_exp_SCT$cluster <- gsub("_Dsec", "", df_sum_exp_SCT$cluster)

#save(anno_cluster, df_sum_exp_SCT, file="Processed_Data/df_sum_exp_SCT.RData")
#load("Processed_Data/df_sum_exp_SCT.RData")

### clean gene list by filtering genome effects ###

df_per_expressed_ref_diff <- df_per_expressed_ref %>%
  tidyr::spread(key="ref", value='pct_exp') %>%
  dplyr::group_by(cluster, species) %>%
  dplyr::mutate(RefMEL_rank_exp=min_rank(-MEL), RefOWN_rank_exp=min_rank(-OWN)) %>%
  dplyr::mutate(rank_diff = abs(RefMEL_rank_exp-RefOWN_rank_exp), rank_diff_ratio=rank_diff/n())

df_per_expressed_ref_diff_filter <- df_per_expressed_ref_diff %>%
  dplyr::filter(rank_diff_ratio < 0.05) %>%
  dplyr::group_by(gene, cluster) %>%
  dplyr::mutate(n_species = n_distinct(species)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_species ==3)

df_per_expressed_ref_diff %>%
  dplyr::filter(cluster == "DOP", species=="Dsec") %>%
  ggplot(.) +
  geom_point(size=0.7, alpha=0.6) +
  aes(x=MEL, y=OWN) +
  theme_bw()

df_per_expressed_ref_diff_filter %>%
  dplyr::filter(cluster == "DOP", species=="Dsec") %>%
  ggplot(.) +
  geom_point(size=0.7, alpha=0.6) +
  aes(x=MEL, y=OWN) +
  theme_bw()

cor(df_per_expressed_ref_diff_filter$MEL, df_per_expressed_ref_diff_filter$OWN)

### scatter plot for top 100 genes example ###

df_top50_gene_cluster_Dmel <- df_sum_exp_SCT %>%
  dplyr::left_join(., dplyr::select(df_per_expressed_ref_diff_filter, gene, species, cluster, n_species), by=c('gene','species','cluster')) %>%
  na.omit() %>%
  dplyr::filter(species == "Dmel") %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
  dplyr::filter(Dmel_rank<=50) %>%
  dplyr::ungroup() %>%
  dplyr::select(-expression, -species)

df_sum_exp_SCT_melrank <- df_sum_exp_SCT %>%
  dplyr::left_join(., df_top50_gene_cluster_Dmel, by=c('gene', 'cluster')) %>%
  na.omit() %>%
  tidyr::spread(key='species', value='expression')

clu_exp_cell <- c("AST","CTX","ENS","PRN","SUB")
clu_exp_cell <- c("OPN",  "PRN", "Clock", "IPC")

plot_celltype_exp_top50_melsec <- df_sum_exp_SCT_melrank %>%
  dplyr::filter(cluster %in% clu_exp_cell) %>%
  ggplot(.) +
  geom_point() +
  geom_smooth(method='lm') +
  aes(x=log2(Dmel), y=log2(Dsec)) +
  theme_bw() +
  facet_wrap(~factor(cluster, levels=clu_exp_cell),  nrow=1, scales='free')

plot_celltype_exp_top50_melsim <- df_sum_exp_SCT_melrank %>%
  dplyr::filter(cluster %in% clu_exp_cell) %>%
  ggplot(.) +
  geom_point() +
  geom_smooth(method='lm') +
  aes(x=log2(Dmel), y=log2(Dsim)) +
  theme_bw() +
  facet_wrap(~factor(cluster, levels=clu_exp_cell),  nrow=1, scales='free')

plot_celltype_exp_top50_simsec <- df_sum_exp_SCT_melrank %>%
  dplyr::filter(cluster %in% clu_exp_cell) %>%
  ggplot(.) +
  geom_point() +
  geom_smooth(method='lm') +
  aes(x=log2(Dsim), y=log2(Dsec)) +
  theme_bw() +
  facet_wrap(~factor(cluster, levels=clu_exp_cell),  nrow=1, scales='free')

plot_celltype_exp_top50_secNoni <- df_sum_exp_SCT_melrank %>%
  dplyr::filter(cluster %in% clu_exp_cell) %>%
  ggplot(.) +
  geom_point() +
  geom_smooth(method='lm') +
  aes(x=log2(Dsec), y=log2(DsecNoni)) +
  theme_bw() +
  facet_wrap(~factor(cluster, levels=clu_exp_cell),  nrow=1, scales='free')

plot_grid(plot_celltype_exp_top50_melsec, plot_celltype_exp_top50_melsim, plot_celltype_exp_top50_simsec, 
          ncol=1, align = 'hv')

plot_celltype_exp_top50_melsec

df_sum_exp_SCT_melrank %>%
  dplyr::filter(cluster %in% clu_exp_cell) %>%
  ggplot(.) +
  geom_point(aes(x=log2(Dmel), y=log2(Dsec)), color='red', alpha =0.7) +
  geom_point(aes(x=log2(Dmel), y=log2(Dsim)), color='blue', alpha =0.7) +
  geom_smooth(method='lm', aes(x=log2(Dmel), y=log2(Dsec)), color='red') +
  #geom_smooth(method='lm', aes(x=log2(Dmel), y=log2(Dsim)), color='blue') +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~cluster, scales='free', nrow=1)

df_sum_exp_SCT_melrank %>%
  dplyr::filter(cluster %in% c("IPC", "PRN", "OPN", "AstA")) %>%
  ggplot(.) +
  geom_point(aes(x=log2(Dmel), y=log2(Dsec)), fill='orange', alpha =0.6, shape=21, stroke=0.3) +
  geom_point(aes(x=log2(Dmel), y=log2(Dsim)), fill='navy', alpha =0.6, shape=21, stroke=0.3) +
  geom_smooth(method='lm', aes(x=log2(Dmel), y=log2(Dsec)), color='orange', fill='orange', alpha =0.2, size=0.5) +
  geom_smooth(method='lm', aes(x=log2(Dmel), y=log2(Dsim)), color='navy', fill='navy', alpha =0.2, size=0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=10, color='black'),
        axis.text = element_text(size=9.5, color='black'),
        strip.text = element_text(size=10, color='black')) +
  labs(x="log2(expression) - Dmel", y= "log2(expression) - Dsim or Dsec") +
  facet_wrap(~cluster, nrow=2, scales='free') 

ggsave("Plots/Manuscript/Fig3a_scatter_example.pdf", width=6, height=6)
  

### permutation based correlation ###

### random permutation ##

df_per_cor_clean <- NULL

gene_sel=50
group=30 ## size of group

set.seed(123)

for (per in 1:400) {
  
  df_top_gene_cluster_Dmel <- df_sum_exp_SCT %>%
    dplyr::left_join(., dplyr::select(df_per_expressed_ref_diff_filter, gene, species, cluster, n_species), by=c('gene','species','cluster')) %>%
    na.omit() %>%
    dplyr::filter(species == "Dmel") %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
    dplyr::filter(Dmel_rank<=gene_sel) %>%
    dplyr::sample_n(group) %>%
    dplyr::ungroup()
  
  for (clu in CellType_order$CellType) {

    df_select_gene <- df_top_gene_cluster_Dmel %>%
      dplyr::filter(cluster==clu)   

    df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
      dplyr::filter(cluster==clu, gene %in% df_select_gene$gene) %>%
      tidyr::spread(key='species', value='expression')
    
    cc_cor_melsec <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method="spearman")
    cc_cor_melsim <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsim, method="spearman")
    cc_cor_simsec <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$Dsim, method="spearman")
    cc_cor_secNoni <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$DsecNoni, method="spearman")
    
    df_per_cor_clean_temp <- data.frame(cluster=clu, 
                                        gene_pop=gene_sel,
                                        trial=per,
                                        group_size=group,
                                        melsec=cc_cor_melsec,
                                        melsim=cc_cor_melsim, 
                                        simsec=cc_cor_simsec, 
                                        secNoni=cc_cor_secNoni)
    
    df_per_cor_clean <- rbind(df_per_cor_clean, df_per_cor_clean_temp)
    
  }
  
}

non_neuronal_cluster <- c("AST","CTX","ENS","PRN","SUB", "Fat", "Blood")

df_per_cor_clean %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(melsec=mean(melsec), melsim=mean(melsim)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(group = ifelse(cluster %in% non_neuronal_cluster, "Non-neuronal", "Neuronal")) %>%
  ggplot(.) +
  geom_abline(intercept = 0, slope = 1, linetype=2, alpha = 0.5) +
  #geom_smooth(aes(x=melsec, y=melsim), method='lm') +
  geom_point(aes(x=melsec, y=melsim, fill=group), shape=21) +
  geom_text_repel(aes(x=melsec, y=melsim, label=cluster), segment.size = 0.1, size = 3) +
  #geom_text_repel(aes(x=melsec, y=melsim, 
  #                    label=factor(cluster, levels =c("PRN", "SUB", "Fat", "Ty", "ort", "AST", "AstA", "fru(Ms)", "Crz", "OCTY", "MON"))) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position=c(0.15,0.9),
        legend.title=element_blank(),
        axis.title = element_text(size=10, color='black', face='italic'),
        axis.text = element_text(size=9.5, color='black'),
        legend.text = element_text(size=9.5, color='black'),
        legend.background = element_blank()) +
  scale_fill_manual(values=c("red", "blue")) +
  lims(x=c(0.35,0.78), y=c(0.35,0.78)) +
  labs(x="Dmel vs Dsec (ρ)", y="Dmel vs Dsim (ρ)")

ggsave(file="Plots/Manuscript/Fig3b_scatter_celltype.pdf", width=5.5, height=5.5)

df_per_cor_clean %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(melsec=mean(melsec), simsec=mean(simsec)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(group = ifelse(cluster %in% non_neuronal_cluster, "Non-neuronal", "Neuronal")) %>%
  ggplot(.) +
  geom_abline(intercept = 0, slope = 1, linetype=2, alpha = 0.5) +
  #geom_smooth(aes(x=melsec, y=melsim), method='lm') +
  geom_point(aes(x=melsec, y=simsec, fill=group), shape=21) +
  geom_text_repel(aes(x=melsec, y=simsec, label=cluster)) +
  #geom_text_repel(aes(x=melsec, y=melsim, 
  #                    label=factor(cluster, levels =c("PRN", "SUB", "Fat", "Ty", "ort", "AST", "AstA", "fru(Ms)", "Crz", "OCTY", "MON"))) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.position=c(0.85,0.1)) +
  scale_fill_manual(values=c("red", "blue")) +
  labs(x="Dmel vs Dsec", y="Dsim vs Dsec")


df_per_cor_clean_gather <- df_per_cor_clean %>%
  tidyr::gather(-cluster, -trial, -group_size, -gene_pop, key='pair', value='cor_rho')

df_per_cor_clean_gather_summary <- df_per_cor_clean_gather %>%
  dplyr::group_by(cluster, pair) %>%
  dplyr::summarise(mean_rho = mean(cor_rho))

#save(df_per_cor_clean, df_per_cor_clean_gather, file="Processed_Data/celltype_cor_permutation.RData")
load("Processed_Data/celltype_cor_permutation.RData")

df_per_cor_clean_gather %>%
  dplyr::filter(cluster %in% anno_cluster, pair != "secNoni") %>%
  ggplot(.) +
  geom_boxplot() +
  #geom_smooth() +
  aes(x=pair, y=cor_rho, color=pair) +
  theme_bw() +
  facet_wrap(~cluster, nrow=4)

df_per_cor_clean_gather %>%
  dplyr::filter(cluster %in% anno_cluster, pair != "secNoni") %>%
  ggplot(.) +
  geom_boxplot() +
  #geom_smooth() +
  aes(x=reorder(cluster,cor_rho), y=cor_rho, color=pair) +
  theme_bw() +
  facet_wrap(~pair, nrow=1)

df_order_melsec <- df_per_cor_clean_gather %>%
  dplyr::filter(pair == "melsec") %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(value=mean(cor_rho)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(value)

df_per_cor_clean_gather %>%
  dplyr::filter(cluster %in% anno_cluster, pair != "secNoni") %>%
  ggplot(.) +
  geom_boxplot(outlier.alpha = 0.7, outlier.size = 0.5, size=0.3, alpha = 0.7) +
  #geom_smooth() +
  aes(x=factor(cluster,levels=df_order_melsec$cluster), y=cor_rho, fill=pair) +
  scale_fill_manual(values=c("orange", "navy",  "brown")) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position=c(0.92,0.15),
        legend.title=element_blank(),
        axis.title.y = element_text(size=10, color='black'),
        axis.text.x = element_text(size=9.5, color='black', angle=45, hjust =1, vjust=1),
        axis.text.y = element_text(size=9.5, color='black'),
        legend.text = element_text(size=9.5, color='black')) +
  labs(y="Spearman's rho (ρ)")

ggsave(file="Plots/Manuscript/Fig3c_boxplot_rho.pdf", width=13, height=5)


df_per_cor_clean_gather %>%
  dplyr::filter(cluster %in% anno_cluster, pair == "melsec") %>%
  ggplot(.) +
  geom_boxplot(fill='blue') +
  #geom_smooth() +
  aes(x=reorder(cluster,cor_rho), y=cor_rho) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  facet_grid(~pair)

df_per_cor_clean_gather %>%
  dplyr::filter(cluster %in% anno_cluster, pair == "melsim") %>%
  ggplot(.) +
  geom_boxplot( fill='red') +
  #geom_smooth() +
  aes(x=reorder(cluster,cor_rho), y=cor_rho) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  facet_grid(~pair)

df_per_cor_clean_gather %>%
  dplyr::filter(cluster %in% anno_cluster, pair == "simsec") %>%
  ggplot(.) +
  geom_boxplot( fill='brown') +
  #geom_smooth() +
  aes(x=reorder(cluster,cor_rho), y=cor_rho) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  facet_grid(~pair)

df_per_cor_clean_gather %>%
  dplyr::filter(cluster %in% anno_cluster, pair == "secNoni") %>%
  ggplot(.) +
  geom_boxplot( fill='purple') +
  #geom_smooth() +
  aes(x=reorder(cluster,cor_rho), y=cor_rho) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  facet_grid(~pair)


df_summary_comp <- CellType_order %>%
  dplyr::rename(cluster=CellType)

df_per_cor_clean_gather_freq <- df_per_cor_clean_gather %>%
  dplyr::left_join(., df_summary_comp, by='cluster')

df_per_cor_clean_gather_freq %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  dplyr::group_by(cluster, pair, mean_freq) %>%
  dplyr::summarise(mean_rho=mean(cor_rho)) %>%
  dplyr::ungroup() %>%
  ggplot(.) +
  geom_smooth(aes(x=log2(mean_freq), y=mean_rho), alpha =0.5, method='lm') +
  geom_point(aes(x=log2(mean_freq), y=mean_rho, fill=cluster), shape=21) +
  geom_text_repel(aes(x=log2(mean_freq), y=mean_rho,label=cluster)) +
  theme_bw() +
  facet_grid(~pair)

## statistical test - with different numbers of permutations ##

df_per_cor_posthoc_celltype <- NULL

for (n_per in c(25,50,100,200,400)) {
  
  df_per_cor_clean_gather_posthoc <- df_per_cor_clean_gather %>%
    dplyr::filter(pair != "secNoni", trial <= n_per) %>%
    dplyr::select(cluster, pair, cor_rho) %>%
    dplyr::group_by(cluster) %>%
    do(multitst = tidy(TukeyHSD(aov(cor_rho ~ pair, data = .))))
  
  for (cl in unique(df_per_cor_clean_gather_posthoc$cluster)) {
    
    df1 <- df_per_cor_clean_gather_posthoc %>%
      dplyr::filter(cluster==cl)
    
    df2 <- df1$multitst[[1]] %>%
      dplyr::mutate(cluster=cl, Perm=n_per)
    
    df_per_cor_posthoc_celltype <- rbind(df_per_cor_posthoc_celltype, df2)
    
  }
  
}

### need to do p-adjustement for multiple testing ### -> p-value is inflated 

df_per_cor_posthoc_celltype_padj <- df_per_cor_posthoc_celltype %>%
  dplyr::rowwise() %>%
  dplyr::mutate(adj2.p.value = p.adjust(adj.p.value, method='bonferroni', n=Perm)) %>%
  dplyr::mutate(adj3.p.value = p.adjust(adj2.p.value, method='bonferroni', n=length(unique(df_per_cor_posthoc_celltype$cluster)))) %>%
  dplyr::ungroup()

df_per_cor_posthoc_celltype_padj %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  ggplot(.) +
  geom_point() +
  #geom_smooth() +
  aes(x=factor(cluster, levels=CellType_order_anno$CellType), y=-log10(adj3.p.value), color=as.factor(Perm)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  facet_grid(~contrast)


df_per_cor_posthoc_celltype_fdr <- df_per_cor_posthoc_celltype %>%
  dplyr::group_by(contrast, Perm) %>%
  dplyr::mutate(adj2.p.value = p.adjust(adj.p.value, method="fdr")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sig = ifelse(adj2.p.value < 0.05, T, F))

df_per_cor_posthoc_celltype_bf <- df_per_cor_posthoc_celltype %>%
  dplyr::group_by(contrast, Perm) %>%
  dplyr::mutate(adj2.p.value = p.adjust(adj.p.value, method="bonferroni")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sig = ifelse(adj2.p.value < 0.05, T, F))

df_per_cor_clean_gather %>%
  dplyr::filter(cluster %in% c("PRN"), pair != "secNoni") %>%
  ggplot(.) +
  geom_boxplot(alpha=0.8, outlier.alpha = 0) +
  geom_jitter(width=0.2) +
  #geom_smooth() +
  aes(x=pair, y=cor_rho, color=pair) +
  theme_bw() +
  facet_wrap(~cluster, nrow=4)



## scatter plot ##
ws=50
sw=10
clu='AST'

df_top_gene_cluster_Dmel <- df_sum_exp_SCT %>%
  dplyr::filter(species == "Dmel", !gene %in% background_deg$gene) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
  dplyr::ungroup()  %>%
  dplyr::filter(Dmel_rank >=sw & Dmel_rank <= sw+ws-1)

df_select_gene <- df_top_gene_cluster_Dmel %>%
  dplyr::filter(cluster==clu)

df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
  dplyr::filter(cluster==clu, gene %in% df_select_gene$gene) %>%
  tidyr::spread(key='species', value='expression')

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_point() +
  geom_smooth() +
  aes(x=log2(Dmel), y=log2(Dsec)) +
  geom_text_repel(aes(label=gene)) + 
  theme_bw() 

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_point() +
  geom_smooth() +
  aes(x=log2(Dmel), y=log2(Dsim)) +
  geom_text_repel(aes(label=gene)) + 
  theme_bw() 

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_point() +
  geom_smooth() +
  aes(x=log2(Dsim), y=log2(Dsec)) +
  geom_text_repel(aes(label=gene)) + 
  theme_bw() 

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_point() +
  geom_smooth() +
  aes(x=log2(DsecNoni), y=log2(Dsec)) +
  geom_text_repel(aes(label=gene)) + 
  theme_bw() 




test2 <- df_deg_cluster_fc_pval_ref %>%
  dplyr::filter(cluster == "PRN", pair=='melsec') %>%
  dplyr::mutate(RefMEL_rank_fc=min_rank(log2FC_MEL), RefOWN_rank_fc=min_rank(log2FC_OWN)) %>%
  dplyr::mutate(rank_diff = abs(RefMEL_rank_fc-RefOWN_rank_fc), rank_diff_ratio=rank_diff/n()) 

nrow(test2)
nrow(dplyr::filter(test2, rank_diff_ratio < 0.1))




df_top_gene_cluster_Dmel <- df_sum_exp_SCT %>%
  dplyr::filter(species == "Dmel") %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
  dplyr::ungroup()  %>%
  dplyr::filter(Dmel_rank <=100)

df_top100_cor <- NULL

for (clu in CellType_order$CellType) {
  
  df_select_gene <- df_top_gene_cluster_Dmel %>%
    dplyr::filter(cluster==clu)
  
  df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
    dplyr::filter(cluster==clu, gene %in% df_select_gene$gene) %>%
    tidyr::spread(key='species', value='expression')
  
  cc_cor_melsec <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec)
  cc_cor_melsim <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsim)
  cc_cor_simsec <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$Dsim)
  cc_cor_secNoni <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$DsecNoni)
  
  df_top100_cor_temp <- data.frame(cluster=clu, 
                               melsec=cc_cor_melsec,
                               melsim=cc_cor_melsim, 
                               simsec=cc_cor_simsec, 
                               secNoni=cc_cor_secNoni)
  
  df_top100_cor <- rbind(df_top100_cor, df_top100_cor_temp)
  
}

df_top100_cor_gather <- df_top100_cor %>%
  tidyr::gather(-cluster, key='pair', value='cor')

df_top100_cor_gather %>%
  ggplot(.) +
  geom_boxplot() +
  aes(x=pair, y=cor) +
  theme_bw()

df_top100_cor %>%
  ggplot(.) +
  geom_point() +
  aes(x=simsec, y=secNoni, color=cluster) +
  geom_text_repel(aes(label=cluster)) +
  theme_bw()

df_top100_cor %>%
  ggplot(.) +
  geom_point() +
  aes(x=simsec, y=melsec, color=cluster) +
  geom_text_repel(aes(label=cluster)) +
  theme_bw()

df_top100_cor %>%
  ggplot(.) +
  geom_point() +
  aes(x=simsec, y=melsim, color=cluster) +
  geom_text_repel(aes(label=cluster)) +
  theme_bw()

df_top100_cor %>%
  ggplot(.) +
  geom_point() +
  aes(x=melsec, y=melsim, color=cluster) +
  geom_text_repel(aes(label=cluster)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw()

df_TF_cluster_top10_Dmel <- df_sum_exp_SCT %>%
  dplyr::filter(species == "Dmel", gene %in% TF_genes) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene, cluster, Dmel_rank) %>%
  dplyr::filter(Dmel_rank <=10)

### top 50 ###

df_top_gene_cluster_Dmel <- df_sum_exp_SCT %>%
  dplyr::filter(species == "Dmel") %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
  dplyr::ungroup()  %>%
  dplyr::filter(Dmel_rank <=50)

df_top50_cor <- NULL

for (clu in CellType_order$CellType) {
  
  df_select_gene <- df_top_gene_cluster_Dmel %>%
    dplyr::filter(cluster==clu)
  
  df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
    dplyr::filter(cluster==clu, gene %in% df_select_gene$gene) %>%
    tidyr::spread(key='species', value='expression')
  
  cc_cor_melsec <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method="spearman")
  cc_cor_melsim <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsim, method="spearman")
  cc_cor_simsec <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$Dsim, method="spearman")
  cc_cor_secNoni <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$DsecNoni, method="spearman")
  
  df_top50_cor_temp <- data.frame(cluster=clu, 
                                   melsec=cc_cor_melsec,
                                   melsim=cc_cor_melsim, 
                                   simsec=cc_cor_simsec, 
                                   secNoni=cc_cor_secNoni)
  
  df_top50_cor <- rbind(df_top50_cor, df_top50_cor_temp)
  
}

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_point() +
  aes(x=Dmel, y=Dsec) +
  geom_text_repel(aes(label=gene)) + 
  theme_bw() +
  lims(x=c(2,8), y=c(0,10))

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_point() +
  aes(x=Dmel, y=Dsim) +
  geom_text_repel(aes(label=gene)) + 
  theme_bw()

df_top50_cor_gather <- df_top50_cor %>%
  tidyr::gather(-cluster, key='pair', value='cor_rho')

df_top50_cor_gather %>%
  ggplot(.) +
  geom_boxplot() +
  aes(x=pair, y=cor_rho) +
  theme_bw()

df_top50_cor %>%
  ggplot(.) +
  geom_point() +
  aes(x=simsec, y=secNoni, color=cluster) +
  geom_text_repel(aes(label=cluster)) +
  theme_bw()

df_top50_cor %>%
  ggplot(.) +
  geom_point(aes(x=simsec, y=melsec, color=cluster)) +
  geom_smooth(aes(x=simsec, y=melsec)) +
  geom_text_repel(aes(x=simsec, y=melsec, label=cluster)) +
  theme_bw()

df_top50_cor %>%
  ggplot(.) +
  geom_point(aes(x=simsec, y=melsim, color=cluster)) +
  geom_smooth(aes(x=simsec, y=melsim)) +
  geom_text_repel(aes(x=simsec, y=melsim, label=cluster)) +
  theme_bw()

df_top50_cor %>%
  ggplot(.) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(aes(x=melsec, y=melsim)) +
  geom_point(aes(x=melsec, y=melsim, color=cluster)) +
  geom_text_repel(aes(x=melsec, y=melsim, label=cluster)) +
  theme_bw()

### top 100 ###

df_top_gene_cluster_Dmel <- df_sum_exp_SCT %>%
  dplyr::filter(species == "Dmel") %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
  dplyr::ungroup()  %>%
  dplyr::filter(Dmel_rank <=100)

df_top100_cor <- NULL

for (clu in CellType_order$CellType) {
  
  df_select_gene <- df_top_gene_cluster_Dmel %>%
    dplyr::filter(cluster==clu)
  
  df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
    dplyr::filter(cluster==clu, gene %in% df_select_gene$gene) %>%
    tidyr::spread(key='species', value='expression')
  
  cc_cor_melsec <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method="spearman")
  cc_cor_melsim <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsim, method="spearman")
  cc_cor_simsec <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$Dsim, method="spearman")
  cc_cor_secNoni <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$DsecNoni, method="spearman")
  
  df_top100_cor_temp <- data.frame(cluster=clu, 
                                  melsec=cc_cor_melsec,
                                  melsim=cc_cor_melsim, 
                                  simsec=cc_cor_simsec, 
                                  secNoni=cc_cor_secNoni)
  
  df_top100_cor <- rbind(df_top100_cor, df_top100_cor_temp)
  
}

df_top100_cor_gather <- df_top100_cor %>%
  tidyr::gather(-cluster, key='pair', value='cor_rho')

df_top100_cor_gather %>%
  ggplot(.) +
  geom_boxplot() +
  aes(x=pair, y=cor_rho) +
  theme_bw()

df_top100_cor %>%
  ggplot(.) +
  geom_point() +
  aes(x=simsec, y=secNoni, color=cluster) +
  geom_text_repel(aes(label=cluster)) +
  theme_bw()

df_top100_cor %>%
  ggplot(.) +
  geom_point(aes(x=simsec, y=melsec, color=cluster)) +
  geom_smooth(aes(x=simsec, y=melsec)) +
  geom_text_repel(aes(x=simsec, y=melsec, label=cluster)) +
  theme_bw()

df_top100_cor %>%
  ggplot(.) +
  geom_point(aes(x=simsec, y=melsim, color=cluster)) +
  geom_smooth(aes(x=simsec, y=melsim)) +
  geom_text_repel(aes(x=simsec, y=melsim, label=cluster)) +
  theme_bw()

df_top100_cor %>%
  ggplot(.) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(aes(x=melsec, y=melsim)) +
  geom_point(aes(x=melsec, y=melsim, color=cluster)) +
  geom_text_repel(aes(x=melsec, y=melsim, label=cluster)) +
  theme_bw()


### top 200 ###

df_top_gene_cluster_Dmel <- df_sum_exp_SCT %>%
  dplyr::filter(species == "Dmel") %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
  dplyr::ungroup()  %>%
  dplyr::filter(Dmel_rank <=200)

df_top200_cor <- NULL

for (clu in CellType_order$CellType) {
  
  df_select_gene <- df_top_gene_cluster_Dmel %>%
    dplyr::filter(cluster==clu)
  
  df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
    dplyr::filter(cluster==clu, gene %in% df_select_gene$gene) %>%
    tidyr::spread(key='species', value='expression')
  
  cc_cor_melsec <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method="spearman")
  cc_cor_melsim <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsim, method="spearman")
  cc_cor_simsec <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$Dsim, method="spearman")
  cc_cor_secNoni <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$DsecNoni, method="spearman")
  
  df_top200_cor_temp <- data.frame(cluster=clu, 
                                   melsec=cc_cor_melsec,
                                   melsim=cc_cor_melsim, 
                                   simsec=cc_cor_simsec, 
                                   secNoni=cc_cor_secNoni)
  
  df_top200_cor <- rbind(df_top200_cor, df_top200_cor_temp)
  
}

df_top200_cor_gather <- df_top200_cor %>%
  tidyr::gather(-cluster, key='pair', value='cor_rho')

df_top200_cor_gather %>%
  ggplot(.) +
  geom_boxplot() +
  aes(x=pair, y=cor_rho) +
  theme_bw()

df_top200_cor %>%
  ggplot(.) +
  geom_point() +
  aes(x=simsec, y=secNoni, color=cluster) +
  geom_text_repel(aes(label=cluster)) +
  theme_bw()

df_top200_cor %>%
  ggplot(.) +
  geom_point(aes(x=simsec, y=melsec, color=cluster)) +
  geom_smooth(aes(x=simsec, y=melsec)) +
  geom_text_repel(aes(x=simsec, y=melsec, label=cluster)) +
  theme_bw()

df_top200_cor %>%
  ggplot(.) +
  geom_point(aes(x=simsec, y=melsim, color=cluster)) +
  geom_smooth(aes(x=simsec, y=melsim)) +
  geom_text_repel(aes(x=simsec, y=melsim, label=cluster)) +
  theme_bw()

df_top200_cor %>%
  ggplot(.) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(aes(x=melsec, y=melsim)) +
  geom_point(aes(x=melsec, y=melsim, color=cluster)) +
  geom_text_repel(aes(x=melsec, y=melsim, label=cluster)) +
  theme_bw()

### top 30 ###

df_top_gene_cluster_Dmel <- df_sum_exp_SCT %>%
  dplyr::filter(species == "Dmel") %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
  dplyr::ungroup()  %>%
  dplyr::filter(Dmel_rank <=30)

df_top30_cor <- NULL

for (clu in CellType_order$CellType) {
  
  df_select_gene <- df_top_gene_cluster_Dmel %>%
    dplyr::filter(cluster==clu)
  
  df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
    dplyr::filter(cluster==clu, gene %in% df_select_gene$gene) %>%
    tidyr::spread(key='species', value='expression')
  
  cc_cor_melsec <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method="spearman")
  cc_cor_melsim <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsim, method="spearman")
  cc_cor_simsec <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$Dsim, method="spearman")
  cc_cor_secNoni <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$DsecNoni, method="spearman")
  
  df_top30_cor_temp <- data.frame(cluster=clu, 
                                   melsec=cc_cor_melsec,
                                   melsim=cc_cor_melsim, 
                                   simsec=cc_cor_simsec, 
                                   secNoni=cc_cor_secNoni)
  
  df_top30_cor <- rbind(df_top30_cor, df_top30_cor_temp)
  
}


df_top30_cor_gather <- df_top30_cor %>%
  tidyr::gather(-cluster, key='pair', value='cor_rho')

df_top30_cor_gather %>%
  ggplot(.) +
  geom_boxplot() +
  aes(x=pair, y=cor_rho) +
  theme_bw()

df_top30_cor %>%
  ggplot(.) +
  geom_point() +
  aes(x=simsec, y=secNoni, color=cluster) +
  geom_text_repel(aes(label=cluster)) +
  theme_bw()

df_top30_cor %>%
  ggplot(.) +
  geom_point(aes(x=simsec, y=melsec, color=cluster)) +
  geom_smooth(aes(x=simsec, y=melsec)) +
  geom_text_repel(aes(x=simsec, y=melsec, label=cluster)) +
  theme_bw()

df_top30_cor %>%
  ggplot(.) +
  geom_point(aes(x=simsec, y=melsim, color=cluster)) +
  geom_smooth(aes(x=simsec, y=melsim)) +
  geom_text_repel(aes(x=simsec, y=melsim, label=cluster)) +
  theme_bw()

df_top30_cor %>%
  ggplot(.) +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(aes(x=melsec, y=melsim)) +
  geom_point(aes(x=melsec, y=melsim, color=cluster)) +
  geom_text_repel(aes(x=melsec, y=melsim, label=cluster)) +
  theme_bw()


### n of genes and cor ###

df_top_gene_cluster_Dmel <- df_sum_exp_SCT %>%
  dplyr::filter(species == "Dmel") %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
  dplyr::ungroup()  %>%
  dplyr::filter(Dmel_rank <=200)

df_topn_cor <- NULL

for (ng in seq(20,400, by=10)) {
  
  df_top_gene_cluster_Dmel <- df_sum_exp_SCT %>%
    dplyr::filter(species == "Dmel") %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
    dplyr::ungroup()  %>%
    dplyr::filter(Dmel_rank <= ng)
  
  for (clu in CellType_order$CellType) {
    
    df_select_gene <- df_top_gene_cluster_Dmel %>%
      dplyr::filter(cluster==clu)
    
    df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
      dplyr::filter(cluster==clu, gene %in% df_select_gene$gene) %>%
      tidyr::spread(key='species', value='expression')
    
    cc_cor_melsec <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method="spearman")
    cc_cor_melsim <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsim, method="spearman")
    cc_cor_simsec <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$Dsim, method="spearman")
    cc_cor_secNoni <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$DsecNoni, method="spearman")
    
    df_topn_cor_temp <- data.frame(cluster=clu, 
                                     n_gene=ng,
                                     melsec=cc_cor_melsec,
                                     melsim=cc_cor_melsim, 
                                     simsec=cc_cor_simsec, 
                                     secNoni=cc_cor_secNoni)
    
    df_topn_cor <- rbind(df_topn_cor, df_topn_cor_temp)
    
  }

}

df_topn_cor_gather <- df_topn_cor %>%
  tidyr::gather(-cluster, -n_gene, key='pair', value='cor_rho')

df_topn_cor_gather %>%
  dplyr::filter(cluster %in% c("OPN", "IPC", "ENS", "AST", "PRN", "CTX", "SUB", "Da6", "Da5/7", "Tbh"), pair == "melsec", n_gene <=200, n_gene >=30) %>%
  ggplot(.) +
  geom_line() +
  aes(x=n_gene, y=cor_rho, color=cluster) +
  theme_bw()

df_topn_cor_gather %>%
  #dplyr::filter(cluster %in% c("OPN", "IPC", "ENS", "AST", "PRN", "CTX", "SUB", "Da6", "Da5/7", "Tbh"), pair == "melsec", n_gene <=200, n_gene >=30) %>%
  ggplot(.) +
  #geom_line() +
  geom_smooth() +
  aes(x=n_gene, y=cor_rho) +
  theme_bw() +
  facet_wrap(~pair, scales='free', ncol=1)





## filtering baseline degs ##

test_degs <- rbind(test_degs_melsec, test_degs_melsim, test_degs_simsec)

background_deg <- test_degs %>%
  dplyr::filter(abs(avg_log2FC) > log2(1.5))

### varying number of genes for input ##

df_topn_cor_clean <- NULL

for (ng in seq(10,200, by=5)) {
  
  df_top_gene_cluster_Dmel <- df_sum_exp_SCT %>%
    dplyr::filter(species == "Dmel", !gene %in% background_deg$gene) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
    dplyr::ungroup()  %>%
    dplyr::filter(Dmel_rank <= ng)
  
  for (clu in CellType_order$CellType) {
    
    df_select_gene <- df_top_gene_cluster_Dmel %>%
      dplyr::filter(cluster==clu)
    
    df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
      dplyr::filter(cluster==clu, gene %in% df_select_gene$gene) %>%
      tidyr::spread(key='species', value='expression')
    
    cc_cor_melsec <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method="spearman")
    cc_cor_melsim <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsim, method="spearman")
    cc_cor_simsec <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$Dsim, method="spearman")
    cc_cor_secNoni <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$DsecNoni, method="spearman")
    
   df_topn_cor_clean_temp <- data.frame(cluster=clu, 
                                   n_gene=ng,
                                   melsec=cc_cor_melsec,
                                   melsim=cc_cor_melsim, 
                                   simsec=cc_cor_simsec, 
                                   secNoni=cc_cor_secNoni)
    
   df_topn_cor_clean <- rbind(df_topn_cor_clean, df_topn_cor_clean_temp)
    
  }
  
}

df_topn_cor_clean_gather <-df_topn_cor_clean %>%
  tidyr::gather(-cluster, -n_gene, key='pair', value='cor_rho')

df_topn_cor_clean_gather %>%
  dplyr::filter(pair != "secNoni") %>%
  #dplyr::filter(cluster %in% c("OPN", "IPC", "ENS", "AST", "PRN", "CTX", "SUB", "Da6", "Da5/7", "Tbh"), pair == "melsec", n_gene <=200, n_gene >=30) %>%
  ggplot(.) +
  #geom_line() +
  geom_smooth() +
  aes(x=n_gene, y=cor_rho, color=pair) +
  theme_bw()

df_topn_cor_clean_gather %>%
  dplyr::filter(pair != "secNoni") %>%
  dplyr::filter(cluster %in% c("OPN", "IPC", "ENS", "AST", "PRN", "CTX", "SUB")) %>%
  ggplot(.) +
  geom_line() +
  #geom_smooth() +
  aes(x=n_gene, y=cor_rho, color=pair) +
  theme_bw() +
  facet_grid(~cluster)

df_topn_cor_clean_gather %>%
  dplyr::filter(cluster %in% anno_cluster, pair != "secNoni") %>%
  dplyr::filter(n_gene <=200) %>%
  ggplot(.) +
  geom_boxplot() +
  #geom_smooth() +
  aes(x=pair, y=cor_rho, color=pair) +
  theme_bw() +
  facet_wrap(~cluster, nrow=4)

df_topn_cor_clean %>%
  dplyr::filter(n_gene ==50) %>%
  ggplot(.) +
  geom_abline(intercept = 0, slope = 1, linetype=2, alpha = 0.5) +
  geom_smooth(aes(x=melsec, y=melsim)) +
  geom_point(aes(x=melsec, y=melsim, color=cluster)) +
  geom_text_repel(aes(x=melsec, y=melsim, label=cluster)) +
  theme_bw()

df_summary_comp <- CellType_order %>%
  dplyr::rename(cluster=CellType)

df_topn_cor_clean_gather_freq <- df_topn_cor_clean_gather %>%
  dplyr::left_join(., df_summary_comp, by='cluster')

df_topn_cor_clean_gather_freq %>%
  dplyr::filter(n_gene==50, cluster %in% anno_cluster) %>%
  ggplot(.) +
  geom_smooth(aes(x=log2(mean_freq), y=cor_rho), alpha =0.5) +
  geom_point(aes(x=log2(mean_freq), y=cor_rho, fill=cluster), shape=21) +
  geom_text_repel(aes(x=log2(mean_freq), y=cor_rho,label=cluster)) +
  theme_bw() +
  facet_grid(~pair)

df_topn_cor_clean_gather_freq %>%
  dplyr::filter(n_gene==50, cluster %in% anno_cluster) %>%
  ggplot(.) +
  geom_smooth(aes(x=log2(mean_freq), y=cor_rho, group=cluster), alpha =0.5) +
  geom_point(aes(x=log2(mean_freq), y=cor_rho, color=cluster, shape=pair)) +
  geom_text_repel(aes(x=log2(mean_freq), y=cor_rho,label=cluster)) +
  theme_bw()

## statistical test ##

df_topn_cor_clean_gather_posthoc <- df_topn_cor_clean_gather %>%
  dplyr::filter(pair != "secNoni") %>%
  dplyr::select(cluster, pair, cor_rho) %>%
  dplyr::group_by(cluster) %>%
  do(multitst = tidy(TukeyHSD(aov(cor_rho ~ pair, data = .))))

df_cor_posthoc_celltype <- NULL

for (cl in unique(df_topn_cor_clean_gather_posthoc$cluster)) {
  
  df1 <- df_topn_cor_clean_gather_posthoc %>%
    dplyr::filter(cluster==cl)
  
  df2 <- df1$multitst[[1]] %>%
    dplyr::mutate(cluster=cl)

  df_cor_posthoc_celltype <- rbind(df_cor_posthoc_celltype, df2)
  
}

df_cor_posthoc_celltype_fdr <- df_cor_posthoc_celltype %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(adj2.p.value = p.adjust(adj.p.value, method="fdr")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sig = ifelse(adj2.p.value < 0.05, T, F))

## scatter plots ##

ng=50
clu='Tk'

df_top_gene_cluster_Dmel <- df_sum_exp_SCT %>%
  dplyr::filter(species == "Dmel", !gene %in% background_deg$gene) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
  dplyr::ungroup()  %>%
  dplyr::filter(Dmel_rank <= ng)

df_select_gene <- df_top_gene_cluster_Dmel %>%
  dplyr::filter(cluster==clu)

df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
  dplyr::filter(cluster==clu, gene %in% df_select_gene$gene) %>%
  tidyr::spread(key='species', value='expression')

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_point() +
  geom_smooth() +
  aes(x=log2(Dmel), y=log2(Dsec)) +
  geom_text_repel(aes(label=gene)) + 
  theme_bw() 

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_point() +
  geom_smooth() +
  aes(x=log2(Dmel), y=log2(Dsim)) +
  geom_text_repel(aes(label=gene)) + 
  theme_bw() 

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_point() +
  geom_smooth() +
  aes(x=log2(Dsim), y=log2(Dsec)) +
  geom_text_repel(aes(label=gene)) + 
  theme_bw() 

df_sum_exp_SCT_select %>%
  #dplyr::filter(!gene %in% c("dnc","Dop2R","CG12239","CG17684","pros","Tk","IA-2", "cpo", "Wdr62", "Rbfox1")) %>%
  ggplot(.) +
  geom_point() +
  aes(x=log2(Dmel), y=log2(Dsec)) +
  geom_text_repel(aes(label=gene)) + 
  theme_bw() 

df_sum_exp_SCT_select %>%
  dplyr::filter(!gene %in% c("dnc","Dop2R","CG12239","CG17684","pros","Tk","IA-2", "cpo", "Wdr62", "Rbfox1")) %>%
  ggplot(.) +
  geom_point() +
  aes(x=log2(Dsim), y=log2(Dsec)) +
  geom_text_repel(aes(label=gene)) + 
  theme_bw() 


### sliding windows ##

df_sw_cor_clean <- NULL

ws=30 ## size of window

for (sw in seq(1,30, by=1)) {
  
  df_top_gene_cluster_Dmel <- df_sum_exp_SCT %>%
    dplyr::filter(species == "Dmel", !gene %in% background_deg$gene) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
    dplyr::ungroup()  %>%
    dplyr::filter(Dmel_rank >=sw & Dmel_rank <= sw+ws-1)
  
  for (clu in CellType_order$CellType) {
    
    df_select_gene <- df_top_gene_cluster_Dmel %>%
      dplyr::filter(cluster==clu)
    
    df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
      dplyr::filter(cluster==clu, gene %in% df_select_gene$gene) %>%
      tidyr::spread(key='species', value='expression')
    
    cc_cor_melsec <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method="spearman")
    cc_cor_melsim <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsim, method="spearman")
    cc_cor_simsec <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$Dsim, method="spearman")
    cc_cor_secNoni <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$DsecNoni, method="spearman")
    
    df_sw_cor_clean_temp <- data.frame(cluster=clu, 
                                         window_start=sw,
                                         window_size=ws,
                                         melsec=cc_cor_melsec,
                                         melsim=cc_cor_melsim, 
                                         simsec=cc_cor_simsec, 
                                         secNoni=cc_cor_secNoni)
    
    df_sw_cor_clean <- rbind(df_sw_cor_clean, df_sw_cor_clean_temp)
    
  }
  
}

df_sw_cor_clean_gather <-df_sw_cor_clean %>%
  tidyr::gather(-cluster, -window_size, -window_start, key='pair', value='cor_rho')

df_sw_cor_clean_gather %>%
  dplyr::filter(pair != "secNoni") %>%
  ggplot(.) +
  geom_point(alpha=0.5, size= 0.7) +
  geom_smooth() +
  aes(x=window_start, y=cor_rho, color=pair) +
  theme_bw()

df_sw_cor_clean_gather %>%
  #dplyr::filter(pair != "secNoni") %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  ggplot(.) +
  geom_line() +
  geom_point() +
  #geom_smooth() +
  aes(x=window_start, y=cor_rho, color=pair) +
  theme_bw() +
  facet_wrap(~cluster, nrow=4)

df_sw_cor_clean_gather %>%
  dplyr::filter(pair != "secNoni") %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  ggplot(.) +
  #geom_line() +
  geom_point(aes(x=window_start, y=cor_rho, color=pair) ) +
  geom_smooth(aes(x=window_start, y=cor_rho) ) +
  theme_bw() +
  facet_wrap(~cluster, nrow=4)

df_sw_cor_clean_gather %>%
  dplyr::filter(cluster %in% anno_cluster, pair != "secNoni") %>%
  ggplot(.) +
  geom_boxplot() +
  #geom_smooth() +
  aes(x=pair, y=cor_rho, color=pair) +
  theme_bw() +
  facet_wrap(~cluster, nrow=4)



## statistical test ##

df_sw_cor_clean_gather_posthoc <- df_sw_cor_clean_gather %>%
  dplyr::filter(pair != "secNoni") %>%
  dplyr::select(cluster, pair, cor_rho) %>%
  dplyr::group_by(cluster) %>%
  do(multitst = tidy(TukeyHSD(aov(cor_rho ~ pair, data = .))))

df_sw_cor_posthoc_celltype <- NULL

for (cl in unique(df_sw_cor_clean_gather_posthoc$cluster)) {
  
  df1 <- df_sw_cor_clean_gather_posthoc %>%
    dplyr::filter(cluster==cl)
  
  df2 <- df1$multitst[[1]] %>%
    dplyr::mutate(cluster=cl)
  
  df_sw_cor_posthoc_celltype <- rbind(df_sw_cor_posthoc_celltype, df2)
  
}

df_sw_cor_posthoc_celltype_fdr <- df_sw_cor_posthoc_celltype %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(adj2.p.value = p.adjust(adj.p.value, method="fdr")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sig = ifelse(adj2.p.value < 0.05, T, F))

df_sw_cor_clean_gather %>%
  dplyr::filter(cluster %in% c("PRN"), pair != "secNoni") %>%
  ggplot(.) +
  geom_boxplot(alpha=0.8, outlier.alpha = 0) +
  geom_jitter(width=0.2) +
  #geom_smooth() +
  aes(x=pair, y=cor_rho, color=pair) +
  theme_bw() +
  facet_wrap(~cluster, nrow=4)



## scatter plot ##
ws=30
sw=10
clu='OPN'

df_top_gene_cluster_Dmel <- df_sum_exp_SCT %>%
  dplyr::filter(species == "Dmel", !gene %in% background_deg$gene) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(Dmel_rank=min_rank(-expression)) %>%
  dplyr::ungroup()  %>%
  dplyr::filter(Dmel_rank >=sw & Dmel_rank <= sw+ws-1)

df_select_gene <- df_top_gene_cluster_Dmel %>%
  dplyr::filter(cluster==clu)

df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
  dplyr::filter(cluster==clu, gene %in% df_select_gene$gene) %>%
  tidyr::spread(key='species', value='expression')

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_point() +
  geom_smooth() +
  aes(x=log2(Dmel), y=log2(Dsec)) +
  geom_text_repel(aes(label=gene)) + 
  theme_bw() 

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_point() +
  geom_smooth() +
  aes(x=log2(Dmel), y=log2(Dsim)) +
  geom_text_repel(aes(label=gene)) + 
  theme_bw() 

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_point() +
  geom_smooth() +
  aes(x=log2(Dsim), y=log2(Dsec)) +
  geom_text_repel(aes(label=gene)) + 
  theme_bw() 


