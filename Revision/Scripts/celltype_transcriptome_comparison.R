library(Seurat)
library(tidyverse)
library(cowplot)
library(combinat)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_per_expressed_ref.RData")
load("Processed_Data/shared_genes.RData")
load("Processed_Data/ClusterLabel_order.RData")

TrioBrain.integrated_slim_labeled_final <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_labeled_final.rds")
TrioBrain.integrated_slim_labeled_final_species <- TrioBrain.integrated_slim_labeled_final
TrioBrain.integrated_slim_labeled_final_species$orig.ident <- gsub("_rep[1-6]", "", TrioBrain.integrated_slim_labeled_final_species$orig.ident)
TrioBrain.integrated_slim_labeled_final_species@meta.data$orig.ident <- factor(TrioBrain.integrated_slim_labeled_final_species@meta.data$orig.ident, 
                                                                               levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))
### averge expression per cluster, species ###

sum_exp <- AverageExpression(TrioBrain.integrated_slim_labeled_final_species, features = shared_genes, slot = 'data', add.ident = 'orig.ident')

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

#save(Anno_idents, df_sum_exp_SCT, file="Processed_Data/df_sum_exp_SCT.RData")
load("Processed_Data/df_sum_exp_SCT.RData")

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

length(unique(df_per_expressed_ref_diff$gene))-length(unique(df_per_expressed_ref_diff_filter$gene))
(length(unique(df_per_expressed_ref_diff$gene))-length(unique(df_per_expressed_ref_diff_filter$gene)))/length(unique(df_per_expressed_ref_diff$gene))

df_per_expressed_ref_diff %>%
  dplyr::filter(cluster == "DOP_1", species=="Dsec") %>%
  ggplot(.) +
  geom_point(size=0.7, alpha=0.6) +
  aes(x=MEL, y=OWN) +
  theme_bw()

df_per_expressed_ref_diff_filter %>%
  dplyr::filter(cluster == "DOP_1", species=="Dsec") %>%
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

PRN_matrix <- df_sum_exp_SCT_melrank %>%
  dplyr::filter(cluster =="PRN") %>%
  dplyr::select(gene, Dmel, Dsim, Dsec)

write.csv(PRN_matrix, file="PRN.csv", quote=F, row.names=F)


clu_exp_cell <- c("AST","CTX","ENS","PRN","SUB")
clu_exp_cell <- c("OPN",  "Poxn_3", "γ-KC", "NPF_AstA")

plot_ClusterLabel_exp_top50_melsec <- df_sum_exp_SCT_melrank %>%
  dplyr::filter(cluster %in% clu_exp_cell) %>%
  ggplot(.) +
  geom_point() +
  geom_smooth(method='lm') +
  aes(x=log2(Dmel), y=log2(Dsec)) +
  theme_bw() +
  facet_wrap(~factor(cluster, levels=clu_exp_cell),  nrow=1, scales='free')

plot_ClusterLabel_exp_top50_melsim <- df_sum_exp_SCT_melrank %>%
  dplyr::filter(cluster %in% clu_exp_cell) %>%
  ggplot(.) +
  geom_point() +
  geom_smooth(method='lm') +
  aes(x=log2(Dmel), y=log2(Dsim)) +
  theme_bw() +
  facet_wrap(~factor(cluster, levels=clu_exp_cell),  nrow=1, scales='free')

plot_ClusterLabel_exp_top50_simsec <- df_sum_exp_SCT_melrank %>%
  dplyr::filter(cluster %in% clu_exp_cell) %>%
  ggplot(.) +
  geom_point() +
  geom_smooth(method='lm') +
  aes(x=log2(Dsim), y=log2(Dsec)) +
  theme_bw() +
  facet_wrap(~factor(cluster, levels=clu_exp_cell),  nrow=1, scales='free')

plot_ClusterLabel_exp_top50_secNoni <- df_sum_exp_SCT_melrank %>%
  dplyr::filter(cluster %in% clu_exp_cell) %>%
  ggplot(.) +
  geom_point() +
  geom_smooth(method='lm') +
  aes(x=log2(Dsec), y=log2(DsecNoni)) +
  theme_bw() +
  facet_wrap(~factor(cluster, levels=clu_exp_cell),  nrow=1, scales='free')

plot_grid(plot_ClusterLabel_exp_top50_melsec, plot_ClusterLabel_exp_top50_melsim, plot_ClusterLabel_exp_top50_simsec, 
          ncol=1, align = 'hv')

plot_ClusterLabel_exp_top50_melsec

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
  dplyr::filter(cluster %in% clu_exp_cell) %>%
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

ggsave("Plots/Manuscript/Fig3a_scatter_example.png", width=6, height=6)

save(df_sum_exp_SCT_melrank, file="Processed_Data/df_sum_exp_SCT_melrank.RData")  

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
  
  for (clu in Anno_idents) {

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

non_neuronal_cluster <- c("AST","CTX","ENS","PRN","SUB")

df_per_cor_clean %>%
  dplyr::filter(cluster %in% Anno_idents) %>%
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
        legend.position=c(0.85,0.1),
        legend.title=element_blank(),
        axis.title = element_text(size=10, color='black', face='italic'),
        axis.text = element_text(size=9.5, color='black'),
        legend.text = element_text(size=9.5, color='black'),
        legend.background = element_blank()) +
  scale_fill_manual(values=c("red", "blue")) +
  lims(x=c(0.38,0.8), y=c(0.38,0.8)) +
  labs(x="Dmel vs Dsec (ρ)", y="Dmel vs Dsim (ρ)")

ggsave(file="Plots/Manuscript/Fig3b_scatter_ClusterLabel.png", width=5.5, height=5.5)

df_per_cor_clean %>%
  dplyr::filter(cluster %in% Anno_idents) %>%
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

#save(df_per_cor_clean, df_per_cor_clean_gather, file="Processed_Data/ClusterLabel_cor_permutation.RData")
load("Processed_Data/ClusterLabel_cor_permutation.RData")

df_per_cor_clean_gather %>%
  dplyr::filter(cluster %in% Anno_idents, pair != "secNoni") %>%
  ggplot(.) +
  geom_boxplot() +
  #geom_smooth() +
  aes(x=pair, y=cor_rho, color=pair) +
  theme_bw() +
  facet_wrap(~cluster, nrow=4)

df_per_cor_clean_gather %>%
  dplyr::filter(cluster %in% Anno_idents, pair != "secNoni") %>%
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
  dplyr::filter(cluster %in% Anno_idents, pair != "secNoni") %>%
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

ggsave(file="Plots/Manuscript/Fig3c_boxplot_rho.png", width=22, height=5)


## PCA analysis

# Pivot the data to wide format: one row per (species, cluster), columns are genes
wide_data <- df_sum_exp_SCT %>%
  dplyr::filter(species != 'DsecNoni') %>%
  pivot_wider(names_from = gene, values_from = expression) %>%
  drop_na()  # Drop rows with missing values (optional)

# Separate species and cluster identifiers
species_cluster <- wide_data %>%
  dplyr::select(species, cluster)

# Extract the expression data for PCA
expression_data <- wide_data %>%
  dplyr::select(-species, -cluster)

# Perform PCA
pca_result <- prcomp(expression_data, center = TRUE, scale. = TRUE)

# Add PC scores to the original data
pc_scores <- as.data.frame(pca_result$x)
pc_scores <- cbind(species_cluster, pc_scores)

# Visualize PCA (PC1 vs PC2 as an example)

# Create a new column for labels, with labels for specified clusters only
glial_clusters <- c("AST", "CTX", "ENS", "PRN", "SUB")
pc_scores <- pc_scores %>%
  mutate(label = ifelse(cluster %in% glial_clusters, as.character(cluster), NA))

# Plot with labels for specified clusters
ggplot(pc_scores, aes(x = PC1, y = PC2, color = species, shape = factor(cluster))) +
  geom_point(size = 3) +
  geom_text(aes(label = label), vjust = -1, hjust = 1, size = 3, na.rm = TRUE) +  # Add text labels for selected clusters
  scale_shape_manual(values = cluster_shapes) +
  labs(title = "PCA of Gene Expression Data (Glial Clusters Highlighted)",
       x = "PC1", y = "PC2",
       shape = "Cluster", color = "Species") +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank())

ggplot(pc_scores, aes(x = PC3, y = PC4, color = species, shape = factor(cluster))) +
  geom_point(size = 3) +
  geom_text(aes(label = label), vjust = -1, hjust = 1, size = 3, na.rm = TRUE) +  # Add text labels for selected clusters
  scale_shape_manual(values = cluster_shapes) +
  labs(title = "PCA of Gene Expression Data (Glial Clusters Highlighted)",
       x = "PC3", y = "PC4",
       shape = "Cluster", color = "Species") +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank())

# Summary of PCA
summary(pca_result)

# Optional: Scree plot for variance explained by PCs
explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
scree_data <- data.frame(PC = seq_along(explained_variance), Variance = explained_variance)

ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Variance Explained (%)") +
  theme_minimal()

# group_specific analysis - (1) glia

wide_data_glia <- df_sum_exp_SCT %>%
  dplyr::filter(species != 'DsecNoni', cluster %in% glial_clusters) %>%
  pivot_wider(names_from = gene, values_from = expression) %>%
  drop_na()  # Drop rows with missing values (optional)

# Separate species and cluster identifiers
species_cluster_glia <- wide_data_glia %>%
  dplyr::select(species, cluster)

# Extract the expression data for PCA
expression_data_glia <- wide_data_glia %>%
  dplyr::select(-species, -cluster) %>%
  dplyr::select(where(~ var(.) > 0))

# Perform PCA
pca_result_glia <- prcomp(expression_data_glia, center = TRUE, scale. = TRUE)

# Add PC scores to the original data
pc_scores_glia <- as.data.frame(pca_result_glia$x)
pc_scores_glia <- cbind(species_cluster_glia, pc_scores_glia)

# Visualize PCA (PC1 vs PC2 as an example)
# Create a new column for labels, with labels for specified clusters only

pc_scores_glia <- pc_scores_glia %>%
  mutate(label = ifelse(cluster %in% glial_clusters, as.character(cluster), NA))

# Plot with labels for specified clusters
ggplot(pc_scores_glia, aes(x = PC1, y = PC2, color = species, shape = factor(cluster))) +
  geom_point(size = 3) +
  geom_text(aes(label = label), vjust = -1, hjust = 1, size = 3, na.rm = TRUE) +  # Add text labels for selected clusters
  scale_shape_manual(values = cluster_shapes) +
  labs(title = "PCA of Gene Expression Data (Glial Clusters Highlighted)",
       x = "PC1", y = "PC2",
       shape = "Cluster", color = "Species") +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank())

ggplot(pc_scores_glia, aes(x = PC3, y = PC4, color = species, shape = factor(cluster))) +
  geom_point(size = 3) +
  geom_text(aes(label = label), vjust = -1, hjust = 1, size = 3, na.rm = TRUE) +  # Add text labels for selected clusters
  scale_shape_manual(values = cluster_shapes) +
  labs(title = "PCA of Gene Expression Data (Glial Clusters Highlighted)",
       x = "PC3", y = "PC4",
       shape = "Cluster", color = "Species") +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank())

ggplot(pc_scores_glia, aes(x = PC5, y = PC6, color = species, shape = factor(cluster))) +
  geom_point(size = 3) +
  geom_text(aes(label = label), vjust = -1, hjust = 1, size = 3, na.rm = TRUE) +  # Add text labels for selected clusters
  scale_shape_manual(values = cluster_shapes) +
  labs(title = "PCA of Gene Expression Data (Glial Clusters Highlighted)",
       x = "PC5", y = "PC6",
       shape = "Cluster", color = "Species") +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank())

