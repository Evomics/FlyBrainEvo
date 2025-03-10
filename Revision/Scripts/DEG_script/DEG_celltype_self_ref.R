library(Seurat)
library(tidyverse)
library(MAST)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/../.."))

TrioBrain.integrated_slim_labeled_final <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_labeled_final.rds")
load("Processed_Data/celltype_order.RData")

test <- TrioBrain.integrated_slim_labeled_final

test@meta.data$orig.ident <- gsub("_rep1", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep2", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep3", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep4", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep5", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep6", "", test@meta.data$orig.ident)

test@meta.data$orig.ident <- factor(test@meta.data$orig.ident, 
                                    levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))

df_deg_clusters_ref_self <- NULL

for (cl in Anno_idents) {
  
  cluster_select <- subset(test, ident=cl)
  
  Idents(cluster_select) <- cluster_select$orig.ident
  
  df_deg_melsec <- FindMarkers(cluster_select, ident.1 = "Dmel", ident.2 = "Dsec", verbose = T, test.use = "MAST", min.pct = 0.05, logfc.threshold = 0) %>%
    tibble::rownames_to_column('gene') %>%
    dplyr::mutate(cluster=cl, pair="melsec")
  
  df_deg_melsim <- FindMarkers(cluster_select, ident.1 = "Dmel", ident.2 = "Dsim", verbose = T, test.use = "MAST", min.pct = 0.05, logfc.threshold = 0) %>%
    tibble::rownames_to_column('gene') %>%
    dplyr::mutate(cluster=cl, pair="melsim")
  
  df_deg_simsec <- FindMarkers(cluster_select, ident.1 = "Dsim", ident.2 = "Dsec", verbose = T, test.use = "MAST", min.pct = 0.05, logfc.threshold = 0) %>%
    tibble::rownames_to_column('gene') %>%
    dplyr::mutate(cluster=cl, pair="simsec")
  
  df_deg_secNoni <- FindMarkers(cluster_select, ident.1 = "Dsec", ident.2 = "DsecNoni", verbose = T, test.use = "MAST", min.pct = 0.05, logfc.threshold = 0) %>%
    tibble::rownames_to_column('gene') %>%
    dplyr::mutate(cluster=cl, pair="secNoni")
  
  df_deg_clusters_ref_self <- rbind(df_deg_clusters_ref_self, df_deg_melsec, df_deg_melsim, df_deg_simsec, df_deg_secNoni)
  
}

df_deg_clusters_ref_self_info <- df_deg_clusters_ref_self %>%
  dplyr::group_by(gene, pair) %>%
  dplyr::mutate(n_hit_cluster = sum(p_val_adj < 0.05)) %>%
  dplyr::ungroup()

save(df_deg_clusters_ref_self, df_deg_clusters_ref_self_info, 
     file="Processed_Data/cluster_specific_DEGs_ref_self.RData")  

