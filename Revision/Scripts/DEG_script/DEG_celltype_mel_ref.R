library(Seurat)
library(tidyverse)
library(MAST)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/../.."))

load(file="Processed_Data/Trio_Dmelref_DF_labeled.RData")
load("Processed_Data/celltype_order.RData")

test2 <- Trio_Dmelref_DF_labeled

test2@meta.data$orig.ident <- gsub("_rep1", "", test2@meta.data$orig.ident)
test2@meta.data$orig.ident <- gsub("_rep2", "", test2@meta.data$orig.ident)
test2@meta.data$orig.ident <- gsub("_rep3", "", test2@meta.data$orig.ident)
test2@meta.data$orig.ident <- gsub("_rep4", "", test2@meta.data$orig.ident)
test2@meta.data$orig.ident <- gsub("_rep5", "", test2@meta.data$orig.ident)
test2@meta.data$orig.ident <- gsub("_rep6", "", test2@meta.data$orig.ident)

test2@meta.data$orig.ident <- factor(test2@meta.data$orig.ident, 
                                    levels=c("Dmel", "Dsim_to_DmelRef", "Dsec_to_DmelRef"))

df_deg_clusters_ref_mel <- NULL

for (cl in Anno_idents) {
  
  cluster_select <- subset(test2, ident=cl)
  
  Idents(cluster_select) <- cluster_select$orig.ident
  
  df_deg_melsec <- FindMarkers(cluster_select, ident.1 = "Dmel", ident.2 = "Dsec_to_DmelRef", verbose = T, test2.use = "MAST", min.pct = 0.05, logfc.threshold = 0) %>%
    tibble::rownames_to_column('gene') %>%
    dplyr::mutate(cluster=cl, pair="melsec")
  
  df_deg_melsim <- FindMarkers(cluster_select, ident.1 = "Dmel", ident.2 = "Dsim_to_DmelRef", verbose = T, test2.use = "MAST", min.pct = 0.05, logfc.threshold = 0) %>%
    tibble::rownames_to_column('gene') %>%
    dplyr::mutate(cluster=cl, pair="melsim")
  
  df_deg_simsec <- FindMarkers(cluster_select, ident.1 = "Dsim_to_DmelRef", ident.2 = "Dsec_to_DmelRef", verbose = T, test2.use = "MAST", min.pct = 0.05, logfc.threshold = 0) %>%
    tibble::rownames_to_column('gene') %>%
    dplyr::mutate(cluster=cl, pair="simsec")

  df_deg_clusters_ref_mel <- rbind(df_deg_clusters_ref_mel, df_deg_melsec, df_deg_melsim, df_deg_simsec)
  
}

df_deg_clusters_ref_mel_info <- df_deg_clusters_ref_mel %>%
  dplyr::group_by(gene, pair) %>%
  dplyr::mutate(n_hit_cluster = sum(p_val_adj < 0.05)) %>%
  dplyr::ungroup()

save(df_deg_clusters_ref_mel, df_deg_clusters_ref_mel_info, 
     file="Processed_Data/cluster_specific_DEGs_refmel.RData")  

