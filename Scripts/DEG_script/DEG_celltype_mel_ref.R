library(Seurat)
library(tidyverse)
library(MAST)

load("/work/FAC/FBM/CIG/rbenton/neuroflies/Daehan/FlyBrainEvo/R/ref_DEG/bin/Trio_Dmelref_DF_labeled.RData")

test <- Trio_Dmelref_DF_labeled

test@meta.data$orig.ident <- gsub("_rep1", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep2", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep3", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep4", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep5", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep6", "", test@meta.data$orig.ident)

test@meta.data$orig.ident <- factor(test@meta.data$orig.ident, 
                                    levels=c("Dmel", "Dsim_to_Dmel", "Dsec_to_Dmel"))

cluster_labels <- as.character(unique(test@active.ident))

df_deg_clusters_ref_mel <- NULL

for (cl in cluster_labels) {
  
  cluster_select <- subset(test, ident=cl)
  
  Idents(cluster_select) <- cluster_select$orig.ident
  
  df_deg_melsec <- FindMarkers(cluster_select, ident.1 = "Dmel", ident.2 = "Dsec_to_Dmel", verbose = T, test.use = "MAST", min.pct = 0.05, logfc.threshold = 0) %>%
    tibble::rownames_to_column('gene') %>%
    dplyr::mutate(cluster=cl, pair="melsec")
  
  df_deg_melsim <- FindMarkers(cluster_select, ident.1 = "Dmel", ident.2 = "Dsim_to_Dmel", verbose = T, test.use = "MAST", min.pct = 0.05, logfc.threshold = 0) %>%
    tibble::rownames_to_column('gene') %>%
    dplyr::mutate(cluster=cl, pair="melsim")
  
  df_deg_simsec <- FindMarkers(cluster_select, ident.1 = "Dsim_to_Dmel", ident.2 = "Dsec_to_Dmel", verbose = T, test.use = "MAST", min.pct = 0.05, logfc.threshold = 0) %>%
    tibble::rownames_to_column('gene') %>%
    dplyr::mutate(cluster=cl, pair="simsec")

  df_deg_clusters_ref_mel <- rbind(df_deg_clusters_ref_mel, df_deg_melsec, df_deg_melsim, df_deg_simsec)
  
}

df_deg_clusters_ref_mel_info <- df_deg_clusters_ref_mel %>%
  dplyr::group_by(gene, pair) %>%
  dplyr::mutate(n_hit_cluster = sum(p_val_adj < 0.05)) %>%
  dplyr::ungroup()

save(df_deg_clusters_ref_mel, df_deg_clusters_ref_mel_info, 
     file="/work/FAC/FBM/CIG/rbenton/neuroflies/Daehan/FlyBrainEvo/R/ref_DEG/out/cluster_specific_DEGs_refmel.RData")  

