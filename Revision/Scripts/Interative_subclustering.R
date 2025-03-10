library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(sctransform)
library(cowplot)
library(parallel)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/shared_genes.RData")

#TrioBrain.integrated_slim <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim.rds")
TrioBrain.integrated_slim<- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim.rds")
TrioBrain.integrated.markers <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated.markers.rds")

df_TrioBrain.integrated_slim_metadata <- TrioBrain.integrated_slim@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))

total_cell_number <- nrow(TrioBrain.integrated_slim@meta.data)

metadata_summary <- TrioBrain.integrated_slim@meta.data %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-n) %>%
  dplyr::mutate(percent = n/total_cell_number*100)

ISub <- function(Seurat_object, cluster, method='rpca', resolution = 0.1, npc=50, genelist=shared_genes, kweight=100) {
  
  Seurat <- subset(Seurat_object, idents = cluster)
  
  DefaultAssay(Seurat) <- "RNA"
  
  Dmel <- subset(Seurat, subset = orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"))
  Dsim <- subset(Seurat, subset = orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"))
  Dsec <- subset(Seurat, subset = orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"))
  DsecNoni <- subset(Seurat, subset = orig.ident %in% c("DsecNoni_rep1", "DsecNoni_rep2", "DsecNoni_rep3", "DsecNoni_rep4", "DsecNoni_rep5", "DsecNoni_rep6"))  
  
  Seurat_list <- list(Dmel, Dsim, Dsec, DsecNoni) 
  
  Seurat_list <- lapply(X = Seurat_list, FUN = SCTransform, method = "glmGamPoi", verbose = FALSE)
  
  Int_features <- SelectIntegrationFeatures(object.list = Seurat_list, nfeatures = 3000, verbose = FALSE)
  
  Seurat_list <- PrepSCTIntegration(object.list = Seurat_list, assay = "SCT", 
                                    anchor.features = Int_features,
                                    verbose = F)
  
  Seurat_list <- lapply(X = Seurat_list, FUN = RunPCA, features = Int_features, verbose = FALSE)
  
  anchors<- FindIntegrationAnchors(object.list = Seurat_list, normalization.method = "SCT",
                                   anchor.features = Int_features, verbose = F, reduction = method, reference=1)
  
  Seurat_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = F,
                                     features.to.integrate = genelist,
                                     k.weight = kweight)
  
  Seurat_integrated <- RunPCA(Seurat_integrated, npcs = npc, verbose = F)
  Seurat_integrated <- FindNeighbors(Seurat_integrated, reduction = "pca", dims = 1:npc, verbose = F)
  Seurat_integrated <- FindClusters(Seurat_integrated, resolution = resolution, verbose = F)
  
  DefaultAssay(Seurat_integrated) <- 'RNA'
  Seurat_integrated <- DietSeurat(Seurat_integrated, assays = c("RNA"))
  
  return(Seurat_integrated)
  
}

DefaultAssay(TrioBrain.integrated_slim) <- 'RNA'
TrioBrain.integrated_slimer <- DietSeurat(TrioBrain.integrated_slim, assays = c("RNA"))
DefaultAssay(TrioBrain.integrated_slim) <- 'SCT'

cluster_range_max <- max(as.numeric(as.vector(dplyr::filter(metadata_summary, percent > 0.5)$seurat_clusters)))

TrioBrain.integrated_slim_ISub_DF <- TrioBrain.integrated_slim

remain1 <- NULL

DF_sd_threshold <- 0.5

for (cl in c(cluster_range_max:7, 5:0)) {
  
  name <- paste0('DF_cluster',cl)
  assign(name, ISub(TrioBrain.integrated_slimer, cluster = cl, resolution = 0.1))

  meta_doublet <- get(name)@meta.data %>%
    mutate(total_cell_count=n(), 
           total_mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
           total_sd_nFeature_RNA = sd(nFeature_RNA, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::mutate(mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
                  sd_nFeature_RNA = sd(nFeature_RNA, na.rm = TRUE),
                  fraction_within_cluter = n()/total_cell_count) %>%
    dplyr::mutate(sd_lower = mean_nFeature_RNA-sd_nFeature_RNA, 
                  sd_upper = mean_nFeature_RNA+sd_nFeature_RNA) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(DF_lower = total_mean_nFeature_RNA-DF_sd_threshold*total_sd_nFeature_RNA, 
                  DF_upper = total_mean_nFeature_RNA+DF_sd_threshold*total_sd_nFeature_RNA) %>%
    dplyr::distinct(seurat_clusters, total_mean_nFeature_RNA, total_sd_nFeature_RNA, 
                    mean_nFeature_RNA, sd_nFeature_RNA, sd_lower, sd_upper,
                    DF_lower, DF_upper, fraction_within_cluter)
  
  meta_doublet_filtered <- meta_doublet %>%
    dplyr::filter(mean_nFeature_RNA>=total_mean_nFeature_RNA+DF_sd_threshold*total_sd_nFeature_RNA)
  
  doublet_clusters <- as.character(meta_doublet_filtered$seurat_clusters)
 
  meta_doublet %>%
      ggplot(.) +
      geom_ribbon(aes(ymin = DF_lower, ymax = DF_upper), fill = "blue", alpha = 0.2) +
      geom_line(aes(y=total_mean_nFeature_RNA), linetype=2, alpha=0.5, color='red') +
      geom_point() +
      geom_errorbar(aes(ymin = sd_lower, ymax = sd_upper), width = 0.1, alpha=0.7) +
      geom_line(alpha=0.5, linewidth=0.3) +
      aes(x=fraction_within_cluter, y=mean_nFeature_RNA) +
      theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x="Fraction in cluster", y="Mean nFeature RNA")
  
  ggsave(glue::glue("Plots/Isub_DF/cluster{cl}_DF_range.png"), width=5, height=5)

  name_DF_cellid <- WhichCells(get(name), idents=doublet_clusters)
  Idents(TrioBrain.integrated_slim_ISub_DF, cells = name_DF_cellid) <- "Doublets"
  
  meta_sum <- get(name)@meta.data %>%
    dplyr::filter(!seurat_clusters %in% doublet_clusters) %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(percent = n/total_cell_number*100)
  
  non_doublet_clusters <- as.character(meta_sum$seurat_clusters)
      
  ## Number of clusters > 0.1%
  
  n_cluster_pass <- meta_sum %>%
    dplyr::filter(percent>0.1) %>%
    nrow()
  
  if (n_cluster_pass >=2) {
    
    # subA - size가 0.5% 이상이면 해당 cluster는 subclustering 진행
    
    df1 <- meta_sum %>%
      dplyr::filter(percent >= 0.5)
    
    for (subA in df1$seurat_clusters) {
      
      name_subA <- paste0('cluster',cl,'_',subA)
      assign(name_subA, ISub(get(name), cluster = subA, resolution = 0.1, kweight=100))
      
      remain1 <- c(remain1, name_subA)
      
    }
    
    # subB - size가 0.1% 이상 0.5% 미만이면 해당 cluster에 대해선 subclustering 종료 및 cell_id 할당
    
    df2 <- meta_sum %>%
      dplyr::filter(percent >= 0.1 & percent <0.5)
    
    for (subB in df2$seurat_clusters) {
      
      name_subB <- paste0('cluster',cl,'_',subB)
      name_subB_cellid <- WhichCells(get(name), idents=subB)
      Idents(TrioBrain.integrated_slim_ISub_DF, cells = name_subB_cellid) <- name_subB
      
    }
    
    # subC - size가 0.1% 미만 cluster는 모 클러스터에 소속, 모 클러스터 cell_id 할당. 추후에 cluster size에 따라 모두 re-order & re-name    
    
    df3 <- meta_sum %>%
      dplyr::filter(percent < 0.1)
    
    for (subC in df3$seurat_clusters) {
      
      name_subC_cellid <- WhichCells(get(name), idents=subC)
      Idents(TrioBrain.integrated_slim_ISub_DF, cells = name_subC_cellid) <- cl
      
    }
    
  }
  
  else if (n_cluster_pass < 2) {
    
    # n_cluster_pass가 하나뿐인 경우, 즉 major subcluster가 안 분리되는 경우에는 전체 subcluster를 하나의 모 cluster-id로 확정 및 cell_id 할당
    
    name_cellid <- WhichCells(subset(get(name), idents= non_doublet_clusters))
    Idents(TrioBrain.integrated_slim_ISub_DF, cells = name_cellid) <- cl
    
  }
  
  rm(list = name)
  
}

for (cl in 6) {
  
  name <- paste0('DF_cluster',cl)
  assign(name, ISub(TrioBrain.integrated_slimer, cluster = cl, resolution = 0.1))
  
  meta_doublet <- get(name)@meta.data %>%
    mutate(total_cell_count=n(), 
           total_mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
           total_sd_nFeature_RNA = sd(nFeature_RNA, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::mutate(mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
                  sd_nFeature_RNA = sd(nFeature_RNA, na.rm = TRUE),
                  fraction_within_cluter = n()/total_cell_count) %>%
    dplyr::mutate(sd_lower = mean_nFeature_RNA-sd_nFeature_RNA, 
                  sd_upper = mean_nFeature_RNA+sd_nFeature_RNA) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(DF_lower = total_mean_nFeature_RNA-DF_sd_threshold*total_sd_nFeature_RNA, 
                  DF_upper = total_mean_nFeature_RNA+DF_sd_threshold*total_sd_nFeature_RNA) %>%
    dplyr::distinct(seurat_clusters, total_mean_nFeature_RNA, total_sd_nFeature_RNA, 
                    mean_nFeature_RNA, sd_nFeature_RNA, sd_lower, sd_upper,
                    DF_lower, DF_upper, fraction_within_cluter)
  
  meta_doublet_filtered <- meta_doublet %>%
    dplyr::filter(mean_nFeature_RNA>=total_mean_nFeature_RNA+DF_sd_threshold*total_sd_nFeature_RNA)
  
  doublet_clusters <- as.character(meta_doublet_filtered$seurat_clusters)
  
  meta_doublet %>%
    ggplot(.) +
    geom_ribbon(aes(ymin = DF_lower, ymax = DF_upper), fill = "blue", alpha = 0.2) +
    geom_line(aes(y=total_mean_nFeature_RNA), linetype=2, alpha=0.5, color='red') +
    geom_point() +
    geom_errorbar(aes(ymin = sd_lower, ymax = sd_upper), width = 0.1, alpha=0.7) +
    geom_line(alpha=0.5, linewidth=0.3) +
    aes(x=fraction_within_cluter, y=mean_nFeature_RNA) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x="Fraction in cluster", y="Mean nFeature RNA")
  
  ggsave(glue::glue("Plots/Isub_DF/cluster{cl}_DF_range.png"), width=5, height=5)
  
  name_DF_cellid <- WhichCells(get(name), idents=doublet_clusters)
  Idents(TrioBrain.integrated_slim_ISub_DF, cells = name_DF_cellid) <- "Doublets"
  
  meta_sum <- get(name)@meta.data %>%
    dplyr::filter(!seurat_clusters %in% doublet_clusters) %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(percent = n/total_cell_number*100)
  
  non_doublet_clusters <- as.character(meta_sum$seurat_clusters)
  
  ## Number of clusters > 0.1%
  
  n_cluster_pass <- meta_sum %>%
    dplyr::filter(percent>0.1) %>%
    nrow()
  
  if (n_cluster_pass >=2) {
    
    # subA - size가 0.5% 이상이면 해당 cluster는 subclustering 진행
    
    df1 <- meta_sum %>%
      dplyr::filter(percent >= 0.5)
    
    for (subA in df1$seurat_clusters) {
      
      name_subA <- paste0('cluster',cl,'_',subA)
      assign(name_subA, ISub(get(name), cluster = subA, resolution = 0.1, kweight=60))
      
      remain1 <- c(remain1, name_subA)
      
    }
    
    # subB - size가 0.1% 이상 0.5% 미만이면 해당 cluster에 대해선 subclustering 종료 및 cell_id 할당
    
    df2 <- meta_sum %>%
      dplyr::filter(percent >= 0.1 & percent <0.5)
    
    for (subB in df2$seurat_clusters) {
      
      name_subB <- paste0('cluster',cl,'_',subB)
      name_subB_cellid <- WhichCells(get(name), idents=subB)
      Idents(TrioBrain.integrated_slim_ISub_DF, cells = name_subB_cellid) <- name_subB
      
    }
    
    # subC - size가 0.1% 미만 cluster는 모 클러스터에 소속, 모 클러스터 cell_id 할당. 추후에 cluster size에 따라 모두 re-order & re-name    
    
    df3 <- meta_sum %>%
      dplyr::filter(percent < 0.1)
    
    for (subC in df3$seurat_clusters) {
      
      name_subC_cellid <- WhichCells(get(name), idents=subC)
      Idents(TrioBrain.integrated_slim_ISub_DF, cells = name_subC_cellid) <- cl
      
    }
    
  }
  
  else if (n_cluster_pass < 2) {
    
    # n_cluster_pass가 하나뿐인 경우, 즉 major subcluster가 안 분리되는 경우에는 전체 subcluster를 하나의 모 cluster-id로 확정 및 cell_id 할당
    
    name_cellid <- WhichCells(subset(get(name), idents= non_doublet_clusters))
    Idents(TrioBrain.integrated_slim_ISub_DF, cells = name_cellid) <- cl
    
  }
  
  rm(list = name)
  
} ## kweight =60

### interative subclustering

for (r in 1:9) {
  
  remainders <- paste0('remain', r)
  n_remainders <- length(get(remainders))
  remainders_new <- paste0('remain', (r+1))
  assign(remainders_new, NULL)
  
  for (j in 1:n_remainders) {
    
    remain_name <- get(remainders)[j]
    remain_seurat <- get(get(remainders)[j])
    
    meta_doublet <- remain_seurat@meta.data %>%
      mutate(total_cell_count=n(), 
             total_mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
             total_sd_nFeature_RNA = sd(nFeature_RNA, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(seurat_clusters) %>%
      dplyr::mutate(mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
                    sd_nFeature_RNA = sd(nFeature_RNA, na.rm = TRUE),
                    fraction_within_cluter = n()/total_cell_count) %>%
      dplyr::mutate(sd_lower = mean_nFeature_RNA-sd_nFeature_RNA, 
                    sd_upper = mean_nFeature_RNA+sd_nFeature_RNA) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(DF_lower = total_mean_nFeature_RNA-DF_sd_threshold*total_sd_nFeature_RNA, 
                    DF_upper = total_mean_nFeature_RNA+DF_sd_threshold*total_sd_nFeature_RNA) %>%
      dplyr::distinct(seurat_clusters, total_mean_nFeature_RNA, total_sd_nFeature_RNA, 
                      mean_nFeature_RNA, sd_nFeature_RNA, sd_lower, sd_upper,
                      DF_lower, DF_upper, fraction_within_cluter)
    
    meta_doublet_filtered <- meta_doublet %>%
      dplyr::filter(mean_nFeature_RNA>=total_mean_nFeature_RNA+DF_sd_threshold*total_sd_nFeature_RNA)
    
    doublet_clusters <- as.character(meta_doublet_filtered$seurat_clusters)
    
    meta_doublet %>%
      ggplot(.) +
      geom_ribbon(aes(ymin = DF_lower, ymax = DF_upper), fill = "blue", alpha = 0.2) +
      geom_line(aes(y=total_mean_nFeature_RNA), linetype=2, alpha=0.5, color='red') +
      geom_point() +
      geom_errorbar(aes(ymin = sd_lower, ymax = sd_upper), width = 0.1, alpha=0.7) +
      geom_line(alpha=0.5, linewidth=0.3) +
      aes(x=fraction_within_cluter, y=mean_nFeature_RNA) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      labs(x="Fraction in cluster", y="Mean nFeature RNA")
    
    ggsave(glue::glue("Plots/Isub_DF/{remain_name}_DF_range.png"), width=5, height=5)
    
    name_DF_cellid <- WhichCells(remain_seurat, idents=doublet_clusters)
    Idents(TrioBrain.integrated_slim_ISub_DF, cells = name_DF_cellid) <- "Doublets"
    
    meta_sum <- remain_seurat@meta.data %>%
      dplyr::filter(!seurat_clusters %in% doublet_clusters) %>%
      dplyr::group_by(seurat_clusters) %>%
      dplyr::summarise(n=n()) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(-n) %>%
      dplyr::mutate(percent = n/total_cell_number*100)
    
    non_doublet_clusters <- as.character(meta_sum$seurat_clusters)

    ## Number of clusters > 0.1%
    
    n_cluster_pass <- meta_sum %>%
      dplyr::filter(percent>0.1) %>%
      nrow()
    
    if (n_cluster_pass >=2) {
      
      # subA - size가 0.5% 이상이면 해당 cluster는 subclustering 진행
      
      df1 <- meta_sum %>%
        dplyr::filter(percent >= 0.5)
      
      for (subA in df1$seurat_clusters) {
        
        remainders_subA <- paste0(remain_name,'_',subA)
        assign(remainders_subA, ISub(remain_seurat, cluster = subA, resolution = 0.1, kweight=80))
        
        assign(remainders_new, c(get(remainders_new), remainders_subA))
        
      }
      
      # subB - size가 0.1% 이상 0.5% 미만이면 해당 cluster에 대해선 subclustering 종료 및 cell_id 할당
      
      df2 <- meta_sum %>%
        dplyr::filter(percent >= 0.1 & percent <0.5)
      
      for (subB in df2$seurat_clusters) {
        
        remainders_subB <- paste0(remain_name,'_',subB)
        remainders_subB_cellid <- WhichCells(remain_seurat, idents=subB)
        Idents(TrioBrain.integrated_slim_ISub_DF, cells = remainders_subB_cellid) <- remainders_subB
        
      }
      
      # subC - size가 0.1% 미만 cluster는 모 클러스터에 소속, 모 클러스터 cell_id 할당. 추후에 cluster size에 따라 모두 re-order & re-name    
      
      df3 <- meta_sum %>%
        dplyr::filter(percent < 0.1)
      
      for (subC in df3$seurat_clusters) {
        
        remainders_subC_cellid <- WhichCells(remain_seurat, idents=subC)
        Idents(TrioBrain.integrated_slim_ISub_DF, cells = remainders_subC_cellid) <- remain_name
        
      }
      
    }
    
    else if (n_cluster_pass < 2) {
      
      # n_cluster_pass가 하나뿐인 경우, 즉 major subcluster가 안 분리되는 경우에는 전체 subcluster를 하나의 모 cluster-id로 확정 및 cell_id 할당
      
      remainders_cellid <- WhichCells(subset(remain_seurat, idents= non_doublet_clusters))
      Idents(TrioBrain.integrated_slim_ISub_DF, cells = remainders_cellid) <- remain_name
      
    }
    
    
  }
  
  
} ## iterative subclustering for clusters > 0.5%

saveRDS(TrioBrain.integrated_slim_ISub_DF, file = "Processed_Data/TrioBrain.integrated_slim_ISub_DF.rds")

unique(Idents(TrioBrain.integrated_slim_ISub_DF))

TrioBrain.integrated_slim_ISub_DF[["CellType"]] <- Idents(object = TrioBrain.integrated_slim_ISub_DF)

TrioBrain.integrated_slim_ISub_DF_metadata <- TrioBrain.integrated_slim_ISub_DF@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))


TrioBrain.integrated_slim_ISub_DF_metadata_summary_rep <- TrioBrain.integrated_slim_ISub_DF_metadata %>%
  dplyr::group_by(species, CellType, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

TrioBrain.integrated_slim_ISub_DF_metadata_summary <- TrioBrain.integrated_slim_ISub_DF_metadata_summary_rep %>%
  dplyr::group_by(species, CellType) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

TrioBrain.integrated_slim_ISub_DF_metadata_summary %>%
  ggplot(.) +
  geom_col(aes(y=CellType, x=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.9)) + 
  geom_errorbar(aes(y=CellType, x=percent_combined,
                    xmin=ifelse(percent_combined-sem<0,0,percent_combined-sem), xmax=percent_combined+sem, group=species), width=0.9, position=position_dodge(width=0.9)) +
  geom_point(data=TrioBrain.integrated_slim_ISub_DF_metadata_summary_rep, size = 1, alpha =0.8,  
             aes(y=CellType, x=percent, group=species), fill='red', shape=21, position=position_dodge(width=0.9)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=CellType, y=n/total*100)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11, color='black'),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=10, color='black'),
        axis.text.y=element_text(size=10, color='black'),
        panel.grid = element_blank()) +
  labs(y="Cluster", x="Percent of cell type (%)", fill="") +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3"))

ggsave("Plots/Isub_DF_clusters_frequencies_all.pdf", width=13, height = 33)
