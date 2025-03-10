library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(sctransform)
library(cowplot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/shared_genes.RData")

Dsim_to_DmelRef_to_DmelRef_decon_combined_DF <- readRDS(file = "Processed_Data/Dsim_to_DmelRef_to_DmelRef_combined_full_intron_decon_DF.rds")
Dsec_to_DmelRef_decon_combined_DF <- readRDS(file = "Processed_Data/Dsec_to_DmelRef_combined_full_intron_decon_DF.rds")
DsecNoni_to_DmelRef_decon_combined_DF <- readRDS(file = "Processed_Data/DsecNoni_to_DmelRef_combined_full_intron_decon_DF.rds")

Dsim_to_DmelRef_Brain_ortho_only_DF <- subset(Dsim_to_DmelRef_decon_combined_DF, features = shared_genes)
Dsec_to_DmelRef_Brain_ortho_only_DF <- subset(Dsec_to_DmelRef_decon_combined_DF, features = shared_genes)
DsecNoni_to_DmelRef_Brain_ortho_only_DF <- subset(DsecNoni_to_DmelRef_decon_combined_DF, features = shared_genes)

DefaultAssay(Dsim_to_DmelRef_Brain_ortho_only_DF) <- "SCT"
DefaultAssay(Dsec_to_DmelRef_Brain_ortho_only_DF) <- "SCT"
DefaultAssay(DsecNoni_to_DmelRef_Brain_ortho_only_DF) <- "SCT"

secsim_DmelRef_Brain_DF.list <- list(Dsim_to_DmelRef_Brain_ortho_only_DF, Dsec_to_DmelRef_Brain_ortho_only_DF, DsecNoni_to_DmelRef_Brain_ortho_only_DF)
secsim_DmelRef_Brain_DF.list <- lapply(X = secsim_DmelRef_Brain_DF.list, FUN = DietSeurat, scale.data = T)

save(secsim_DmelRef_Brain_DF.list, file="Processed_Data/secsim_DmelRef_Brain_DF.list.RData")

resolution=0.1
npc=50
DF_sd_threshold <- 0.5

sample_list <- c("Dsim_to_DmelRef","Dsec_to_DmelRef", "DsecNoni_to_DmelRef")

ISub_sp <- function(Seurat_object, cluster,  resolution = 0.1, npc=50) {
  
  Seurat_subset <- subset(Seurat_object, idents = cluster)
  
  Seurat_subset <- SCTransform(Seurat_subset, method = "glmGamPoi", verbose = FALSE)
  Seurat_subset <- RunPCA(Seurat_subset, npcs = npc, verbose = F)
  Seurat_subset <- FindNeighbors(Seurat_subset, reduction = "pca", dims = 1:npc, verbose = F)
  Seurat_subset <- FindClusters(Seurat_subset, resolution = resolution, verbose = F)
  
  return(Seurat_subset)
  
}

ISub_DF_sp <- function(Seurat_object, n_iteration = 10, resolution = 0.1, npc=50, sample_name) {
  
  Seurat_object <- SCTransform(Seurat_object, method = "glmGamPoi", verbose = FALSE)
  Seurat_object <- RunPCA(Seurat_object, features = VariableFeatures(object = Seurat_object), verbose = F)
  
  Seurat_object <- FindNeighbors(Seurat_object, reduction = "pca", dims = 1:npc, verbose = F)
  Seurat_object <- FindClusters(Seurat_object, resolution = resolution, verbose = F)
  Seurat_object <- RunUMAP(Seurat_object, dims = 1:npc, verbose = F)
  Seurat_object <- RunTSNE(Seurat_object, dims = 1:npc, verbose = F)
  
  total_cell_number <- nrow(Seurat_object@meta.data)
  
  metadata_summary <- Seurat_object@meta.data %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(percent = n/total_cell_number*100)
  
  cluster_range_max <- max(as.numeric(as.vector(dplyr::filter(metadata_summary, percent > 0.5)$seurat_clusters)))
  
  Seurat_object_Isub_DF <- Seurat_object
  
  remain1 <- NULL
  
  for (cl in cluster_range_max:0) {
    
    name <- paste0('DF_cluster',cl)
    assign(name, ISub_sp(Seurat_object_Isub_DF, cluster = cl, resolution = 0.1))
    
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
    
    ggsave(glue::glue("Plots/Isub_DF_sp/{sample_name}_cluster{cl}_DF_range.png"), width=5, height=5)
    
    name_DF_cellid <- WhichCells(get(name), idents=doublet_clusters)
    Idents(Seurat_object_Isub_DF, cells = name_DF_cellid) <- "Doublets"
    
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
        assign(name_subA, ISub_sp(get(name), cluster = subA, resolution = 0.1))
        
        remain1 <- c(remain1, name_subA)
        
      }
      
      # subB - size가 0.1% 이상 0.5% 미만이면 해당 cluster에 대해선 subclustering 종료 및 cell_id 할당
      
      df2 <- meta_sum %>%
        dplyr::filter(percent >= 0.1 & percent <0.5)
      
      for (subB in df2$seurat_clusters) {
        
        name_subB <- paste0('cluster',cl,'_',subB)
        name_subB_cellid <- WhichCells(get(name), idents=subB)
        Idents(Seurat_object_Isub_DF, cells = name_subB_cellid) <- name_subB
        
      }
      
      # subC - size가 0.1% 미만 cluster는 모 클러스터에 소속, 모 클러스터 cell_id 할당. 추후에 cluster size에 따라 모두 re-order & re-name    
      
      df3 <- meta_sum %>%
        dplyr::filter(percent < 0.1)
      
      for (subC in df3$seurat_clusters) {
        
        name_subC_cellid <- WhichCells(get(name), idents=subC)
        Idents(Seurat_object_Isub_DF, cells = name_subC_cellid) <- cl
        
      }
      
    }
    
    else if (n_cluster_pass < 2) {
      
      # n_cluster_pass가 하나뿐인 경우, 즉 major subcluster가 안 분리되는 경우에는 전체 subcluster를 하나의 모 cluster-id로 확정 및 cell_id 할당
      
      name_cellid <- WhichCells(subset(get(name), idents= non_doublet_clusters))
      Idents(Seurat_object_Isub_DF, cells = name_cellid) <- cl
      
    }
    
    rm(list = name)
    
  }
  
  for (r in 1:n_iteration) {
    
    remainders <- paste0('remain', r)
    n_remainders <- length(get(remainders))
    remainders_new <- paste0('remain', (r+1))
    assign(remainders_new, NULL)
    
    if (length(get(remainders))>=1) {
      
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
        
        ggsave(glue::glue("Plots/Isub_DF_sp/{sample_name}_{remain_name}_DF_range.png"), width=5, height=5)
        
        name_DF_cellid <- WhichCells(remain_seurat, idents=doublet_clusters)
        Idents(Seurat_object_Isub_DF, cells = name_DF_cellid) <- "Doublets"
        
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
            assign(remainders_subA, ISub_sp(remain_seurat, cluster = subA, resolution = 0.1))
            
            assign(remainders_new, c(get(remainders_new), remainders_subA))
            
          }
          
          # subB - size가 0.1% 이상 0.5% 미만이면 해당 cluster에 대해선 subclustering 종료 및 cell_id 할당
          
          df2 <- meta_sum %>%
            dplyr::filter(percent >= 0.1 & percent <0.5)
          
          for (subB in df2$seurat_clusters) {
            
            remainders_subB <- paste0(remain_name,'_',subB)
            remainders_subB_cellid <- WhichCells(remain_seurat, idents=subB)
            Idents(Seurat_object_Isub_DF, cells = remainders_subB_cellid) <- remainders_subB
            
          }
          
          # subC - size가 0.1% 미만 cluster는 모 클러스터에 소속, 모 클러스터 cell_id 할당. 추후에 cluster size에 따라 모두 re-order & re-name    
          
          df3 <- meta_sum %>%
            dplyr::filter(percent < 0.1)
          
          for (subC in df3$seurat_clusters) {
            
            remainders_subC_cellid <- WhichCells(remain_seurat, idents=subC)
            Idents(Seurat_object_Isub_DF, cells = remainders_subC_cellid) <- remain_name
            
          }
          
        }
        
        else if (n_cluster_pass < 2) {
          
          # n_cluster_pass가 하나뿐인 경우, 즉 major subcluster가 안 분리되는 경우에는 전체 subcluster를 하나의 모 cluster-id로 확정 및 cell_id 할당
          
          remainders_cellid <- WhichCells(subset(remain_seurat, idents= non_doublet_clusters))
          Idents(Seurat_object_Isub_DF, cells = remainders_cellid) <- remain_name
          
        }
      
      
    }
    
    }
    
    
  } ## iterative subclustering for clusters > 0.5%
  
  return(Seurat_object_Isub_DF)
  
}

Dsim_to_DmelRef_ISub_DF <- ISub_DF_sp(Seurat_object=secsim_DmelRef_Brain_DF.list[[1]], n_iteration = 10, sample_name=sample_list[1])
Dsec_to_DmelRef_ISub_DF <- ISub_DF_sp(Seurat_object=secsim_DmelRef_Brain_DF.list[[2]], n_iteration = 10, sample_name=sample_list[2])
DsecNoni_to_DmelRef_ISub_DF <- ISub_DF_sp(Seurat_object=secsim_DmelRef_Brain_DF.list[[3]], n_iteration = 10, sample_name=sample_list[3])

Dsim_to_DmelRef_ISub_DF[["cluster_id"]] <- Idents(object = Dsim_to_DmelRef_ISub_DF)
Dsec_to_DmelRef_ISub_DF[["cluster_id"]] <- Idents(object = Dsec_to_DmelRef_ISub_DF)
DsecNoni_to_DmelRef_ISub_DF[["cluster_id"]] <- Idents(object = DsecNoni_to_DmelRef_ISub_DF)

secsim_DmelRef_ISub_DF_list <- list(Dsim_to_DmelRef_ISub_DF, Dsec_to_DmelRef_ISub_DF, DsecNoni_to_DmelRef_ISub_DF)

saveRDS(secsim_DmelRef_ISub_DF_list, file = "Processed_Data/secsim_DmelRef_ISub_DF_list.rds")
secsim_DmelRef_ISub_DF_list <- readRDS(file = "Processed_Data/secsim_DmelRef_ISub_DF_list.rds")


#colnames(secsim_DmelRef_ISub_DF_list[[4]]@meta.data) <- gsub("clustser_id", "cluster_id", colnames(secsim_DmelRef_ISub_DF_list[[4]]@meta.data))

## Annotation ##

secsim_DmelRef_ISub_DF_marker_list <- lapply(X = secsim_DmelRef_ISub_DF_list, FUN = FindAllMarkers, only.pos = TRUE, min.pct = 0.15, 
                                   logfc.threshold = 0.25, test.use = "MAST")

saveRDS(secsim_DmelRef_ISub_DF_marker_list, file = "Processed_Data/secsim_DmelRef_ISub_DF_marker_list.rds")
#secsim_DmelRef_ISub_DF_marker_list <- readRDS(file = "Processed_Data/secsim_DmelRef_ISub_DF_marker_list.rds")

### annotate based on markeres ###

sig_markers <- function(df_marker) {
  
  df_sigs <- df_marker %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::select(gene, cluster,p_val_adj) %>%
    tidyr::spread(key=gene, value=p_val_adj) %>%
    tidyr::gather(-cluster, key='gene', value='p_val_adj')
  
  return(df_sigs)
  
}

secsim_DmelRef_ISub_DF_sig_marker_list <- lapply(X=secsim_DmelRef_ISub_DF_marker_list, FUN=sig_markers)

dataset <- c("Dsim_to_DmelRef","Dsec_to_DmelRef","DsecNoni_to_DmelRef")

annotate_celltype_individual <- function(seurat_object,
                              cell_type, df_markers, exclude='NA', 
                              metadata) {
  
  anno_cluster <- NULL
  
  celltype_clusters <- df_markers %>%
    dplyr::filter(celltype==cell_type & !cluster %in% exclude)
  
  cluster_ids <- as.character(celltype_clusters$cluster)
  
  if(length(cluster_ids)>1) {
    
    cluster_size <- metadata %>%
      dplyr::filter(cluster_id %in% cluster_ids) %>%
      dplyr::arrange(-percent) %>%
      tibble::rownames_to_column() %>%
      dplyr::mutate(cluster_name=paste(cell_type,rowname,sep="_"))
    
    for (i in 1:nrow(cluster_size)) {
      
      Idents(seurat_object, 
             cells = WhichCells(seurat_object, idents=cluster_size$cluster_id[i])) <- cluster_size$cluster_name[i]
      
      anno_cluster <- c(anno_cluster, as.character(cluster_size$cluster_id[i]))
      
    }
    
  }
  
  else if(length(cluster_ids)==1) {
    
    Idents(seurat_object, 
           cells = WhichCells(seurat_object, idents=cluster_ids)) <- cell_type
    
    anno_cluster <- cluster_ids
    
  }
  
  return_list <- list(seurat_object, anno_cluster)
  
  return(return_list)
  
}

secsim_DmelRef_ISub_DF_list_labeled <- secsim_DmelRef_ISub_DF_list
anno_complete_list <- list()

## annotate glia ##

glial <- c("alrm", "Eaat1", "Gat", "axo","Gs2","Tret1-1","dve", "sog", "baz","CG40470")
glial_celltypes <- c("ENS","AST","PRN","SUB", "CTX")

for (i in 1:3) {
  
  seurat_anno <- secsim_DmelRef_ISub_DF_list_labeled[[i]]
  
  total_cell_number <- nrow(seurat_anno@meta.data)
  
  ISub_DF_metadata_summary <- seurat_anno@meta.data %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(percent = n/total_cell_number*100)
  
  anno_complete_list[[i]] <- NA
  
  df_glial_markers <- secsim_DmelRef_ISub_DF_sig_marker_list[[i]] %>%
    na.omit() %>%
    dplyr::filter(gene %in% c(glial)) %>%
    dplyr::mutate(gene = ifelse(gene == "Tret1-1", "Tret1_1", gene)) %>%
    tidyr::spread(gene, p_val_adj) %>%
    dplyr::mutate(celltype=ifelse(!is.na(axo) & !is.na(Gs2) & is.na(sog), "ENS",
                                  ifelse(is.na(axo) & !is.na(Gs2) & !is.na(alrm), "AST",
                                         ifelse(!is.na(Tret1_1) & !is.na(dve) & !is.na(sog), "PRN",
                                                ifelse(!is.na(Gs2) & !is.na(baz), "SUB",
                                                       ifelse(!is.na(axo) & !is.na(CG40470), "CTX", "others")))))) %>%
    dplyr::arrange(celltype)
  
  ## annotate glial clusters ##
  

  for (ct in glial_celltypes) {
    
    out_annotate_celltype <- annotate_celltype_individual(seurat_object=secsim_DmelRef_ISub_DF_list_labeled[[i]],
                                               cell_type=ct, df_markers = df_glial_markers, metadata=ISub_DF_metadata_summary)
    
    secsim_DmelRef_ISub_DF_list_labeled[[i]] <- out_annotate_celltype[[1]]
    
    anno_complete_list[[i]] <- c(anno_complete_list[[i]], out_annotate_celltype[[2]])
    
  }
  
}

## annotate KC populations ##

KC_markers <- c("crb","CG32204","Pka-C1","mamo","sNPF")
KC_celltypes <- c("αβ-KC","γ-KC","α'β'-KC")

for (i in 1:3) {
  
  seurat_anno <- secsim_DmelRef_ISub_DF_list_labeled[[i]]
  
  total_cell_number <- nrow(seurat_anno@meta.data)
  
  ISub_DF_metadata_summary <- seurat_anno@meta.data %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(percent = n/total_cell_number*100)
  
  df_KC_markers <- secsim_DmelRef_ISub_DF_sig_marker_list[[i]] %>%
    na.omit() %>%
    dplyr::filter(gene %in% c(KC_markers)) %>%
    dplyr::mutate(gene = ifelse(gene == "Pka-C1", "Pka_C1", gene)) %>%
    tidyr::spread(gene, p_val_adj) %>%
    dplyr::mutate(celltype=ifelse(!is.na(crb) & !is.na(CG32204) & !is.na(Pka_C1) & !is.na(sNPF) & is.na(mamo), "αβ-KC",
                                  ifelse(!is.na(mamo) & !is.na(CG32204) & !is.na(Pka_C1) & !is.na(sNPF) & is.na(crb), "γ-KC",
                                         ifelse(!is.na(mamo) & !is.na(CG32204) & !is.na(Pka_C1) & is.na(sNPF) & is.na(crb), "α'β'-KC", "others")))) %>%
    dplyr::arrange(celltype)
  
  ## annotate KC clusters ##
  
  
  for (ct in KC_celltypes) {
    
    out_annotate_celltype <- annotate_celltype_individual(seurat_object=secsim_DmelRef_ISub_DF_list_labeled[[i]],
                                                          cell_type=ct, df_markers = df_KC_markers, metadata=ISub_DF_metadata_summary, exclude=anno_complete_list[[i]])
    
    secsim_DmelRef_ISub_DF_list_labeled[[i]] <- out_annotate_celltype[[1]]
    
    anno_complete_list[[i]] <- c(anno_complete_list[[i]], out_annotate_celltype[[2]])
    
  }
  
}

## annotate Monoaminergic neurons ##

MA_markers <- c("Vmat", "ple", "DAT","Trhn","SerT","Tdc2","Tbh")
MA_celltypes <- c("MON","TY","SER","DOP", "OCTY")

for (i in 1:3) {
  
  seurat_anno <- secsim_DmelRef_ISub_DF_list_labeled[[i]]
  
  total_cell_number <- nrow(seurat_anno@meta.data)
  
  ISub_DF_metadata_summary <- seurat_anno@meta.data %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(percent = n/total_cell_number*100)
  
  df_MA_markers <- secsim_DmelRef_ISub_DF_sig_marker_list[[i]] %>%
    na.omit() %>%
    dplyr::filter(gene %in% c(MA_markers)) %>%
    tidyr::spread(gene, p_val_adj) %>%
    dplyr::filter(!is.na(Vmat)) %>%
    dplyr::mutate(celltype=ifelse(!is.na(Tdc2) &is.na(SerT) & is.na(Tbh), "TY",
                                  ifelse(!is.na(Tdc2) & !is.na(Tbh) & is.na(SerT), "OCTY", 
                                         ifelse(!is.na(Trhn), "SER",
                                                ifelse(!is.na(ple) & is.na(Trhn) & is.na(Tdc2), "DOP",
                                                       ifelse(!is.na(ple) | !is.na(Trhn) | !is.na(Tdc2), "MON","others")))))) %>%
    dplyr::arrange(celltype)
  
  ## annotate MA clusters ##
  
  
  for (ct in MA_celltypes) {
    
    out_annotate_celltype <- annotate_celltype_individual(seurat_object=secsim_DmelRef_ISub_DF_list_labeled[[i]],
                                                          cell_type=ct, df_markers = df_MA_markers, metadata=ISub_DF_metadata_summary, exclude=anno_complete_list[[i]])
    
    secsim_DmelRef_ISub_DF_list_labeled[[i]] <- out_annotate_celltype[[1]]
    
    anno_complete_list[[i]] <- c(anno_complete_list[[i]], out_annotate_celltype[[2]])
    
  }
  
}

## annotate Fru, clock, Poxn, OPN

fru_markers <- c('bru3','fru','tei','Ca-alpha1T','CG2269','jeb')
clock_markers <- c('tim','Dh44')
OPN_markers <- c('SiaT','otp','acj6')
Poxn_markers <- c('Poxn','CG14687')

KD_celltypes <- c("Fru","OPN","Poxn","Clock")

for (i in 1:3) {
  
  seurat_anno <- secsim_DmelRef_ISub_DF_list_labeled[[i]]
  
  total_cell_number <- nrow(seurat_anno@meta.data)
  
  ISub_DF_metadata_summary <- seurat_anno@meta.data %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(percent = n/total_cell_number*100)
  
  df_KD_markers <- secsim_DmelRef_ISub_DF_sig_marker_list[[i]] %>%
    na.omit() %>%
    na.omit() %>%
    filter(gene %in% c(fru_markers,OPN_markers,Poxn_markers,clock_markers)) %>%
    tidyr::spread(key = gene, value = p_val_adj) %>%
    rowwise() %>%  # Operate on each row individually
    mutate(celltype = case_when(
      sum(!is.na(c_across(all_of(fru_markers)))) == length(fru_markers) ~ "Fru",  # At least 3 markers must be non-NA
      sum(!is.na(c_across(all_of(OPN_markers)))) == length(OPN_markers) ~ "OPN",
      sum(!is.na(c_across(all_of(Poxn_markers)))) == length(Poxn_markers) ~ "Poxn",
      sum(!is.na(c_across(all_of(clock_markers)))) == length(clock_markers) ~ "Clock",
      TRUE ~ "others"
    )) %>%
    ungroup() %>%
    dplyr::filter(celltype != "others") %>%
    arrange(celltype) 
  
  ## annotate KD clusters ##
  
  
  for (ct in KD_celltypes) {
    
    out_annotate_celltype <- annotate_celltype_individual(seurat_object=secsim_DmelRef_ISub_DF_list_labeled[[i]],
                                                          cell_type=ct, df_markers = df_KD_markers, metadata=ISub_DF_metadata_summary, exclude=anno_complete_list[[i]])
    
    secsim_DmelRef_ISub_DF_list_labeled[[i]] <- out_annotate_celltype[[1]]
    
    anno_complete_list[[i]] <- c(anno_complete_list[[i]], out_annotate_celltype[[2]])
    
  }
  
}

## annotate based on neuropeptides ##

# load neuropeptide genes #

Dmel_gene_name <- read.table(file="genelist/Dmel_ID_genename.tsv") %>% dplyr::distinct(V1,V2) %>%
  dplyr::rename(FBID=V1, gene=V2)

df_neuropep <- data.table::fread("genelist/list_neuropeptide_FBgg0000179.txt", header = F) %>%
  dplyr::mutate(class="Neuropeptide") %>%
  dplyr::select(FBID=V1, class) %>%
  dplyr::left_join(., Dmel_gene_name, by="FBID") %>% ### few genes without matching gene_name
  na.omit()

neuropep_genes <- df_neuropep$gene

get_marker_nps <- function(seurat_object, df_marker_sigs, frac_threshold=0.05){
  
  ISub_DF_metadata_summary <- seurat_anno@meta.data %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(percent = n/total_cell_number*100)
  
  df_np_markers_specificity <- df_marker_sigs %>%
    na.omit() %>%
    dplyr::filter(gene %in% neuropep_genes) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(n_cluster=n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(n_cluster)
  
  specificity_threshold = nrow(ISub_DF_metadata_summary)*frac_threshold
  
  np_specific <- df_np_markers_specificity %>%
    dplyr::filter(n_cluster <= specificity_threshold) %>%
    dplyr::pull(gene) 
  
  return(np_specific)
  
}

for (i in 1:3) {
  
  seurat_anno <- secsim_DmelRef_ISub_DF_list_labeled[[i]]
  
  total_cell_number <- nrow(seurat_anno@meta.data)
  
  ISub_DF_metadata_summary <- seurat_anno@meta.data %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(percent = n/total_cell_number*100)
  
  np_markers <- get_marker_nps(seurat_anno, secsim_DmelRef_ISub_DF_sig_marker_list[[i]], frac_threshold=0.05)
  
  df_np_markers <- secsim_DmelRef_ISub_DF_sig_marker_list[[i]] %>%
    na.omit() %>%
    dplyr::filter(gene %in% np_markers) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(n_sp_np = n()) %>%
    dplyr::mutate(celltype=paste(gene, collapse = "_")) %>%
    dplyr::distinct(cluster, celltype)
  
  NP_celltypes <- unique(df_np_markers$celltype)
  
  ## annotate NP clusters ##
  
  for (ct in NP_celltypes) {
    
    out_annotate_celltype <- annotate_celltype_individual(seurat_object=secsim_DmelRef_ISub_DF_list_labeled[[i]],
                                                          cell_type=ct, df_markers = df_np_markers, metadata=ISub_DF_metadata_summary, exclude=anno_complete_list[[i]])
    
    secsim_DmelRef_ISub_DF_list_labeled[[i]] <- out_annotate_celltype[[1]]
    
    anno_complete_list[[i]] <- c(anno_complete_list[[i]], out_annotate_celltype[[2]])
    
  }
  
}

## annotate remianing clusters by neurotransmitter ##

NT_genes <- c("ChAT", "VAChT","VGlut", "Gad1", "VGAT")

for (i in 1:3) {
  
  seurat_anno <- secsim_DmelRef_ISub_DF_list_labeled[[i]]
  
  total_cell_number <- nrow(seurat_anno@meta.data)
  
  ISub_DF_metadata_summary <- seurat_anno@meta.data %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-n) %>%
    dplyr::mutate(percent = n/total_cell_number*100)
  
  df_NT_markers <- secsim_DmelRef_ISub_DF_sig_marker_list[[i]] %>%
    na.omit() %>%
    dplyr::filter(gene %in% c(NT_genes), !cluster %in% "Doublets") %>%
    tidyr::spread(gene, p_val_adj) %>%
    dplyr::mutate(celltype=ifelse(!is.na(ChAT) | !is.na(VAChT), "Ach",
                                  ifelse(!is.na(VGlut), "Glu",
                                         ifelse(!is.na(Gad1) | !is.na(VGAT), "GABA", "others")))) %>%
    dplyr::arrange(celltype)
  
  NT_celltypes <- unique(df_NT_markers$celltype)
  
  ## annotate NT clusters ##
  
  for (ct in NT_celltypes) {
    
    out_annotate_celltype <- annotate_celltype_individual(seurat_object=secsim_DmelRef_ISub_DF_list_labeled[[i]],
                                                          cell_type=ct, df_markers = df_NT_markers, metadata=ISub_DF_metadata_summary, exclude=anno_complete_list[[i]])
    
    secsim_DmelRef_ISub_DF_list_labeled[[i]] <- out_annotate_celltype[[1]]
    
    anno_complete_list[[i]] <- c(anno_complete_list[[i]], out_annotate_celltype[[2]])
    
  }
  
}

DimPlot(secsim_DmelRef_ISub_DF_list_labeled[[1]], label=T, reduction = 'tsne', repel=T) + NoLegend() + NoAxes()

add_cluster_label <- function(seurat_object) {
  
  seurat_object[["ClusterLabel"]] <- Idents(object = seurat_object)
  
  return(seurat_object)
  
}

secsim_DmelRef_ISub_DF_list_labeled <- lapply(X=secsim_DmelRef_ISub_DF_list_labeled, FUN=add_cluster_label)

#saveRDS(secsim_DmelRef_ISub_DF_list_labeled, file = "Processed_Data/secsim_DmelRef_ISub_DF_list_labeled.rds")
secsim_DmelRef_ISub_DF_list_labeled <- readRDS(file = "Processed_Data/secsim_DmelRef_ISub_DF_list_labeled.rds")

DimPlot(subset(secsim_DmelRef_ISub_DF_list_labeled[[1]], ClusterLabel != 'Doublets'), label=T, reduction = 'tsne', repel=T, group.by='ClusterLabel') + NoLegend() + NoAxes()
DimPlot(subset(secsim_DmelRef_ISub_DF_list_labeled[[2]], ClusterLabel != 'Doublets'), label=T, reduction = 'tsne', repel=T, group.by='ClusterLabel') + NoLegend() + NoAxes()

## cross-species cluster comparisons ##

TrioBrain.integrated_slim_ISub_DF_labeled <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_ISub_DF_labeled.rds")
Trio_ISub_DF_list_labeled <- readRDS(file = "Processed_Data/Trio_ISub_DF_list_labeled.rds")

cell_names_integrated <- rownames(TrioBrain.integrated_slim_ISub_DF_labeled@meta.data)
cell_names_Dmel <- rownames(Trio_ISub_DF_list_labeled[[1]]@meta.data)
cell_names_Dsim <- rownames(Trio_ISub_DF_list_labeled[[2]]@meta.data)
cell_names_Dsec <- rownames(Trio_ISub_DF_list_labeled[[3]]@meta.data)
cell_names_DsecNoni <- rownames(Trio_ISub_DF_list_labeled[[4]]@meta.data)
cell_names_Dsim_to_DmelRef <- rownames(secsim_DmelRef_ISub_DF_list_labeled[[1]]@meta.data)
cell_names_Dsec_to_DmelRef <- rownames(secsim_DmelRef_ISub_DF_list_labeled[[2]]@meta.data)
cell_names_DsecNoni_to_DmelRef <- rownames(secsim_DmelRef_ISub_DF_list_labeled[[3]]@meta.data)

cell_clusters_integrted <- TrioBrain.integrated_slim_ISub_DF_labeled@meta.data$ClusterLabel
cell_clusters_Dmel <- Trio_ISub_DF_list_labeled[[1]]@meta.data$ClusterLabel
cell_clusters_Dsim <- Trio_ISub_DF_list_labeled[[2]]@meta.data$ClusterLabel
cell_clusters_Dsec <- Trio_ISub_DF_list_labeled[[3]]@meta.data$ClusterLabel
cell_clusters_DsecNoni <- Trio_ISub_DF_list_labeled[[4]]@meta.data$ClusterLabel
cell_clusters_Dsim_to_DmelRef <- secsim_DmelRef_ISub_DF_list_labeled[[1]]@meta.data$ClusterLabel
cell_clusters_Dsec_to_DmelRef <- secsim_DmelRef_ISub_DF_list_labeled[[2]]@meta.data$ClusterLabel
cell_clusters_DsecNoni_to_DmelRef <- secsim_DmelRef_ISub_DF_list_labeled[[3]]@meta.data$ClusterLabel

cell_origin <- TrioBrain.integrated_slim_ISub_DF_labeled@meta.data$orig.ident
cell_type <- TrioBrain.integrated_slim_ISub_DF_labeled@meta.data$orig.ident

cell_info_integrated <- data.frame(Cell = cell_names_integrated, Cluster = cell_clusters_integrted, orig= cell_origin) %>%
  dplyr::mutate(dataset='integrated')
cell_info_Dmel <- data.frame(Cell = cell_names_Dmel, Cluster = cell_clusters_Dmel) %>%
  dplyr::mutate(dataset='Dmel')
cell_info_Dsim <- data.frame(Cell = cell_names_Dsim, Cluster = cell_clusters_Dsim) %>%
  dplyr::mutate(dataset='Dsim')
cell_info_Dsec <- data.frame(Cell = cell_names_Dsec, Cluster = cell_clusters_Dsec) %>%
  dplyr::mutate(dataset='Dsec')
cell_info_DsecNoni <- data.frame(Cell = cell_names_DsecNoni, Cluster = cell_clusters_DsecNoni) %>%
  dplyr::mutate(dataset='DsecNoni')
cell_info_Dsim_to_DmelRef <- data.frame(Cell = cell_names_Dsim_to_DmelRef, Cluster = cell_clusters_Dsim_to_DmelRef) %>%
  dplyr::mutate(dataset='Dsim_to_DmelRef')
cell_info_Dsec_to_DmelRef <- data.frame(Cell = cell_names_Dsec_to_DmelRef, Cluster = cell_clusters_Dsec_to_DmelRef) %>%
  dplyr::mutate(dataset='Dsec_to_DmelRef')
cell_info_DsecNoni_to_DmelRef <- data.frame(Cell = cell_names_DsecNoni_to_DmelRef, Cluster = cell_clusters_DsecNoni_to_DmelRef) %>%
  dplyr::mutate(dataset='DsecNoni_to_DmelRef')

cell_info_integrated$Cell <- gsub("_[1-4]", "", cell_info_integrated$Cell)
cell_info_integrated$orig <- gsub("_rep[1-6]", "", cell_info_integrated$orig)

process_cell_info <- function(data, orig_value, cell_info) {
  data %>%
    dplyr::filter(orig == orig_value) %>%
    dplyr::left_join(., cell_info, by = 'Cell') %>%
    dplyr::group_by(Cluster.x) %>%
    dplyr::mutate(n_Cluster.x = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Cluster.y) %>%
    dplyr::mutate(n_Cluster.y = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Cluster.x, Cluster.y) %>%
    dplyr::mutate(n_Cluster.xy = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(xy_x = n_Cluster.xy / n_Cluster.x,
                  xy_y = n_Cluster.xy / n_Cluster.y)
}

# Apply the function to each subset

cell_info_join_Dmel <- process_cell_info(cell_info_integrated, 'Dmel', cell_info_Dmel)
cell_info_join_Dsim <- process_cell_info(cell_info_integrated, 'Dsim', cell_info_Dsim)
cell_info_join_Dsec <- process_cell_info(cell_info_integrated, 'Dsec', cell_info_Dsec)
cell_info_join_DsecNoni <- process_cell_info(cell_info_integrated, 'DsecNoni', cell_info_DsecNoni)
cell_info_join_Dsim_to_DmelRef <- process_cell_info(cell_info_integrated, 'Dsim_to_DmelRef', cell_info_Dsim_to_DmelRef)
cell_info_join_Dsec_to_DmelRef <- process_cell_info(cell_info_integrated, 'Dsec_to_DmelRef', cell_info_Dsec_to_DmelRef)
cell_info_join_DsecNoni_to_DmelRef <- process_cell_info(cell_info_integrated, 'DsecNoni_to_DmelRef', cell_info_DsecNoni_to_DmelRef)

# Combine the results

cell_info_join_all <- dplyr::bind_rows(cell_info_join_Dmel, cell_info_join_Dsim, cell_info_join_Dsec, cell_info_join_DsecNoni,
                                       cell_info_join_Dsim_to_DmelRef, cell_info_join_Dsec_to_DmelRef, cell_info_join_DsecNoni_to_DmelRef)

save(cell_info_join_all, file="Processed_Data/cell_info_join_all.RData")




## frequency check

test_metadata <- test@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_to_DmelRef_rep1", "Dsim_to_DmelRef_rep2", "Dsim_to_DmelRef_rep3", "Dsim_to_DmelRef_rep4", "Dsim_to_DmelRef_rep5", "Dsim_to_DmelRef_rep6"), "Dsim_to_DmelRef",
                                      ifelse(orig.ident %in% c("Dsec_to_DmelRef_rep1", "Dsec_to_DmelRef_rep2", "Dsec_to_DmelRef_rep3", "Dsec_to_DmelRef_rep4", "Dsec_to_DmelRef_rep5", "Dsec_to_DmelRef_rep6"), "Dsec",
                                             "DsecNoni_to_DmelRef")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim_to_DmelRef", "Dsec","DsecNoni_to_DmelRef")))


test_metadata_summary_rep <- test_metadata %>%
  dplyr::group_by(species, CellType, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

test_metadata_summary <- test_metadata_summary_rep %>%
  dplyr::group_by(species, CellType) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

test_metadata_summary %>%
  ggplot(.) +
  geom_col(aes(y=CellType, x=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.9)) + 
  geom_errorbar(aes(y=CellType, x=percent_combined,
                    xmin=ifelse(percent_combined-sem<0,0,percent_combined-sem), xmax=percent_combined+sem, group=species), width=0.9, position=position_dodge(width=0.9)) +
  geom_point(data=test_metadata_summary_rep, size = 1, alpha =0.8,  
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


