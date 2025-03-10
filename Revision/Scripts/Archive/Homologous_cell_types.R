library(ggplot2)
library(dplyr)
library(ggsankey) #devtools::install_github("davidsjoberg/ggsankey")
library(readr)
library(ggplot2)
library(ggrepel)
library(sctransform)
library(cowplot)
library(ggsankey)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/shared_genes.RData")
load("Processed_Data/cell_info_join_Trio.RData")

TrioBrain.integrated_slim_ISub_DF_labeled <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_ISub_DF_labeled.rds")
TrioBrain.integrated_slim_ISub_DF_labeled_species <- TrioBrain.integrated_slim_ISub_DF_labeled
TrioBrain.integrated_slim_ISub_DF_labeled_species$orig.ident <- gsub("_rep[1-6]", "", TrioBrain.integrated_slim_ISub_DF_labeled_species$orig.ident)

Trio_ISub_DF_list_labeled <- readRDS(file = "Processed_Data/Trio_ISub_DF_list_labeled.rds")

species_list <- c("Dmel", "Dsim", "Dsec", "DsecNoni")
pair_list <- list(c("Dmel","Dsec"),c("Dmel","Dsim"), c("Dsim","Dsec"), c("Dsec", "DsecNoni"))

### Functions ###

add_cluster_label <- function(seurat_object) {
  
  seurat_object[["ClusterLabel"]] <- Idents(object = seurat_object)
  
  return(seurat_object)
  
}

Get_subcluster <- function(Seurat_object, cluster, method='rpca', resolution = 0.1, npc=50, genelist=shared_genes, kweight=100) {
  
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
  Seurat_integrated <- RunUMAP(Seurat_integrated, dims = 1:npc, verbose = F)
  Seurat_integrated <- RunTSNE(Seurat_integrated, dims = 1:npc, verbose = F)
  
  return(Seurat_integrated)
  
}
Get_subcluster_sp <- function(Seurat_object, cluster,  resolution = 0.1, npc=50) {
  
  Seurat_subset <- subset(Seurat_object, idents = cluster)
  
  Seurat_subset <- SCTransform(Seurat_subset, method = "glmGamPoi", verbose = FALSE)
  Seurat_subset <- RunPCA(Seurat_subset, npcs = npc, verbose = F)
  Seurat_subset <- FindNeighbors(Seurat_subset, reduction = "pca", dims = 1:npc, verbose = F)
  Seurat_subset <- FindClusters(Seurat_subset, resolution = resolution, verbose = F)
  Seurat_subset <- RunUMAP(Seurat_subset, dims = 1:npc, verbose = F)
  Seurat_subset <- RunTSNE(Seurat_subset, dims = 1:npc, verbose = F)
  
  return(Seurat_subset)
  
}
Get_subcluster_celltypes <- function(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                     individual=Trio_ISub_DF_list_labeled, celltype=MA_celltypes, res=0.1) {
  
  cluster_name_integrated <- as.character(unique(Idents(integrated)))
  
  filter_target_cluster_integrated <- as.character(matrix(cluster_name_integrated, 
                                                          ncol = length(cluster_name_integrated), 
                                                          nrow = length(celltype), 
                                                          byrow = TRUE)[t(sapply(celltype, function(ct) grepl(ct, cluster_name_integrated)))])
  
  subcluster_list <- list()
  
  for (i in 1:4) {
    
    cluster_name <- as.character(unique(Idents(individual[[i]])))
    
    filter_target_cluster <- unique(as.character(matrix(cluster_name, 
                                                        ncol = length(cluster_name), 
                                                        nrow = length(celltype), 
                                                        byrow = TRUE)[t(sapply(celltype, function(ct) grepl(ct, cluster_name)))]))
    
    subcluster_list[[i]] <- Get_subcluster_sp(individual[[i]], cluster=filter_target_cluster, resolution = res)
    
  }
  
  
  subcluster_list[[5]] <- Get_subcluster(integrated, cluster=filter_target_cluster_integrated, resolution = res)
  
  return(subcluster_list)
  
}
Get_subcluster_celltypes_from_Int <- function(integrated=TrioBrain.integrated_slim_ISub_DF_labeled,
                                              celltype=MA_celltypes, res=0.1) {
  
  cluster_name_integrated <- as.character(unique(Idents(integrated)))
  
  filter_target_cluster_integrated <- as.character(matrix(cluster_name_integrated, 
                                                          ncol = length(cluster_name_integrated), 
                                                          nrow = length(celltype), 
                                                          byrow = TRUE)[t(sapply(celltype, function(ct) grepl(ct, cluster_name_integrated)))])
  
  subcluster_list <- list()
  
  for (i in 1:4) {
    
    sp_rep <- unique(integrated$orig.ident)[grepl(species_list[i], unique(integrated$orig.ident))]
    
    subcluster_list[[i]] <- Get_subcluster_sp(subset(integrated, orig.ident %in% sp_rep), cluster=filter_target_cluster_integrated, resolution = res)
    
  }
  
  
  subcluster_list[[5]] <- Get_subcluster(integrated, cluster=filter_target_cluster_integrated, resolution = res)
  
  return(subcluster_list)
  
}



process_subset_cell_info <- function(data, orig_value, subset_cell_info) {
  data %>%
    dplyr::filter(orig == orig_value) %>%
    dplyr::left_join(., subset_cell_info, by = 'Cell') %>%
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

generate_sankey_celltype <- function(df=cell_info_join_Trio, target_ct=glial_celltypes, target_sp=c("Dmel","Dsec"), filter_threshold =0.05) {
  
  df_target_prep <- df %>%
    filter(orig %in% target_sp) %>%
    filter(
      sapply(target_ct, function(ct) grepl(ct, Cluster.x)) %>% rowSums() > 0 |
        sapply(target_ct, function(ct) grepl(ct, Cluster.y)) %>% rowSums() > 0
    ) %>%
    filter(xy_x>=filter_threshold & xy_y>=filter_threshold) %>%
    dplyr::select(Cell, Cluster.x, Cluster.y, orig) %>%
    tidyr::spread(key='orig', value='Cluster.y') %>%
    dplyr::select(sp1=target_sp[1], Integrated=Cluster.x, sp2=target_sp[2])
  
  df_target_cluster_sp1 <- as.character(na.omit(unique(df_target_prep$sp1)))
  df_target_cluster_sp2 <- as.character(na.omit(unique(df_target_prep$sp2)))
  df_target_cluster_Integrated <- as.character(na.omit(unique(df_target_prep$Integrated)))
  
  df_target <- df %>%
    filter(orig %in% target_sp) %>%
    dplyr::select(Cell, Cluster.x, Cluster.y, orig) %>%
    tidyr::spread(key='orig', value='Cluster.y') %>%
    dplyr::select(sp1=target_sp[1], Integrated=Cluster.x, sp2=target_sp[2]) %>%
    #dplyr::filter((sp1 %in% df_target_cluster_sp1 & Integrated %in% df_target_cluster_Integrated)
    #              | (sp2 %in% df_target_cluster_sp2& Integrated %in% df_target_cluster_Integrated)) %>%
    dplyr::filter((sp1 %in% df_target_cluster_sp1)
                  | (sp2 %in% df_target_cluster_sp2)) 
  
  
  filter_cluster_sp1 <- as.character(matrix(df_target_cluster_sp1, 
                                            ncol = length(df_target_cluster_sp1), 
                                            nrow = length(target_ct), 
                                            byrow = TRUE)[t(sapply(target_ct, function(ct) grepl(ct, df_target_cluster_sp1)))])
  
  filter_cluster_Integrated <- as.character(matrix(df_target_cluster_Integrated, 
                                                   ncol = length(df_target_cluster_Integrated), 
                                                   nrow = length(target_ct), 
                                                   byrow = TRUE)[t(sapply(target_ct, function(ct) grepl(ct, df_target_cluster_Integrated)))])
  
  df_plot_target <- df_target %>%
    make_long(sp1,Integrated, 
              sp2) %>%
    dplyr::filter(if_any(c(node, next_x, next_node), ~ !is.na(.))) %>%
    dplyr::mutate(next_node=ifelse((x=='sp1' & next_x=='Integrated' & !next_node %in% filter_cluster_Integrated) | 
                                     x=='Integrated' & next_x=='sp2' & !node %in% filter_cluster_Integrated,
                                   'others', next_node)) %>%
    dplyr::filter(next_node != 'others' | is.na(next_node)) %>%
    dplyr::filter(x!='sp1' | !is.na(node)) 
  
  df_plot_target$x <- factor(df_plot_target$x, levels = c('sp1','Integrated','sp2'), labels=c(target_sp[1],'Integrated',target_sp[2]))
  df_plot_target$next_x <- factor(df_plot_target$next_x, levels = c('sp1','Integrated','sp2'), labels=c(target_sp[1],'Integrated',target_sp[2]))
  
  pl_target <- ggplot(df_plot_target, aes(x = x,                        
                                          next_x = next_x,                                     
                                          node = node,
                                          next_node = next_node,        
                                          fill = factor(node))) +
    geom_sankey(flow.alpha = 0.5,          #This Creates the transparency of your node 
                node.color = 'black',     # This is your node color  - set NA for transparency       
                space = 0,
                show.legend = F,
                na.rm=F)  +
    geom_sankey_text(aes(label = node), 
                     space = 0,
                     size = 3,
                     show.legend = F,
                     na.rm=F) +  # Adjust size and vjust as needed
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          panel.grid = element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_text(size=12, color='black'))
  
  pl_target
  
  return(pl_target)
  
}
sankey_plots_trio <- function(list_celltypes, celltypes, pairlist=pair_list) {
  
  list_celltypes <- lapply(X=list_celltypes, FUN=add_cluster_label)
  
  cell_names_integrated <- rownames(list_celltypes[[5]]@meta.data)
  cell_names_Dmel <- rownames(list_celltypes[[1]]@meta.data)
  cell_names_Dsim <- rownames(list_celltypes[[2]]@meta.data)
  cell_names_Dsec <- rownames(list_celltypes[[3]]@meta.data)
  
  cell_clusters_integrted <- list_celltypes[[5]]@meta.data$ClusterLabel
  cell_clusters_Dmel <- list_celltypes[[1]]@meta.data$ClusterLabel
  cell_clusters_Dsim <- list_celltypes[[2]]@meta.data$ClusterLabel
  cell_clusters_Dsec <- list_celltypes[[3]]@meta.data$ClusterLabel
  
  cell_origin <- list_celltypes[[5]]@meta.data$orig.ident
  cell_type <- list_celltypes[[5]]@meta.data$orig.ident
  
  cell_info_integrated <- data.frame(Cell = cell_names_integrated, Cluster = cell_clusters_integrted, orig= cell_origin) %>%
    dplyr::mutate(dataset='integrated')
  cell_info_Dmel <- data.frame(Cell = cell_names_Dmel, Cluster = cell_clusters_Dmel) %>%
    dplyr::mutate(dataset='Dmel')
  cell_info_Dsim <- data.frame(Cell = cell_names_Dsim, Cluster = cell_clusters_Dsim) %>%
    dplyr::mutate(dataset='Dsim')
  cell_info_Dsec <- data.frame(Cell = cell_names_Dsec, Cluster = cell_clusters_Dsec) %>%
    dplyr::mutate(dataset='Dsec')
  
  cell_info_integrated$Cell <- gsub("_[1-4]", "", cell_info_integrated$Cell)
  cell_info_integrated$orig <- gsub("_rep[1-6]", "", cell_info_integrated$orig)
  
  # Apply the function to each subset
  cell_info_join_Dmel <- process_subset_cell_info(cell_info_integrated, 'Dmel', cell_info_Dmel)
  cell_info_join_Dsim <- process_subset_cell_info(cell_info_integrated, 'Dsim', cell_info_Dsim)
  cell_info_join_Dsec <- process_subset_cell_info(cell_info_integrated, 'Dsec', cell_info_Dsec)
  
  # Combine the results
  cell_info_join_Trio <- dplyr::bind_rows(cell_info_join_Dmel, cell_info_join_Dsim, cell_info_join_Dsec)
  
  plot_list <- list()
  
  for (i in 1:3) {
    
    pair <- pairlist[[i]]
    
    sp1 <- pair[1]
    sp2 <- pair[2]
    
    plot_list[[i]] <- generate_sankey_celltype(df=cell_info_join_Trio, target_ct=celltypes, target_sp=pair, filter_threshold =0)
    
  }
  
  return(plot_list)
  
}
plot_sp_dimplot_merge <- function(plot_list) {
  
  dimplot_list <- list()
  
  for (i in 1:3) {
    
    dimplot_list[[i]] <- DimPlot(plot_list[[i]], label=T) + NoAxes() + ggtitle(species_list[[i]]) + NoLegend()
    
  }
  
  dimplot_list[[5]] <- DimPlot(plot_list[[5]], label=T) + NoAxes() + ggtitle('Integrated') + NoLegend()
  
  return(dimplot_list)
  
}
plot_sp_dotplot_merge <- function(plot_list, markers) {
  
  dotplot_list <- list()
  
  for (i in 1:3) {
    
    dotplot_list[[i]] <- DotPlot(plot_list[[i]], features = markers) + ggtitle(species_list[[i]]) + 
      theme(axis.title = element_blank())
  }
  
  dotplot_list[[5]] <- DotPlot(plot_list[[5]], features = markers) +  ggtitle('Integrated') + 
    theme(axis.title = element_blank())
  
  return(dotplot_list)
  
}
plot_cell_type_groups <- function(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                  individual=Trio_ISub_DF_list_labeled, celltype=glial_celltypes,
                                  color='darkorange') {
  
  plot_list <- list()
  
  cluster_name <- as.character(unique(Idents(integrated)))
  
  filter_target_cluster <- as.character(matrix(cluster_name, 
                                               ncol = length(cluster_name), 
                                               nrow = length(celltype), 
                                               byrow = TRUE)[t(sapply(celltype, function(ct) grepl(ct, cluster_name)))])
  
  highlight_cells <- WhichCells(integrated, idents=filter_target_cluster)
  
  for (i in 1:3) {
    
    highlight_cells_id_sp <- WhichCells(subset(integrated, orig.ident == species_list[i]), idents=filter_target_cluster)
    highlight_cells_id_sp <- gsub("_[1-4]", "", highlight_cells_id_sp)
    
    plot_list[[i]] <- DimPlot(individual[[i]], label=F,shuffle=T,
                              cells.highlight= list(highlight=highlight_cells_id_sp), 
                              reduction = "tsne", 
                              pt.size = 0.01,sizes.highlight = 0.01,
                              cols.highlight = color, cols= "grey40", raster=FALSE) + NoLegend() + NoAxes() + ggtitle(species_list[i])
    
  }
  
  plot_list[[4]] <- DimPlot(integrated, label=F,shuffle=T,
                            cells.highlight= list(highlight=highlight_cells), 
                            reduction = "tsne", 
                            pt.size = 0.01,sizes.highlight = 0.01,
                            cols.highlight = color, cols= "grey40", raster=FALSE) + NoLegend() + NoAxes() + ggtitle('Integrated')
  
  return(plot_list)
  
}
plot_grid_merge <- function(plot_list, int=5) {
  
  plot_out <- cowplot::plot_grid(plot_list[[int]], plot_list[[1]],plot_list[[2]],plot_list[[3]], ncol=2)
  
  return(plot_out)
  
}

### Functions end ###

## celltype group plots ###

glial_celltypes <- c("ENS","AST","PRN","SUB", "CTX")
glial_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                     individual=Trio_ISub_DF_list_labeled, celltype=glial_celltypes,
                                     color='darkorange')
plot_grid_merge(glial_plots)
ggsave("Plots/celltype_groups/glial_groups.png", width=10, height=11)

SUB_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                     individual=Trio_ISub_DF_list_labeled, celltype='SUB',
                                     color='purple')
plot_grid_merge(SUB_plots)
ggsave("Plots/celltype_groups/SUBs.png", width=10, height=11)

 ## manual inspection ##

list_subcluster_celltypes_glia <- Get_subcluster_celltypes(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                         individual=Trio_ISub_DF_list_labeled, celltype=glial_celltypes, res=0.5)

Glia_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_glia[[5]])

dimplot_list_glia <- plot_sp_dimplot_merge(list_subcluster_celltypes_glia)
plot_grid_merge(dimplot_list_glia)
ggsave("Plots/celltype_groups/glia_sp_subtypes.png", width=10, height=11)

glial_markers <- c("alrm", "Eaat1", "Gat", "axo","Gs2", "scro", "Prip", "Tret1-1","dve", "sog", "baz","CG40470")

dotplot_list_glia <- plot_sp_dotplot_merge(list_subcluster_celltypes_glia, glial_markers)
plot_grid_merge(dotplot_list_glia)
ggsave("Plots/celltype_groups/glia_sp_subtypes_dotplots.png", width=18, height=10)

## Re-annotate Glial clusters ##

## Dmel - AST : 1, ENS_1 : 0, ENS_2 : 4, PRN : 2, CTX : 3, SUB : 5

AST_Dmel <- WhichCells(list_subcluster_celltypes_glia[[1]], idents=1)
Idents(list_subcluster_celltypes_glia[[1]], cells = AST_Dmel) <- 'AST'
ENS_1_Dmel <- WhichCells(list_subcluster_celltypes_glia[[1]], idents=0)
Idents(list_subcluster_celltypes_glia[[1]], cells = ENS_1_Dmel ) <- 'ENS_1'
ENS_2_Dmel <- WhichCells(list_subcluster_celltypes_glia[[1]], idents=4)
Idents(list_subcluster_celltypes_glia[[1]], cells = ENS_2_Dmel ) <- 'ENS_2'
PRN_Dmel <- WhichCells(list_subcluster_celltypes_glia[[1]], idents=2)
Idents(list_subcluster_celltypes_glia[[1]], cells = PRN_Dmel) <- 'PRN'
CTX_Dmel <- WhichCells(list_subcluster_celltypes_glia[[1]], idents=3)
Idents(list_subcluster_celltypes_glia[[1]], cells = CTX_Dmel) <- 'CTX'
SUB_Dmel <- WhichCells(list_subcluster_celltypes_glia[[1]], idents=5)
Idents(list_subcluster_celltypes_glia[[1]], cells = SUB_Dmel) <- 'SUB'

DimPlot(list_subcluster_celltypes_glia[[1]], label=T) + NoLegend() + NoAxes()

## Dsim - AST : 0,6, ENS_1 : 1,5, ENS_2 : 3, PRN : 2, CTX : 4

AST_Dsim <- WhichCells(list_subcluster_celltypes_glia[[2]], idents=c(0,6))
Idents(list_subcluster_celltypes_glia[[2]], cells = AST_Dsim) <- 'AST'
ENS_1_Dsim <- WhichCells(list_subcluster_celltypes_glia[[2]], idents=c(1,5))
Idents(list_subcluster_celltypes_glia[[2]], cells = ENS_1_Dsim ) <- 'ENS_1'
ENS_2_Dsim <- WhichCells(list_subcluster_celltypes_glia[[2]], idents=3)
Idents(list_subcluster_celltypes_glia[[2]], cells = ENS_2_Dsim ) <- 'ENS_2'
PRN_Dsim <- WhichCells(list_subcluster_celltypes_glia[[2]], idents=2)
Idents(list_subcluster_celltypes_glia[[2]], cells = PRN_Dsim) <- 'PRN'
CTX_Dsim <- WhichCells(list_subcluster_celltypes_glia[[2]], idents=4)
Idents(list_subcluster_celltypes_glia[[2]], cells = CTX_Dsim) <- 'CTX'

DimPlot(list_subcluster_celltypes_glia[[2]], label=T) + NoLegend() + NoAxes()

## Dsec - AST : 1,6,7, ENS_1 : 0, ENS_2 : 3, PRN : 2, CTX : 4, SUB : 5

AST_Dsec <- WhichCells(list_subcluster_celltypes_glia[[3]], idents=c(1,6,7))
Idents(list_subcluster_celltypes_glia[[3]], cells = AST_Dsec) <- 'AST'
ENS_1_Dsec <- WhichCells(list_subcluster_celltypes_glia[[3]], idents=0)
Idents(list_subcluster_celltypes_glia[[3]], cells = ENS_1_Dsec ) <- 'ENS_1'
ENS_2_Dsec <- WhichCells(list_subcluster_celltypes_glia[[3]], idents=3)
Idents(list_subcluster_celltypes_glia[[3]], cells = ENS_2_Dsec ) <- 'ENS_2'
PRN_Dsec <- WhichCells(list_subcluster_celltypes_glia[[3]], idents=2)
Idents(list_subcluster_celltypes_glia[[3]], cells = PRN_Dsec) <- 'PRN'
CTX_Dsec <- WhichCells(list_subcluster_celltypes_glia[[3]], idents=4)
Idents(list_subcluster_celltypes_glia[[3]], cells = CTX_Dsec) <- 'CTX'
SUB_Dsec <- WhichCells(list_subcluster_celltypes_glia[[3]], idents=5)
Idents(list_subcluster_celltypes_glia[[3]], cells = SUB_Dsec) <- 'SUB'

DimPlot(list_subcluster_celltypes_glia[[3]], label=T) + NoLegend() + NoAxes()

## Integrted - AST : 0,8,10,13, ENS_1 : 1,2,6, ENS_2 : 4,  ENS_3 : 7, PRN : 3, CTX : 5, SUB : 9

AST_Int <- WhichCells(list_subcluster_celltypes_glia[[5]], idents=c(0,8,10,13))
Idents(list_subcluster_celltypes_glia[[5]], cells = AST_Int) <- 'AST'
ENS_1_Int <- WhichCells(list_subcluster_celltypes_glia[[5]], idents=c(1,2,6))
Idents(list_subcluster_celltypes_glia[[5]], cells = ENS_1_Int ) <- 'ENS_1'
ENS_2_Int <- WhichCells(list_subcluster_celltypes_glia[[5]], idents=4)
Idents(list_subcluster_celltypes_glia[[5]], cells = ENS_2_Int ) <- 'ENS_2'
ENS_3_Int <- WhichCells(list_subcluster_celltypes_glia[[5]], idents=7)
Idents(list_subcluster_celltypes_glia[[5]], cells = ENS_3_Int ) <- 'ENS_3'
PRN_Int <- WhichCells(list_subcluster_celltypes_glia[[5]], idents=3)
Idents(list_subcluster_celltypes_glia[[5]], cells = PRN_Int) <- 'PRN'
CTX_Int <- WhichCells(list_subcluster_celltypes_glia[[5]], idents=5)
Idents(list_subcluster_celltypes_glia[[5]], cells = CTX_Int) <- 'CTX'
SUB_Int <- WhichCells(list_subcluster_celltypes_glia[[5]], idents=9)
Idents(list_subcluster_celltypes_glia[[5]], cells = SUB_Int) <- 'SUB'

## Merge ENS ##

ENS_Dmel <- WhichCells(list_subcluster_celltypes_glia[[1]], idents=c('ENS_1','ENS_2'))
Idents(list_subcluster_celltypes_glia[[1]], cells = ENS_Dmel ) <- 'ENS'

ENS_Dsim <- WhichCells(list_subcluster_celltypes_glia[[2]], idents=c('ENS_1','ENS_2'))
Idents(list_subcluster_celltypes_glia[[2]], cells = ENS_Dsim ) <- 'ENS'

ENS_Dsec <- WhichCells(list_subcluster_celltypes_glia[[3]], idents=c('ENS_1','ENS_2'))
Idents(list_subcluster_celltypes_glia[[3]], cells = ENS_Dsec ) <- 'ENS'

ENS_Int <- WhichCells(list_subcluster_celltypes_glia[[5]], idents=c('ENS_1','ENS_2','ENS_3'))
Idents(list_subcluster_celltypes_glia[[5]], cells = ENS_Int ) <- 'ENS'

DimPlot(list_subcluster_celltypes_glia[[5]], label=T) + NoLegend() + NoAxes()

### Final Plot ###

dimplot_list_glia_annotated <- plot_sp_dimplot_merge(list_subcluster_celltypes_glia)
plot_grid_merge(dimplot_list_glia_annotated)
ggsave("Plots/celltype_groups/glia_sp_subtypes_annotated.png", width=10, height=11)

glia_simple_marker <- c("Eaat1", "axo", "Tret1-1", "CG40470", "sog","baz")

dotplot_list_glia <- plot_sp_dotplot_merge(list_subcluster_celltypes_glia, glia_simple_marker)
plot_grid_merge(dotplot_list_glia)
ggsave("Plots/celltype_groups/glia_sp_subtypes_annotated_dotplots_final.png", width=15, height=10)

## Generate saneky plot ##

glial_cell_sankey_final <- sankey_plots_trio(list_celltypes=list_subcluster_celltypes_glia, celltypes = glial_celltypes)

cowplot::plot_grid(plotlist=glial_cell_sankey_final, nrow=1)

ggsave("Plots/celltype_groups/glia_sp_subtypes_annotated_sankey_final.png", width=20, height=8)

## KC manual insepction ##

KC_celltypes <- c("αβ-KC","γ-KC","α'β'-KC")
KC_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                   individual=Trio_ISub_DF_list_labeled, celltype=KC_celltypes,
                                   color='red')
plot_grid_merge(KC_plots)
ggsave("Plots/celltype_groups/KC_groups.png", width=10, height=11)


list_subcluster_celltypes_KC <- Get_subcluster_celltypes(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                           individual=Trio_ISub_DF_list_labeled, celltype=KC_celltypes, res=0.2)

KC_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_KC[[5]])

dimplot_list_KC <- plot_sp_dimplot_merge(list_subcluster_celltypes_KC)
plot_grid_merge(dimplot_list_KC)
ggsave("Plots/celltype_groups/KC_sp_subtypes.png", width=10, height=11)

KC_markers <- c("crb","CG32204","Pka-C1","mamo","sNPF", "rn")

dotplot_list_KC <- plot_sp_dotplot_merge(list_subcluster_celltypes_KC, KC_markers)
plot_grid_merge(dotplot_list_KC)
ggsave("Plots/celltype_groups/KC_sp_subtypes_dotplots.png", width=15, height=8)

## Re-annotate KC clusters ##

## Dmel - αβ-KC_1 : 0, αβ-KC_2 : 4, γ-KC : 1,3, α'β'-KC : 2 

abKC_1_Dmel <- WhichCells(list_subcluster_celltypes_KC[[1]], idents=0)
Idents(list_subcluster_celltypes_KC[[1]], cells = abKC_1_Dmel) <- 'αβ-KC_1'
abKC_2_Dmel <- WhichCells(list_subcluster_celltypes_KC[[1]], idents=4)
Idents(list_subcluster_celltypes_KC[[1]], cells = abKC_2_Dmel) <- 'αβ-KC_2'
rKC_Dmel <- WhichCells(list_subcluster_celltypes_KC[[1]], idents=c(1,3))
Idents(list_subcluster_celltypes_KC[[1]], cells = rKC_Dmel) <- 'γ-KC'
apbpKC_Dmel <- WhichCells(list_subcluster_celltypes_KC[[1]], idents=2)
Idents(list_subcluster_celltypes_KC[[1]], cells = apbpKC_Dmel) <- "α'β'-KC"

DimPlot(list_subcluster_celltypes_KC[[1]], label=T) + NoLegend() + NoAxes()

## Dsim - αβ-KC_1 : 0, αβ-KC_2 : 3, γ-KC : 1,  α'β'-KC : 2 

abKC_1_Dsim <- WhichCells(list_subcluster_celltypes_KC[[2]], idents=0)
Idents(list_subcluster_celltypes_KC[[2]], cells = abKC_1_Dsim) <- 'αβ-KC_1'
abKC_2_Dsim <- WhichCells(list_subcluster_celltypes_KC[[2]], idents=3)
Idents(list_subcluster_celltypes_KC[[2]], cells = abKC_2_Dsim) <- 'αβ-KC_2'
rKC_Dsim <- WhichCells(list_subcluster_celltypes_KC[[2]], idents=1)
Idents(list_subcluster_celltypes_KC[[2]], cells = rKC_Dsim) <- 'γ-KC'
apbpKC_Dsim <- WhichCells(list_subcluster_celltypes_KC[[2]], idents=2)
Idents(list_subcluster_celltypes_KC[[2]], cells = apbpKC_Dsim) <- "α'β'-KC"

DimPlot(list_subcluster_celltypes_KC[[2]], label=T) + NoLegend() + NoAxes()

## Dsec - αβ-KC_1 : 0, αβ-KC_2 : 4, γ-KC : 1,3,  α'β'-KC : 2 

abKC_1_Dsec <- WhichCells(list_subcluster_celltypes_KC[[3]], idents=0)
Idents(list_subcluster_celltypes_KC[[3]], cells = abKC_1_Dsec) <- 'αβ-KC_1'
abKC_2_Dsec <- WhichCells(list_subcluster_celltypes_KC[[3]], idents=4)
Idents(list_subcluster_celltypes_KC[[3]], cells = abKC_2_Dsec) <- 'αβ-KC_2'
rKC_Dsec <- WhichCells(list_subcluster_celltypes_KC[[3]], idents=c(1,3))
Idents(list_subcluster_celltypes_KC[[3]], cells = rKC_Dsec) <- 'γ-KC'
apbpKC_Dsec <- WhichCells(list_subcluster_celltypes_KC[[3]], idents=2)
Idents(list_subcluster_celltypes_KC[[3]], cells = apbpKC_Dsec) <- "α'β'-KC"

DimPlot(list_subcluster_celltypes_KC[[3]], label=T) + NoLegend() + NoAxes()

## Integrted - αβ-KC_1 : 0,4, αβ-KC_2 : 5, γ-KC : 1,3,  α'β'-KC : 2 

abKC_1_Int <- WhichCells(list_subcluster_celltypes_KC[[5]], idents=c(0,4))
Idents(list_subcluster_celltypes_KC[[5]], cells = abKC_1_Int) <- 'αβ-KC_1'
abKC_2_Int <- WhichCells(list_subcluster_celltypes_KC[[5]], idents=5)
Idents(list_subcluster_celltypes_KC[[5]], cells = abKC_2_Int) <- 'αβ-KC_2'
rKC_Int <- WhichCells(list_subcluster_celltypes_KC[[5]], idents=c(1,3))
Idents(list_subcluster_celltypes_KC[[5]], cells = rKC_Int) <- 'γ-KC'
apbpKC_Int <- WhichCells(list_subcluster_celltypes_KC[[5]], idents=2)
Idents(list_subcluster_celltypes_KC[[5]], cells = apbpKC_Int) <- "α'β'-KC"

DimPlot(list_subcluster_celltypes_KC[[5]], label=T) + NoLegend() + NoAxes()


### Final Plot ###

dimplot_list_KC_annotated <- plot_sp_dimplot_merge(list_subcluster_celltypes_KC)
plot_grid_merge(dimplot_list_KC_annotated)
ggsave("Plots/celltype_groups/KC_sp_subtypes_annotated.png", width=10, height=11)

dotplot_list_KC <- plot_sp_dotplot_merge(list_subcluster_celltypes_KC, KC_markers)
plot_grid_merge(dotplot_list_KC)
ggsave("Plots/celltype_groups/KC_sp_subtypes_annotated_dotplots_final.png", width=15, height=10)

KC_cell_sankey_final <- sankey_plots_trio(list_celltypes=list_subcluster_celltypes_KC, celltypes = KC_celltypes)

cowplot::plot_grid(plotlist=KC_cell_sankey_final, nrow=1)

ggsave("Plots/celltype_groups/KC_sp_subtypes_annotated_sankey_final.png", width=20, height=8)




### Manual inspection of MA cells ###

MA_celltypes <- c("MON","TY","SER","DOP", "OCTY")
MA_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                  individual=Trio_ISub_DF_list_labeled, celltype=MA_celltypes,
                                  color='turquoise')
plot_grid_merge(MA_plots)
ggsave("Plots/celltype_groups/MA_groups.png", width=10, height=11)

list_subcluster_celltypes_MA <- Get_subcluster_celltypes(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                 individual=Trio_ISub_DF_list_labeled, celltype=MA_celltypes)

MA_markers <- c("Vmat", "ple", "DAT","Trhn","SerT","Tdc2","Tbh")

dimplot_list_ma <- plot_sp_dimplot_merge(list_subcluster_celltypes_MA)
plot_grid_merge(dimplot_list_ma)
ggsave("Plots/celltype_groups/MA_sp_subtypes.png", width=10, height=11)

dotplot_list_ma <- plot_sp_dotplot_merge(list_subcluster_celltypes_MA, MA_markers)
plot_grid_merge(dotplot_list_ma)
ggsave("Plots/celltype_groups/MA_sp_subtypes_dotplots.png", width=15, height=10)

## Re-annotate MON clusters ##

 ## Dmel - 0 : DOP/SER, 1 : DOP, 2 : TY_1, 3: OCTY, 4: TY_2, 5: Tbh, 6:SER

DOP_SER_Dmel <- WhichCells(list_subcluster_celltypes_MA[[1]], idents=0)
Idents(list_subcluster_celltypes_MA[[1]], cells = DOP_SER_Dmel) <- 'DOP_SER'

DOP_Dmel <- WhichCells(list_subcluster_celltypes_MA[[1]], idents=1)
Idents(list_subcluster_celltypes_MA[[1]], cells = DOP_Dmel) <- 'DOP'

TY_1_Dmel <- WhichCells(list_subcluster_celltypes_MA[[1]], idents=2)
Idents(list_subcluster_celltypes_MA[[1]], cells = TY_1_Dmel) <- 'TY_1'

OCTY_Dmel <- WhichCells(list_subcluster_celltypes_MA[[1]], idents=3)
Idents(list_subcluster_celltypes_MA[[1]], cells = OCTY_Dmel) <- 'OCTY'

TY_2_Dmel <- WhichCells(list_subcluster_celltypes_MA[[1]], idents=4)
Idents(list_subcluster_celltypes_MA[[1]], cells = TY_2_Dmel) <- 'TY_2'

Tbh_Dmel <- WhichCells(list_subcluster_celltypes_MA[[1]], idents=5)
Idents(list_subcluster_celltypes_MA[[1]], cells = Tbh_Dmel) <- 'Tbh'

SER_Dmel <- WhichCells(list_subcluster_celltypes_MA[[1]], idents=6)
Idents(list_subcluster_celltypes_MA[[1]], cells = SER_Dmel) <- 'SER'

DimPlot(list_subcluster_celltypes_MA[[1]], label=T) + NoLegend() + NoAxes()

 ## Dsim - 0 : DOP/SER, 1 : DOP, 2 : TY_1, 3: OCTY, 4: TY_2

DOP_SER_Dsim <- WhichCells(list_subcluster_celltypes_MA[[2]], idents=0)
Idents(list_subcluster_celltypes_MA[[2]], cells = DOP_SER_Dsim) <- 'DOP_SER'

DOP_Dsim <- WhichCells(list_subcluster_celltypes_MA[[2]], idents=1)
Idents(list_subcluster_celltypes_MA[[2]], cells = DOP_Dsim) <- 'DOP'

TY_1_Dsim <- WhichCells(list_subcluster_celltypes_MA[[2]], idents=2)
Idents(list_subcluster_celltypes_MA[[2]], cells = TY_1_Dsim) <- 'TY_1'

OCTY_Dsim <- WhichCells(list_subcluster_celltypes_MA[[2]], idents=3)
Idents(list_subcluster_celltypes_MA[[2]], cells = OCTY_Dsim) <- 'OCTY'

TY_2_Dsim <- WhichCells(list_subcluster_celltypes_MA[[2]], idents=4)
Idents(list_subcluster_celltypes_MA[[2]], cells = TY_2_Dsim) <- 'TY_2'

DimPlot(list_subcluster_celltypes_MA[[2]], label=T) + NoLegend() + NoAxes()

## Dsec - 0 : DOP/SER, 1 : DOP, 2 : TY_1, 3: OCTY, 4: TY_2

DOP_SER_Dsec <- WhichCells(list_subcluster_celltypes_MA[[3]], idents=0)
Idents(list_subcluster_celltypes_MA[[3]], cells = DOP_SER_Dsec) <- 'DOP_SER'

DOP_Dsec <- WhichCells(list_subcluster_celltypes_MA[[3]], idents=1)
Idents(list_subcluster_celltypes_MA[[3]], cells = DOP_Dsec) <- 'DOP'

TY_1_Dsec <- WhichCells(list_subcluster_celltypes_MA[[3]], idents=2)
Idents(list_subcluster_celltypes_MA[[3]], cells = TY_1_Dsec) <- 'TY_1'

OCTY_Dsec <- WhichCells(list_subcluster_celltypes_MA[[3]], idents=3)
Idents(list_subcluster_celltypes_MA[[3]], cells = OCTY_Dsec) <- 'OCTY'

TY_2_Dsec <- WhichCells(list_subcluster_celltypes_MA[[3]], idents=4)
Idents(list_subcluster_celltypes_MA[[3]], cells = TY_2_Dsec) <- 'TY_2'

DimPlot(list_subcluster_celltypes_MA[[3]], label=T) + NoLegend() + NoAxes()

## Integrated - 0 : DOP_1, 1: DOP/SER_1, 2: DOP_2, 3: Doublet_1, 4: TY_1, 5: TY_2, 6: Doublet_2, 7: DOP_3, 8: OCTY, 
##               9: Doublet_3, 10: SER, 11: DOP/SER_2, 12: Tbh_1, 13: Tbh_2, 14 : DOP/SER_3, 15: DOP_4, 16: SER_2

DOP_1_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=0)
Idents(list_subcluster_celltypes_MA[[5]], cells = DOP_1_Int) <- 'DOP_1'
DOP_2_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=2)
Idents(list_subcluster_celltypes_MA[[5]], cells = DOP_2_Int) <- 'DOP_2'
DOP_3_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=7)
Idents(list_subcluster_celltypes_MA[[5]], cells = DOP_3_Int) <- 'DOP_3'
DOP_4_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=15)
Idents(list_subcluster_celltypes_MA[[5]], cells = DOP_4_Int) <- 'DOP_4'

DOP_SER_1_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=1)
Idents(list_subcluster_celltypes_MA[[5]], cells = DOP_SER_1_Int) <- 'DOP_SER_1'
DOP_SER_2_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=11)
Idents(list_subcluster_celltypes_MA[[5]], cells = DOP_SER_2_Int) <- 'DOP_SER_2'
DOP_SER_3_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=14)
Idents(list_subcluster_celltypes_MA[[5]], cells = DOP_SER_3_Int) <- 'DOP_SER_3'

TY_1_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=4)
Idents(list_subcluster_celltypes_MA[[5]], cells = TY_1_Int) <- 'TY_1'
TY_2_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=5)
Idents(list_subcluster_celltypes_MA[[5]], cells = TY_2_Int) <- 'TY_2'

OCTY_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=8)
Idents(list_subcluster_celltypes_MA[[5]], cells = OCTY_Int) <- 'OCTY'

Tbh_1_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=12)
Idents(list_subcluster_celltypes_MA[[5]], cells = Tbh_1_Int) <- 'Tbh_1'

Tbh_2_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=13)
Idents(list_subcluster_celltypes_MA[[5]], cells = Tbh_2_Int) <- 'Tbh_2'

SER_1_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=10)
Idents(list_subcluster_celltypes_MA[[5]], cells = SER_1_Int) <- 'SER_1'

SER_2_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=16)
Idents(list_subcluster_celltypes_MA[[5]], cells = SER_2_Int) <- 'SER_2'

Doublet_1_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=3)
Idents(list_subcluster_celltypes_MA[[5]], cells = Doublet_1_Int) <- 'Doublet_1'
Doublet_2_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=6)
Idents(list_subcluster_celltypes_MA[[5]], cells = Doublet_2_Int) <- 'Doublet_2'
Doublet_3_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=9)
Idents(list_subcluster_celltypes_MA[[5]], cells = Doublet_3_Int) <- 'Doublet_3'
Doublet_4_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=17)
Idents(list_subcluster_celltypes_MA[[5]], cells = Doublet_4_Int) <- 'Doublet_4'

DimPlot(list_subcluster_celltypes_MA[[5]], label=T) + NoLegend() + NoAxes()

dimplot_list_ma_annotated <- plot_sp_dimplot_merge(list_subcluster_celltypes_MA)
plot_grid_merge(dimplot_list_ma_annotated)
ggsave("Plots/celltype_groups/MA_sp_subtypes_annotated.png", width=10, height=11)

list_subcluster_celltypes_MA <- lapply(X=list_subcluster_celltypes_MA, FUN=add_cluster_label)

MA_cell_names_integrated <- rownames(list_subcluster_celltypes_MA[[5]]@meta.data)
MA_cell_names_Dmel <- rownames(list_subcluster_celltypes_MA[[1]]@meta.data)
MA_cell_names_Dsim <- rownames(list_subcluster_celltypes_MA[[2]]@meta.data)
MA_cell_names_Dsec <- rownames(list_subcluster_celltypes_MA[[3]]@meta.data)

MA_cell_clusters_integrted <- list_subcluster_celltypes_MA[[5]]@meta.data$ClusterLabel
MA_cell_clusters_Dmel <- list_subcluster_celltypes_MA[[1]]@meta.data$ClusterLabel
MA_cell_clusters_Dsim <- list_subcluster_celltypes_MA[[2]]@meta.data$ClusterLabel
MA_cell_clusters_Dsec <- list_subcluster_celltypes_MA[[3]]@meta.data$ClusterLabel

MA_cell_origin <- list_subcluster_celltypes_MA[[5]]@meta.data$orig.ident
MA_cell_type <- list_subcluster_celltypes_MA[[5]]@meta.data$orig.ident

MA_cell_info_integrated <- data.frame(Cell = MA_cell_names_integrated, Cluster = MA_cell_clusters_integrted, orig= MA_cell_origin) %>%
  dplyr::mutate(dataset='integrated')
MA_cell_info_Dmel <- data.frame(Cell = MA_cell_names_Dmel, Cluster = MA_cell_clusters_Dmel) %>%
  dplyr::mutate(dataset='Dmel')
MA_cell_info_Dsim <- data.frame(Cell = MA_cell_names_Dsim, Cluster = MA_cell_clusters_Dsim) %>%
  dplyr::mutate(dataset='Dsim')
MA_cell_info_Dsec <- data.frame(Cell = MA_cell_names_Dsec, Cluster = MA_cell_clusters_Dsec) %>%
  dplyr::mutate(dataset='Dsec')

MA_cell_info_integrated$Cell <- gsub("_[1-4]", "", MA_cell_info_integrated$Cell)
MA_cell_info_integrated$orig <- gsub("_rep[1-6]", "", MA_cell_info_integrated$orig)

# Apply the function to each subset
MA_cell_info_join_Dmel <- process_subset_cell_info(MA_cell_info_integrated, 'Dmel', MA_cell_info_Dmel)
MA_cell_info_join_Dsim <- process_subset_cell_info(MA_cell_info_integrated, 'Dsim', MA_cell_info_Dsim)
MA_cell_info_join_Dsec <- process_subset_cell_info(MA_cell_info_integrated, 'Dsec', MA_cell_info_Dsec)

# Combine the results
MA_cell_info_join_Trio <- dplyr::bind_rows(MA_cell_info_join_Dmel, MA_cell_info_join_Dsim, MA_cell_info_join_Dsec)

generate_sankey_celltype(df=MA_cell_info_join_Trio, target_ct=c(MA_celltypes, 'Doublet', 'Tbh'), target_sp=c("Dmel","Dsec"), filter_threshold =0)

for (i in 1:3) {
  
  pair <- pair_list[[i]]
  
  sp1 <- pair[1]
  sp2 <- pair[2]
  
  generate_sankey_celltype(df=MA_cell_info_join_Trio, target_ct=c(MA_celltypes, 'Doublet', 'Tbh'), target_sp=pair, filter_threshold =0)
  
  ggsave(glue::glue("Plots/sankey/pl_MA_reannotaed_{sp1}{sp2}.png"), width=12, height=12)
  
}

### Poxn cell types ###

Poxn_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                   individual=Trio_ISub_DF_list_labeled, celltype='Poxn',
                                   color='red')
plot_grid_merge(Poxn_plots)
ggsave("Plots/celltype_groups/Poxn_plots.png", width=10, height=11)

list_subcluster_celltypes_Poxn <- Get_subcluster_celltypes(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                         individual=Trio_ISub_DF_list_labeled, celltype='Poxn', res = 1) ## default = 0.1
DimPlot(list_subcluster_celltypes_Poxn[[4]], label=T)

dimplot_list_Poxn <- plot_sp_dimplot_merge(list_subcluster_celltypes_Poxn)
plot_grid_merge(dimplot_list_Poxn)
ggsave("Plots/celltype_groups/Poxn_sp_subtypes.png", width=10, height=11)

#Poxn_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Poxn[[4]])
Poxn_Dmel_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Poxn[[1]])
Poxn_Dsim_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Poxn[[2]])

Poxn_markers <- c("Poxn", "Nos", "caup", "Imp", "dpr8", "klg", "trol", "DIP-beta", "zld", "ara", "CG14687")

dotplot_list_Poxn <- plot_sp_dotplot_merge(list_subcluster_celltypes_Poxn, Poxn_markers)
plot_grid_merge(dotplot_list_Poxn)
ggsave("Plots/celltype_groups/Poxn_sp_subtypes_dotplots.png", width=15, height=10)

DotPlot(list_subcluster_celltypes_Poxn[[4]], features=Poxn_markers)

## Re-annotate Poxn clusters ##

## Dmel - 0: Poxn_1, 1: Poxn_2, 2: Poxn_3, 3: Poxn_2, 4: Poxn_4, 5: Poxn_5 

Poxn_1_Dmel <- WhichCells(list_subcluster_celltypes_Poxn[[1]], idents=0)
Idents(list_subcluster_celltypes_Poxn[[1]], cells = Poxn_1_Dmel) <- 'Poxn_1'
Poxn_2_Dmel <- WhichCells(list_subcluster_celltypes_Poxn[[1]], idents=c(1,3))
Idents(list_subcluster_celltypes_Poxn[[1]], cells = Poxn_2_Dmel) <- 'Poxn_2'
Poxn_3_Dmel <- WhichCells(list_subcluster_celltypes_Poxn[[1]], idents=2)
Idents(list_subcluster_celltypes_Poxn[[1]], cells = Poxn_3_Dmel) <- 'Poxn_3'
Poxn_4_Dmel <- WhichCells(list_subcluster_celltypes_Poxn[[1]], idents=4)
Idents(list_subcluster_celltypes_Poxn[[1]], cells = Poxn_4_Dmel) <- 'Poxn_4'
Poxn_5_Dmel <- WhichCells(list_subcluster_celltypes_Poxn[[1]], idents=5)
Idents(list_subcluster_celltypes_Poxn[[1]], cells = Poxn_5_Dmel) <- 'Poxn_5'
DimPlot(list_subcluster_celltypes_Poxn[[1]], label=T) + NoLegend() + NoAxes()

## Dsim - Poxn_1: 0, Poxn_2: 2, Poxn_3: 1, Poxn_4 : 3, Poxn_5 : 4

Poxn_1_Dsim <- WhichCells(list_subcluster_celltypes_Poxn[[2]], idents=0)
Idents(list_subcluster_celltypes_Poxn[[2]], cells = Poxn_1_Dsim) <- 'Poxn_1'
Poxn_2_Dsim <- WhichCells(list_subcluster_celltypes_Poxn[[2]], idents=2)
Idents(list_subcluster_celltypes_Poxn[[2]], cells = Poxn_2_Dsim) <- 'Poxn_2'
Poxn_3_Dsim <- WhichCells(list_subcluster_celltypes_Poxn[[2]], idents=1)
Idents(list_subcluster_celltypes_Poxn[[2]], cells = Poxn_3_Dsim) <- 'Poxn_3'
Poxn_4_Dsim <- WhichCells(list_subcluster_celltypes_Poxn[[2]], idents=3)
Idents(list_subcluster_celltypes_Poxn[[2]], cells = Poxn_4_Dsim) <- 'Poxn_4'
Poxn_5_Dsim <- WhichCells(list_subcluster_celltypes_Poxn[[2]], idents=4)
Idents(list_subcluster_celltypes_Poxn[[2]], cells = Poxn_5_Dsim) <- 'Poxn_5'
DimPlot(list_subcluster_celltypes_Poxn[[2]], label=T) + NoLegend() + NoAxes()

## Dsec - Poxn_1: 1, Poxn_2: 2 & 3, Poxn_3: 0, Poxn_4 : 4, Poxn_5 : 5

Poxn_1_Dsec <- WhichCells(list_subcluster_celltypes_Poxn[[3]], idents=1)
Idents(list_subcluster_celltypes_Poxn[[3]], cells = Poxn_1_Dsec) <- 'Poxn_1'
Poxn_2_Dsec <- WhichCells(list_subcluster_celltypes_Poxn[[3]], idents=c(2,3))
Idents(list_subcluster_celltypes_Poxn[[3]], cells = Poxn_2_Dsec) <- 'Poxn_2'
Poxn_3_Dsec <- WhichCells(list_subcluster_celltypes_Poxn[[3]], idents=0)
Idents(list_subcluster_celltypes_Poxn[[3]], cells = Poxn_3_Dsec) <- 'Poxn_3'
Poxn_4_Dsec <- WhichCells(list_subcluster_celltypes_Poxn[[3]], idents=4)
Idents(list_subcluster_celltypes_Poxn[[3]], cells = Poxn_4_Dsec) <- 'Poxn_4'
Poxn_5_Dsec <- WhichCells(list_subcluster_celltypes_Poxn[[3]], idents=5)
Idents(list_subcluster_celltypes_Poxn[[3]], cells = Poxn_5_Dsec) <- 'Poxn_5'
DimPlot(list_subcluster_celltypes_Poxn[[3]], label=T) + NoLegend() + NoAxes()

## DsecNoni - Poxn_1: 0, Poxn_2 : 2,4, Poxn_3 : 1, Poxn_4 : 3, Poxn_5 : 5

Poxn_1_DsecNoni <- WhichCells(list_subcluster_celltypes_Poxn[[4]], idents=0)
Idents(list_subcluster_celltypes_Poxn[[4]], cells = Poxn_1_DsecNoni) <- 'Poxn_1'
Poxn_2_DsecNoni <- WhichCells(list_subcluster_celltypes_Poxn[[4]], idents=c(2,4))
Idents(list_subcluster_celltypes_Poxn[[4]], cells = Poxn_2_DsecNoni) <- 'Poxn_2'
Poxn_3_DsecNoni <- WhichCells(list_subcluster_celltypes_Poxn[[4]], idents=1)
Idents(list_subcluster_celltypes_Poxn[[4]], cells = Poxn_3_DsecNoni) <- 'Poxn_3'
Poxn_4_DsecNoni <- WhichCells(list_subcluster_celltypes_Poxn[[4]], idents=3)
Idents(list_subcluster_celltypes_Poxn[[4]], cells = Poxn_4_DsecNoni) <- 'Poxn_4'
Poxn_5_DsecNoni <- WhichCells(list_subcluster_celltypes_Poxn[[4]], idents=5)
Idents(list_subcluster_celltypes_Poxn[[4]], cells = Poxn_5_DsecNoni) <- 'Poxn_5'
DimPlot(list_subcluster_celltypes_Poxn[[4]], label=T) + NoLegend() + NoAxes()

## Integrted - Poxn_1: 1&2, Poxn_2: 3&4&6, Poxn_3: 0, Poxn_4 : 7, Poxn_5 : 5

Poxn_1_Int <- WhichCells(list_subcluster_celltypes_Poxn[[5]], idents=c(1,2))
Idents(list_subcluster_celltypes_Poxn[[5]], cells = Poxn_1_Int) <- 'Poxn_1'
Poxn_2_Int <- WhichCells(list_subcluster_celltypes_Poxn[[5]], idents=c(3,4,6))
Idents(list_subcluster_celltypes_Poxn[[5]], cells = Poxn_2_Int) <- 'Poxn_2'
Poxn_3_Int <- WhichCells(list_subcluster_celltypes_Poxn[[5]], idents=c(0,7))
Idents(list_subcluster_celltypes_Poxn[[5]], cells = Poxn_3_Int) <- 'Poxn_3'
Poxn_4_Int <- c(paste0(Poxn_4_Dmel,"_1"), paste0(Poxn_4_Dsim,"_2"), paste0(Poxn_4_Dsec,"_3"),paste0(Poxn_4_DsecNoni,"_4"))

Poxn_5_Int <- WhichCells(list_subcluster_celltypes_Poxn[[5]], idents=5)
Idents(list_subcluster_celltypes_Poxn[[5]], cells = Poxn_5_Int) <- 'Poxn_5'
DimPlot(list_subcluster_celltypes_Poxn[[5]], label=T) + NoLegend() + NoAxes()

dimplot_list_Poxn_annotated <- plot_sp_dimplot_merge(list_subcluster_celltypes_Poxn)
plot_grid_merge(dimplot_list_Poxn_annotated)
ggsave("Plots/celltype_groups/Poxn_sp_subtypes_annotated.png", width=10, height=11)

DimPlot(TrioBrain.integrated_slim_ISub_DF_labeled, label=F,shuffle=T,
        cells.highlight= list(highlight=c(paste0(Poxn_4_Dsec,"_3"),paste0(Poxn_4_Dsim,"_2"))), 
        reduction = "tsne", 
        pt.size = 0.01,sizes.highlight = 1,
        cols.highlight = 'red', cols= "grey40", raster=FALSE) + NoLegend() + NoAxes()

## 


## Generate saneky plot ##

list_subcluster_celltypes_Poxn <- lapply(X=list_subcluster_celltypes_Poxn, FUN=add_cluster_label)

Poxn_cell_names_integrated <- rownames(list_subcluster_celltypes_Poxn[[5]]@meta.data)
Poxn_cell_names_Dmel <- rownames(list_subcluster_celltypes_Poxn[[1]]@meta.data)
Poxn_cell_names_Dsim <- rownames(list_subcluster_celltypes_Poxn[[2]]@meta.data)
Poxn_cell_names_Dsec <- rownames(list_subcluster_celltypes_Poxn[[3]]@meta.data)

Poxn_cell_clusters_integrted <- list_subcluster_celltypes_Poxn[[5]]@meta.data$ClusterLabel
Poxn_cell_clusters_Dmel <- list_subcluster_celltypes_Poxn[[1]]@meta.data$ClusterLabel
Poxn_cell_clusters_Dsim <- list_subcluster_celltypes_Poxn[[2]]@meta.data$ClusterLabel
Poxn_cell_clusters_Dsec <- list_subcluster_celltypes_Poxn[[3]]@meta.data$ClusterLabel

Poxn_cell_origin <- list_subcluster_celltypes_Poxn[[5]]@meta.data$orig.ident
Poxn_cell_type <- list_subcluster_celltypes_Poxn[[5]]@meta.data$orig.ident

Poxn_cell_info_integrated <- data.frame(Cell = Poxn_cell_names_integrated, Cluster = Poxn_cell_clusters_integrted, orig= Poxn_cell_origin) %>%
  dplyr::mutate(dataset='integrated')
Poxn_cell_info_Dmel <- data.frame(Cell = Poxn_cell_names_Dmel, Cluster = Poxn_cell_clusters_Dmel) %>%
  dplyr::mutate(dataset='Dmel')
Poxn_cell_info_Dsim <- data.frame(Cell = Poxn_cell_names_Dsim, Cluster = Poxn_cell_clusters_Dsim) %>%
  dplyr::mutate(dataset='Dsim')
Poxn_cell_info_Dsec <- data.frame(Cell = Poxn_cell_names_Dsec, Cluster = Poxn_cell_clusters_Dsec) %>%
  dplyr::mutate(dataset='Dsec')

Poxn_cell_info_integrated$Cell <- gsub("_[1-4]", "", Poxn_cell_info_integrated$Cell)
Poxn_cell_info_integrated$orig <- gsub("_rep[1-6]", "", Poxn_cell_info_integrated$orig)

# Apply the function to each subset
Poxn_cell_info_join_Dmel <- process_subset_cell_info(Poxn_cell_info_integrated, 'Dmel', Poxn_cell_info_Dmel)
Poxn_cell_info_join_Dsim <- process_subset_cell_info(Poxn_cell_info_integrated, 'Dsim', Poxn_cell_info_Dsim)
Poxn_cell_info_join_Dsec <- process_subset_cell_info(Poxn_cell_info_integrated, 'Dsec', Poxn_cell_info_Dsec)

# Combine the results
Poxn_cell_info_join_Trio <- dplyr::bind_rows(Poxn_cell_info_join_Dmel, Poxn_cell_info_join_Dsim, Poxn_cell_info_join_Dsec)

for (i in 1:3) {
  
  pair <- pair_list[[i]]
  
  sp1 <- pair[1]
  sp2 <- pair[2]
  
  generate_sankey_celltype(df=Poxn_cell_info_join_Trio, target_ct='Poxn', target_sp=pair, filter_threshold =0)
  
  ggsave(glue::glue("Plots/sankey/pl_Poxn_reannotaed_{sp1}{sp2}.png"), width=12, height=12)
  
}


### Clock cell types ###

Clock_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                     individual=Trio_ISub_DF_list_labeled, celltype='Clock',
                                     color='red')
plot_grid_merge(Clock_plots, int=4)
ggsave("Plots/celltype_groups/Clock_groups.png", width=10, height=11)

list_subcluster_celltypes_Clock <- Get_subcluster_celltypes(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                           individual=Trio_ISub_DF_list_labeled, celltype='Clock', res=0.5)

for (i in 1:5) {
  
  DefaultAssay(list_subcluster_celltypes_Clock[[i]]) <- "SCT"
  list_subcluster_celltypes_Clock[[i]] <- PrepSCTFindMarkers(list_subcluster_celltypes_Clock[[i]])
  
}

Clock_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Clock[[5]], only.pos = T)
Clock_Dmel_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Clock[[1]], only.pos = T)

dimplot_list_Clock <- plot_sp_dimplot_merge(list_subcluster_celltypes_Clock)
plot_grid_merge(dimplot_list_Clock)
ggsave("Plots/celltype_groups/Clock_sp_subtypes.png", width=10, height=11)

Clock_markers <- c("Clk","tim", "pros","Imp", "AstC", "Dh44")

FeaturePlot(list_subcluster_celltypes_Clock[[1]], feature=c("pros", "Imp", "Dh44", "AstC", "tim", "Clk"), ncol=3, order=T)
FeaturePlot(list_subcluster_celltypes_Clock[[2]], feature=c("pros", "Imp", "Dh44", "AstC"), ncol=2)
FeaturePlot(list_subcluster_celltypes_Clock[[3]], feature=c("pros", "Imp", "Dh44", "AstC"), ncol=2)
FeaturePlot(list_subcluster_celltypes_Clock[[5]], feature=c("pros", "Imp", "Dh44", "AstC", "tim"), order=T)

dotplot_list_Clock <- plot_sp_dotplot_merge(list_subcluster_celltypes_Clock, Clock_markers)
plot_grid_merge(dotplot_list_Clock)
ggsave("Plots/celltype_groups/Clock_sp_subtypes_dotplots.png", width=15, height=8)

## Re-annotate Clock clusters ##

## Dmel, Dsim, Dsec - Clock_pros : 0, Clock_Imp : 1

Clock_pros_Dmel <- WhichCells(list_subcluster_celltypes_Clock[[1]], idents=0)
Idents(list_subcluster_celltypes_Clock[[1]], cells = Clock_pros_Dmel) <- 'Clock_pros'
Clock_Imp_Dmel <- WhichCells(list_subcluster_celltypes_Clock[[1]], idents=1)
Idents(list_subcluster_celltypes_Clock[[1]], cells = Clock_Imp_Dmel) <- 'Clock_Imp'

Clock_pros_Dsim <- WhichCells(list_subcluster_celltypes_Clock[[2]], idents=0)
Idents(list_subcluster_celltypes_Clock[[2]], cells = Clock_pros_Dsim) <- 'Clock_pros'
Clock_Imp_Dsim <- WhichCells(list_subcluster_celltypes_Clock[[2]], idents=1)
Idents(list_subcluster_celltypes_Clock[[2]], cells = Clock_Imp_Dsim) <- 'Clock_Imp'

Clock_pros_Dsec <- WhichCells(list_subcluster_celltypes_Clock[[3]], idents=0)
Idents(list_subcluster_celltypes_Clock[[3]], cells = Clock_pros_Dsec) <- 'Clock_pros'
Clock_Imp_Dsec <- WhichCells(list_subcluster_celltypes_Clock[[3]], idents=1)
Idents(list_subcluster_celltypes_Clock[[3]], cells = Clock_Imp_Dsec) <- 'Clock_Imp'

## Int - Clock_pros_1 : 0,1, Clock_pros_2 : 3,6, Clock_pros_3 : 2, Clock_Imp_1 : 5,8, Clock_Imp_2 : 4,7,9,10

Clock_pros_1_Int <- WhichCells(list_subcluster_celltypes_Clock[[5]], idents=c(0,1))
Idents(list_subcluster_celltypes_Clock[[5]], cells = Clock_pros_1_Int) <- 'Clock_pros_1'
Clock_pros_2_Int <- WhichCells(list_subcluster_celltypes_Clock[[5]], idents=c(3,6))
Idents(list_subcluster_celltypes_Clock[[5]], cells = Clock_pros_2_Int) <- 'Clock_pros_2'
Clock_pros_3_Int <- WhichCells(list_subcluster_celltypes_Clock[[5]], idents=2)
Idents(list_subcluster_celltypes_Clock[[5]], cells = Clock_pros_3_Int) <- 'Clock_pros_3'
Clock_Imp_1_Int <- WhichCells(list_subcluster_celltypes_Clock[[5]], idents=c(5,8))
Idents(list_subcluster_celltypes_Clock[[5]], cells = Clock_Imp_1_Int) <- 'Clock_Imp_1'
Clock_Imp_2_Int <- WhichCells(list_subcluster_celltypes_Clock[[5]], idents=c(4,7,9,10))
Idents(list_subcluster_celltypes_Clock[[5]], cells = Clock_Imp_2_Int) <- 'Clock_Imp_2'

### Final Plot ###

dimplot_list_Clock_annotated <- plot_sp_dimplot_merge(list_subcluster_celltypes_Clock)
plot_grid_merge(dimplot_list_Clock_annotated)
ggsave("Plots/celltype_groups/Clock_sp_subtypes_annotated.png", width=10, height=11)

dotplot_list_Clock <- plot_sp_dotplot_merge(list_subcluster_celltypes_Clock, Clock_markers)
plot_grid_merge(dotplot_list_Clock)
ggsave("Plots/celltype_groups/Clock_sp_subtypes_annotated_dotplots_final.png", width=15, height=10)

Clock_cell_sankey_final <- sankey_plots_trio(list_celltypes=list_subcluster_celltypes_Clock, celltypes = 'Clock')

cowplot::plot_grid(plotlist=Clock_cell_sankey_final, nrow=1)

ggsave("Plots/celltype_groups/Clock_sp_subtypes_annotated_sankey_final.png", width=20, height=8)


### Fru+ cell types ###

Fru_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                   individual=Trio_ISub_DF_list_labeled, celltype='Fru',
                                   color='red')
plot_grid_merge(Fru_plots)
ggsave("Plots/celltype_groups/Fru_groups.png", width=10, height=11)

list_subcluster_celltypes_Fru_from_Int <- Get_subcluster_celltypes_from_Int(integrated=TrioBrain.integrated_slim_ISub_DF_labeled,
                                              celltype='Fru', res=1)

dimplot_list_Fru <- plot_sp_dimplot_merge(list_subcluster_celltypes_Fru_from_Int)
plot_grid_merge(dimplot_list_Fru)
ggsave("Plots/celltype_groups/Fru_sp_subtypes_from_Int.png", width=10, height=11)

for (i in 1:5) {
  
  DefaultAssay(list_subcluster_celltypes_Fru_from_Int[[i]]) <- "SCT"
  list_subcluster_celltypes_Fru_from_Int[[i]] <- PrepSCTFindMarkers(list_subcluster_celltypes_Fru_from_Int[[i]])
  
}

Fru_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Fru_from_Int[[5]], only.pos = T)
Fru_Dmel_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Fru_from_Int[[1]], only.pos = T)

Fru_Int_subtype_markers_816 <- Fru_Int_subtype_markers %>%
  dplyr::filter(cluster %in% c(8,16)) %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(n()==2)

Fru_markers <- c("pros", "Imp","Ms", "foxo", "VGlut", "dve", "br", "Wnt4", "CG43689")

FeaturePlot(list_subcluster_celltypes_Fru_from_Int[[1]], feature=Fru_markers, ncol=3, order=T)
FeaturePlot(list_subcluster_celltypes_Fru_from_Int[[5]], feature=Fru_markers, ncol=3, order=T)

dotplot_list_Fru <- plot_sp_dotplot_merge(list_subcluster_celltypes_Fru_from_Int, Fru_markers)
plot_grid_merge(dotplot_list_Fru)
ggsave("Plots/celltype_groups/Fru_sp_subtypes_dotplots.png", width=15, height=10)

### Manual insepction of NP cells ###

## neuropeptide ##

Dmel_gene_name <- read.table(file="genelist/Dmel_ID_genename.tsv") %>% dplyr::distinct(V1,V2) %>%
  dplyr::rename(FBID=V1, gene=V2)

df_neuropep <- data.table::fread("genelist/list_neuropeptide_FBgg0000179.txt", header = F) %>%
  dplyr::mutate(class="Neuropeptide") %>%
  dplyr::select(FBID=V1, class) %>%
  dplyr::left_join(., Dmel_gene_name, by="FBID") %>% ### few genes without matching gene_name
  na.omit()

neuropep_genes <- df_neuropep$gene

np_celltype_list <- NULL

for (i in 1:length(neuropep_genes)) {
  
  clusterx <- unique(cell_info_join_Trio$Cluster.x)
  clustery <- unique(cell_info_join_Trio$Cluster.y)
  
  np_celltype_list <- c(np_celltype_list, as.character(clusterx[grepl(neuropep_genes[i],clusterx)]), as.character(clustery[grepl(neuropep_genes[i],clustery)]))
  
}

NP_celltypes <- unique(np_celltype_list)

## Dividing them by specific NP groups ##

## Proc ##

Proc_celltypes <- NP_celltypes[grepl('Proc', NP_celltypes)]

Proc_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                    individual=Trio_ISub_DF_list_labeled, celltype=Proc_celltypes,
                                    color='red')
plot_grid_merge(Proc_plots, int=4)
ggsave("Plots/celltype_groups/Proc_groups.png", width=10, height=11)

list_subcluster_celltypes_Proc <- Get_subcluster_celltypes(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                           individual=Trio_ISub_DF_list_labeled, celltype=Proc_celltypes, res=0.1)

Proc_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Proc[[5]])

dimplot_list_Proc <- plot_sp_dimplot_merge(list_subcluster_celltypes_Proc)
plot_grid_merge(dimplot_list_Proc)
ggsave("Plots/celltype_groups/Proc_sp_subtypes.png", width=10, height=11)

Proc_markers <- c("crb","CG32204","Pka-C1","mamo","sNPF", "rn")

dotplot_list_Proc <- plot_sp_dotplot_merge(list_subcluster_celltypes_Proc, Proc_markers)
plot_grid_merge(dotplot_list_Proc)
ggsave("Plots/celltype_groups/Proc_sp_subtypes_dotplots.png", width=15, height=8)




