library(Seurat)
library(ggplot2)
library(dplyr)
library(ggsankey) #devtools::install_github("davidsjoberg/ggsankey")
library(readr)
library(ggrepel)
library(sctransform)
library(cowplot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/shared_genes.RData")
#load("Processed_Data/cell_info_join_Trio.RData")

#TrioBrain.integrated_slim_ISub_DF_labeled <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_ISub_DF_labeled.rds")
TrioBrain.integrated_slim_ISub_DF_labeled_species <- TrioBrain.integrated_slim_ISub_DF_labeled
TrioBrain.integrated_slim_ISub_DF_labeled_species$orig.ident <- gsub("_rep[1-6]", "", TrioBrain.integrated_slim_ISub_DF_labeled_species$orig.ident)

species_list <- c("Dmel", "Dsim", "Dsec", "DsecNoni")
pair_list <- list(c("Dmel","Dsec"),c("Dmel","Dsim"), c("Dsim","Dsec"), c("Dsec", "DsecNoni"))

Trio_ISub_DF_list <- readRDS(file = "Processed_Data/Trio_ISub_DF_list.rds")

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
    
    integrated_subset <- subset(integrated, orig.ident %in% sp_rep)
    
    filter_target_cluster_integrated_sp <- intersect(unique(integrated_subset$ClusterLabel), filter_target_cluster_integrated)
    
    subcluster_list[[i]] <- Get_subcluster_sp(integrated_subset, cluster=filter_target_cluster_integrated_sp, resolution = res)
    
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
plot_sp_dimplot_merge <- function(plot_list, pts=1) {
  
  dimplot_list <- list()
  
  for (i in 1:3) {
    
    dimplot_list[[i]] <- DimPlot(plot_list[[i]], label=T, pt.size=pts) + NoAxes() + ggtitle(species_list[[i]]) + NoLegend()
    
  }
  
  dimplot_list[[5]] <- DimPlot(plot_list[[5]], label=T, pt.size=pts/2) + NoAxes() + ggtitle('Integrated') + NoLegend()
  
  return(dimplot_list)
  
}
plot_sp_featureplot_merge <- function(plot_list, markers=marker_gene, n_col=ncol, ptsize=0.5) {
  
  featureplot_list <- list()
  
  for (i in 1:4) {
    
    DefaultAssay(plot_list[[i]]) <- "SCT"
    
    featureplot_list[[i]] <- FeaturePlot(plot_list[[i]], features=markers, ncol=n_col, order=T, pt.size=ptsize) * NoAxes()  * 
      theme(legend.position = 'bottom', legend.key.width = unit(0.15, "in"), legend.key.height = unit(0.05, "in"), 
            legend.text = element_text(size=8.5, color='black', vjust=3),
            plot.title = element_text(hjust = 0.5, face='bold.italic'))  * 
      scale_color_gradient(breaks=c(0,1,2,3,4,5), low="lightgrey", high="blue")
    
  }
  
  DefaultAssay(plot_list[[5]]) <- "SCT"
  featureplot_list[[5]] <- FeaturePlot(plot_list[[5]], features=markers, ncol=n_col, order=T, pt.size=ptsize/2) * NoAxes() *
    theme(legend.position = 'bottom', legend.key.width = unit(0.15, "in"), legend.key.height = unit(0.05, "in"), 
          legend.text = element_text(size=8.5, color='black', vjust=3),
          plot.title = element_text(hjust = 0.5, face='bold.italic')) * 
    scale_color_gradient(breaks=c(0,1,2,3,4,5), low="lightgrey", high="blue")
  
  return(featureplot_list)
  
}
plot_sp_dotplot_merge <- function(plot_list, markers) {
  
  dotplot_list <- list()
  
  for (i in 1:3) {
    
    DefaultAssay(plot_list[[i]]) <- "SCT"
    
    dotplot_list[[i]] <- DotPlot(plot_list[[i]], features = markers) + ggtitle(species_list[[i]]) + 
      theme(axis.title = element_blank())
  }
  
  DefaultAssay(plot_list[[5]]) <- "SCT"
  
  dotplot_list[[5]] <- DotPlot(plot_list[[5]], features = markers) +  ggtitle('Integrated') + 
    theme(axis.title = element_blank())
  
  return(dotplot_list)
  
}
plot_cell_type_groups <- function(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                  individual=Trio_ISub_DF_list, celltype=glial_celltypes,
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
  
  plot_list[[5]] <- DimPlot(integrated, label=F,shuffle=T,
                            cells.highlight= list(highlight=highlight_cells), 
                            reduction = "tsne", 
                            pt.size = 0.01,sizes.highlight = 0.01,
                            cols.highlight = color, cols= "grey40", raster=FALSE) + NoLegend() + NoAxes() + ggtitle('Integrated')
  
  return(plot_list)
  
}
plot_grid_merge <- function(plot_list, int=5, n_col=2) {
  
  plot_out <- cowplot::plot_grid(plot_list[[int]], plot_list[[1]],plot_list[[2]],plot_list[[3]], ncol=n_col)
  
  return(plot_out)
  
}
plot_celltype_summary <- function(plotlist, marker_list, ncols=3, pt_size=0.01, labelsize=5, exclude = 'Unannotated', r_width = c(2,3)) {
  
  plotlist <- lapply(plotlist, FUN = add_cluster_label)
  
  if (!is.na(exclude)) {
    
    for (i in 1:5) {
      
      plotlist[[i]] <- subset(plotlist[[i]], idents = exclude, invert = T)
      
    }
    
  }
  
  metadata_labeled <- plotlist[[5]]@meta.data %>%
    tibble::rownames_to_column('Cell') %>%
    dplyr::select(Cell, orig.ident, ClusterLabel)

  featureplot_list <- plot_sp_featureplot_merge(plot_list=plotlist, markers=marker_list, n_col=ncols, ptsize=pt_size)
  
  plot_output_list <- list()
  
  for (i in 1:4) {
    
    Seurat_obj <- plotlist[[i]]
    
    Seurat_obj@meta.data <- Seurat_obj@meta.data %>%
      tibble::rownames_to_column('Cell') %>%
      dplyr::left_join(., metadata_labeled, by=c('Cell','orig.ident')) %>%
      tibble::column_to_rownames('Cell')
    
    Dimplot <- DimPlot(Seurat_obj, 
                                 group.by='ClusterLabel.y', label=T, label.size=labelsize, repel=T, pt.size=pt_size) + NoLegend() + NoAxes() + ggtitle(species_list[[i]]) +
      theme(plot.title = element_text(hjust = 0.5, face='bold.italic'))
    
    plot_output_list[[i]]  <- cowplot::plot_grid(Dimplot, featureplot_list[[i]], ncol=2, rel_widths = r_width)
    
  }

  Dimplot_Int <- DimPlot(plotlist[[5]], label=T, label.size=labelsize, repel=T, pt.size=pt_size/2) + NoLegend() + NoAxes() + ggtitle('Integrated')+
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_output_list[[5]]  <- cowplot::plot_grid(Dimplot_Int, featureplot_list[[5]], ncol=2, rel_widths = r_width)
  
  return(plot_output_list)
  
}

### Functions end ###

## celltype group plots ###

glial_celltypes <- c("ENS","AST","PRN","SUB", "CTX")
glial_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                     individual=Trio_ISub_DF_list, celltype=glial_celltypes,
                                     color='darkorange')
plot_grid_merge(glial_plots)
ggsave("Plots/celltype_groups/glial_groups.png", width=10, height=11)

SUB_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                     individual=Trio_ISub_DF_list_labeled, celltype='SUB',
                                     color='purple')
plot_grid_merge(SUB_plots)
ggsave("Plots/celltype_groups/SUBs.png", width=10, height=11)

 ## manual inspection ##

list_subcluster_celltypes_glia <- Get_subcluster_celltypes_from_Int(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                         celltype=glial_celltypes, res=0.1)

Glia_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_glia[[5]], only.pos = T)

dimplot_list_glia <- plot_sp_dimplot_merge(list_subcluster_celltypes_glia)
plot_grid_merge(dimplot_list_glia)
ggsave("Plots/celltype_groups/glia_sp_subtypes.png", width=10, height=11)

glial_markers <- c("axo", "trol", "Eaat1", "Tret1-1", "CG40470", "baz", "Hml")

dotplot_list_glia <- plot_sp_dotplot_merge(list_subcluster_celltypes_glia, glial_markers)
plot_grid_merge(dotplot_list_glia)
ggsave("Plots/celltype_groups/glia_sp_subtypes_dotplots.png", width=15, height=8)

## Re-annotate Glial clusters ##

## Integrted - ENS : 0,4, AST : 1,6, PRN : 2, CTX : 3, SUB : 5

ENS_Int <- WhichCells(list_subcluster_celltypes_glia[[5]], idents=c(0,4)) 
Idents(list_subcluster_celltypes_glia[[5]], cells = ENS_Int ) <- 'ENS'
AST_Int <- WhichCells(list_subcluster_celltypes_glia[[5]], idents=c(1,6))
Idents(list_subcluster_celltypes_glia[[5]], cells = AST_Int) <- 'AST'
PRN_Int <- WhichCells(list_subcluster_celltypes_glia[[5]], idents=2)
Idents(list_subcluster_celltypes_glia[[5]], cells = PRN_Int) <- 'PRN'
CTX_Int <- WhichCells(list_subcluster_celltypes_glia[[5]], idents=3)
Idents(list_subcluster_celltypes_glia[[5]], cells = CTX_Int) <- 'CTX'
SUB_Int <- WhichCells(list_subcluster_celltypes_glia[[5]], idents=5)
Idents(list_subcluster_celltypes_glia[[5]], cells = SUB_Int) <- 'SUB'
Unanno_glia_Int <- WhichCells(list_subcluster_celltypes_glia[[5]], idents=7)
Idents(list_subcluster_celltypes_glia[[5]], cells = Unanno_glia_Int) <- 'Unanno_glia'

Idents(list_subcluster_celltypes_glia[[1]], cells = Unanno_glia_Int) <- 'Unanno_glia'
Idents(list_subcluster_celltypes_glia[[2]], cells = Unanno_glia_Int) <- 'Unanno_glia'
Idents(list_subcluster_celltypes_glia[[3]], cells = Unanno_glia_Int) <- 'Unanno_glia'
Idents(list_subcluster_celltypes_glia[[4]], cells = Unanno_glia_Int) <- 'Unanno_glia'

glia_anno_plot_merge <- plot_celltype_summary(plotlist=list_subcluster_celltypes_glia, marker_list=glial_markers[1:6], 
                                            ncols=3, pt_size=0.01, labelsize=5,
                                            exclude = 'Unanno_glia') 
cowplot::plot_grid(glia_anno_plot_merge[[5]],glia_anno_plot_merge[[1]],glia_anno_plot_merge[[2]],glia_anno_plot_merge[[3]],
                   ncol=1)
ggsave("Plots/celltype_groups/glia_sp_subtypes_summary.png", width=20, height=30)

## merge all glia ##

Seurat_obj_trio <- TrioBrain.integrated_slim_ISub_DF_labeled

Idents(Seurat_obj_trio, cells = c(ENS_Int,AST_Int,PRN_Int,CTX_Int,SUB_Int)) <- 'glia'
Seurat_obj_trio[["ClusterLabel"]] <- Idents(object = Seurat_obj_trio)

Seurat_obj_trio_metadata <- Seurat_obj_trio@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))

Seurat_obj_trio_metadata_summary_rep <- Seurat_obj_trio_metadata %>%
  dplyr::group_by(species, ClusterLabel, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

Seurat_obj_trio_metadata_summary <- Seurat_obj_trio_metadata_summary_rep %>%
  dplyr::group_by(species, ClusterLabel) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

Seurat_obj_trio_metadata_summary %>%
  dplyr::filter(ClusterLabel=='glia') %>%
  ggplot(.) +
  geom_col(aes(x=species, y=percent_combined, fill=species), alpha = 1, position=position_dodge(width=0.8)) + 
  geom_errorbar(aes(x=species, y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,0,percent_combined-sem), ymax=percent_combined+sem, group=species), width=0.8, position=position_dodge(width=0.8)) +
  geom_point(data=dplyr::filter(Seurat_obj_trio_metadata_summary_rep, ClusterLabel=='glia'),
             size = 1, alpha =0.8,  
             aes(x=species, y=percent, group=species,fill=species), shape=21, position=position_dodge(width=0.8)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=ClusterLabel, y=n/total*100)) +
  theme_bw() +
  theme(axis.title.y = element_text(size=14, color='black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x=element_text(size=11, color='black', angle=45, hjust=1, vjust=1, face='italic'),
        legend.text = element_text(size=13, face='italic', color='black'),
        panel.grid = element_blank(),
        strip.text = element_text(size=12, color='black'),
        legend.position='none') +
  labs(x="Cluster", y="Frequency (% of cells)", fill="")  +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8","#984EA3")) +
  facet_wrap(~ClusterLabel, scales='free', ncol=8)

Seurat_obj_trio$orig.ident <- gsub("_rep[1-6]", "", Seurat_obj_trio$orig.ident)

DotPlot(subset(Seurat_obj_trio, idents='glia'), feature='repo', group.by='orig.ident')


## KC manual insepction ##

KC_celltypes <- c("αβ-KC","γ-KC","α'β'-KC")
KC_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                   individual=Trio_ISub_DF_list, celltype=KC_celltypes,
                                   color='red')
plot_grid_merge(KC_plots)
ggsave("Plots/celltype_groups/KC_groups.png", width=10, height=11)

list_subcluster_celltypes_KC <- Get_subcluster_celltypes_from_Int(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                                    celltype=KC_celltypes, res=0.1)

KC_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_KC[[5]], only.pos = T)

dimplot_list_KC <- plot_sp_dimplot_merge(list_subcluster_celltypes_KC)
plot_grid_merge(dimplot_list_KC)
ggsave("Plots/celltype_groups/KC_sp_subtypes.png", width=10, height=11)

KC_markers <- c("crb","CG32204","Pka-C1","mamo","sNPF", "rn")

dotplot_list_KC <- plot_sp_dotplot_merge(list_subcluster_celltypes_KC, KC_markers)
plot_grid_merge(dotplot_list_KC)
ggsave("Plots/celltype_groups/KC_sp_subtypes_dotplots.png", width=15, height=8)

## Re-annotate KC clusters ##

## Integrted - αβ-KC_1 : 0, αβ-KC_2 : 3, γ-KC : 1,  α'β'-KC : 2 

abKC_1_Int <- WhichCells(list_subcluster_celltypes_KC[[5]], idents=0)
Idents(list_subcluster_celltypes_KC[[5]], cells = abKC_1_Int) <- 'αβ-KC_1'
abKC_2_Int <- WhichCells(list_subcluster_celltypes_KC[[5]], idents=3)
Idents(list_subcluster_celltypes_KC[[5]], cells = abKC_2_Int) <- 'αβ-KC_2'
rKC_Int <- WhichCells(list_subcluster_celltypes_KC[[5]], idents=1)
Idents(list_subcluster_celltypes_KC[[5]], cells = rKC_Int) <- 'γ-KC'
apbpKC_Int <- WhichCells(list_subcluster_celltypes_KC[[5]], idents=2)
Idents(list_subcluster_celltypes_KC[[5]], cells = apbpKC_Int) <- "α'β'-KC"

KC_anno_plot_merge <- plot_celltype_summary(plotlist=list_subcluster_celltypes_KC, marker_list=c("crb", "mamo", "sNPF", "rn"), 
                                              ncols=2, pt_size=0.01, labelsize=5,
                                              exclude = NA) 
cowplot::plot_grid(KC_anno_plot_merge[[5]],KC_anno_plot_merge[[1]],KC_anno_plot_merge[[2]],KC_anno_plot_merge[[3]],
                   ncol=1)
ggsave("Plots/celltype_groups/KC_sp_subtypes_summary.png", width=20, height=30)


### Manual inspection of MA cells ###

MA_celltypes <- "MON"

MA_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                  individual=Trio_ISub_DF_list, celltype='MON',
                                  color='blue')
plot_grid_merge(MA_plots)
ggsave("Plots/celltype_groups/MA_groups.png", width=10, height=11)

list_subcluster_celltypes_MA <- Get_subcluster_celltypes_from_Int(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                                  celltype=MA_celltypes, res=0.1)

for (i in 1:5) {
  
  DefaultAssay(list_subcluster_celltypes_MA[[i]]) <- "SCT"
  list_subcluster_celltypes_MA[[i]] <- PrepSCTFindMarkers(list_subcluster_celltypes_MA[[i]])
  
}

MA_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_MA[[5]], only.pos = T)

FeaturePlot(list_subcluster_celltypes_MA[[3]], feature=c('pros', 'Imp'), order=T)

DotPlot(list_subcluster_celltypes_MA[[5]], feature='Vmat')

dimplot_list_MA <- plot_sp_dimplot_merge(list_subcluster_celltypes_MA)
plot_grid_merge(dimplot_list_MA)
ggsave("Plots/celltype_groups/MA_sp_subtypes.png", width=10, height=11)

MA_markers <- c("ple", "DAT","Trhn","SerT","Tdc2","Tbh")

DimPlot(list_subcluster_celltypes_MA[[5]], cells.highlight = WhichCells(list_subcluster_celltypes_MA[[5]], idents = 6))

dotplot_list_MA <- plot_sp_dotplot_merge(list_subcluster_celltypes_MA, MA_markers)
plot_grid_merge(dotplot_list_MA)
ggsave("Plots/celltype_groups/MA_sp_subtypes_dotplots.png", width=15, height=8)

plot_sp_featureplot_merge(list_subcluster_celltypes_MA, markers=MA_markers, n_col=3, ptsize=0.5) 

## Re-annotate MA clusters ##

## Integrted - DOP_1 : 0,2,3, DOP_2 : 6,7, DOP/SER : 1,11,14, SER : 10, TY_1 : 4, TY_2 : 5, TY_3 : 9, OCTY : 8, Tbh : 12

DOP_1_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=c(0,2,3))
Idents(list_subcluster_celltypes_MA[[5]], cells = DOP_1_Int) <- 'DOP_1'
DOP_2_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=c(6,7))
Idents(list_subcluster_celltypes_MA[[5]], cells = DOP_2_Int) <- 'DOP_2'
DOPSER_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=c(1,11,14))
Idents(list_subcluster_celltypes_MA[[5]], cells = DOPSER_Int) <- 'DOP/SER'
SER_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=10)
Idents(list_subcluster_celltypes_MA[[5]], cells = SER_Int) <- 'SER'
TY_1_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=4)
Idents(list_subcluster_celltypes_MA[[5]], cells = TY_1_Int) <- 'TY_1'
TY_2_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=5)
Idents(list_subcluster_celltypes_MA[[5]], cells = TY_2_Int) <- 'TY_2'
TY_3_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=9)
Idents(list_subcluster_celltypes_MA[[5]], cells = TY_3_Int) <- 'TY_3'
OCTY_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=8)
Idents(list_subcluster_celltypes_MA[[5]], cells = OCTY_Int) <- 'OCTY'
Tbh_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=12)
Idents(list_subcluster_celltypes_MA[[5]], cells = Tbh_Int) <- 'Tbh'
Unanno_MA_Int <- WhichCells(list_subcluster_celltypes_MA[[5]], idents=c(13,15,16,17))
Idents(list_subcluster_celltypes_MA[[5]], cells = Unanno_MA_Int) <- 'Unanno_MA'

Idents(list_subcluster_celltypes_MA[[1]], cells = Unanno_MA_Int) <- 'Unanno_MA'
Idents(list_subcluster_celltypes_MA[[2]], cells = Unanno_MA_Int) <- 'Unanno_MA'
Idents(list_subcluster_celltypes_MA[[3]], cells = Unanno_MA_Int) <- 'Unanno_MA'
Idents(list_subcluster_celltypes_MA[[4]], cells = Unanno_MA_Int) <- 'Unanno_MA'

MA_anno_plot_merge <- plot_celltype_summary(plotlist=list_subcluster_celltypes_MA, marker_list=c("pros", "Imp","DAT","SerT","Tdc2","Tbh"), 
                                            ncols=3, pt_size=0.01, labelsize=5,
                                            exclude = 'Unanno_MA') 
cowplot::plot_grid(MA_anno_plot_merge[[5]],MA_anno_plot_merge[[1]],MA_anno_plot_merge[[2]],MA_anno_plot_merge[[3]],
                   ncol=1)
ggsave("Plots/celltype_groups/MA_sp_subtypes_summary.png", width=20, height=30)

### Clock cell types ###

Clock_celltypes <- "Clock"

Clock_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                     individual=Trio_ISub_DF_list, celltype='Clock',
                                     color='red')
plot_grid_merge(Clock_plots)
ggsave("Plots/celltype_groups/Clock_groups.png", width=10, height=11)

list_subcluster_celltypes_Clock <- Get_subcluster_celltypes_from_Int(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                                  celltype=Clock_celltypes, res=0.1)

for (i in 1:5) {
  
  DefaultAssay(list_subcluster_celltypes_Clock[[i]]) <- "SCT"
  list_subcluster_celltypes_Clock[[i]] <- PrepSCTFindMarkers(list_subcluster_celltypes_Clock[[i]])
  
}

Clock_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Clock[[5]], only.pos = T)
FeaturePlot(list_subcluster_celltypes_Clock[[5]], feature=c('pros', 'Imp', 'Gad1', 'CNMa', 'CCHa1'), order=T)
FeaturePlot(list_subcluster_celltypes_Clock[[2]], feature=c('pros', 'Imp', 'Gad1', 'CNMa', 'CCHa1'), order=T)
FeaturePlot(list_subcluster_celltypes_Clock[[1]], feature=c('pros', 'Imp', 'Gad1', 'CNMa', 'CCHa1'), order=T)

dimplot_list_Clock <- plot_sp_dimplot_merge(list_subcluster_celltypes_Clock)
plot_grid_merge(dimplot_list_Clock)
ggsave("Plots/celltype_groups/Clock_sp_subtypes.png", width=10, height=11)

Clock_markers <- c('pros', 'Imp', 'Gad1')

dotplot_list_Clock <- plot_sp_dotplot_merge(list_subcluster_celltypes_Clock, Clock_markers)
plot_grid_merge(dotplot_list_Clock)
ggsave("Plots/celltype_groups/Clock_sp_subtypes_dotplots.png", width=12, height=8)

## Re-annotate Clock clusters ##

## Integrted - Clock_1 : 0,1,5 ; Clock_2 : 2,3,7,10,11,13 ; Clock_3 : 4,6,8,9,12

Clock_1_Int <- WhichCells(list_subcluster_celltypes_Clock[[5]], idents=c(0,1,5))
Idents(list_subcluster_celltypes_Clock[[5]], cells = Clock_1_Int) <- 'Clock_1'
Clock_2_Int <- WhichCells(list_subcluster_celltypes_Clock[[5]], idents=c(2,3,7,10,11,13))
Idents(list_subcluster_celltypes_Clock[[5]], cells = Clock_2_Int) <- 'Clock_2'
Clock_3_Int <- WhichCells(list_subcluster_celltypes_Clock[[5]], idents=c(4,6,8,9,12))
Idents(list_subcluster_celltypes_Clock[[5]], cells = Clock_3_Int) <- 'Clock_3'

Clock_anno_plot_merge <- plot_celltype_summary(plotlist=list_subcluster_celltypes_Clock, marker_list=Clock_markers, 
                                            ncols=3, pt_size=1, labelsize=5,
                                            exclude = NA, r_width=c(1,2)) 

cowplot::plot_grid(Clock_anno_plot_merge[[5]],Clock_anno_plot_merge[[1]],Clock_anno_plot_merge[[2]],Clock_anno_plot_merge[[3]],
                   ncol=1)
ggsave("Plots/celltype_groups/Clock_sp_subtypes_summary.png", width=12, height=12)


### OPN cell types ###

OPN_celltypes <- "OPN"

OPN_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                     individual=Trio_ISub_DF_list, celltype='OPN',
                                     color='red')
plot_grid_merge(OPN_plots)
ggsave("Plots/celltype_groups/OPN_groups.png", width=10, height=11)

list_subcluster_celltypes_OPN <- Get_subcluster_celltypes_from_Int(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                                     celltype=OPN_celltypes, res=0.1)

dimplot_list_OPN <- plot_sp_dimplot_merge(list_subcluster_celltypes_OPN)
plot_grid_merge(dimplot_list_OPN)
ggsave("Plots/celltype_groups/OPN_sp_subtypes.png", width=10, height=11)

OPN_markers <- c('pros', 'Imp', 'Gad1', 'ChAT')

dotplot_list_OPN <- plot_sp_dotplot_merge(list_subcluster_celltypes_OPN, OPN_markers)
plot_grid_merge(dotplot_list_OPN)
ggsave("Plots/celltype_groups/OPN_sp_subtypes_dotplots.png", width=12, height=8)

for (i in 1:5) {
  
  DefaultAssay(list_subcluster_celltypes_OPN[[i]]) <- "SCT"
  list_subcluster_celltypes_OPN[[i]] <- PrepSCTFindMarkers(list_subcluster_celltypes_OPN[[i]])
  
}

FeaturePlot(list_subcluster_celltypes_OPN[[5]], feature=c('pros', 'Imp'), order=T)
FeaturePlot(list_subcluster_celltypes_OPN[[3]], feature=c('pros', 'Imp'), order=T)
FeaturePlot(list_subcluster_celltypes_OPN[[1]], feature=c('pros', 'Imp'), order=T)

## Re-annotate OPN clusters ##

## Integrted - all into one cluster

OPN_Int <- WhichCells(list_subcluster_celltypes_OPN[[5]], idents=0:6)
Idents(list_subcluster_celltypes_OPN[[5]], cells = OPN_Int) <- 'OPN'

DimPlot(list_subcluster_celltypes_OPN[[5]])

OPN_anno_plot_merge <- plot_celltype_summary(plotlist=list_subcluster_celltypes_OPN, marker_list=c('pros', 'Imp'), 
                                               ncols=3, pt_size=0.5, labelsize=5,
                                               exclude = 'Unannotated_OPN', r_width=c(1,2)) 

cowplot::plot_grid(OPN_anno_plot_merge[[5]],OPN_anno_plot_merge[[1]],OPN_anno_plot_merge[[2]],OPN_anno_plot_merge[[3]],
                   ncol=1)
ggsave("Plots/celltype_groups/OPN_sp_subtypes_summary.png", width=12, height=12)


### Poxn cell types ###

Poxn_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                   individual=Trio_ISub_DF_list, celltype='Poxn',
                                   color='red')
plot_grid_merge(Poxn_plots)
ggsave("Plots/celltype_groups/Poxn_plots.png", width=10, height=11)

list_subcluster_celltypes_Poxn <- Get_subcluster_celltypes_from_Int(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                                   celltype='Poxn', res=0.1)

dimplot_list_Poxn <- plot_sp_dimplot_merge(list_subcluster_celltypes_Poxn)
plot_grid_merge(dimplot_list_Poxn)
ggsave("Plots/celltype_groups/Poxn_sp_subtypes.png", width=10, height=11)

for (i in 1:5) {
  
  DefaultAssay(list_subcluster_celltypes_Poxn[[i]]) <- "SCT"
  list_subcluster_celltypes_Poxn[[i]] <- PrepSCTFindMarkers(list_subcluster_celltypes_Poxn[[i]])
  
}

Poxn_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Poxn[[5]], only.pos = T)
Poxn_Dmel_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Poxn[[1]], only.pos = T)

FeaturePlot(list_subcluster_celltypes_Poxn[[5]], feature=c('klg', 'dpr8', 'trol', 'ara', 'Nos','Corin', 'vg'), order=T)

Poxn_Int_subtype_markers_256 <- Poxn_Int_subtype_markers %>%
  dplyr::filter(cluster %in% c(2,5,6)) %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(n()==3)

Poxn_markers <- c('Lmpt','AstC','Ptx1','DIP-alpha','CCHa2-R','trol')

FeaturePlot(list_subcluster_celltypes_Poxn[[5]], feature=Poxn_markers, order=T)

dotplot_list_Poxn <- plot_sp_dotplot_merge(list_subcluster_celltypes_Poxn, Poxn_markers)
plot_grid_merge(dotplot_list_Poxn)
ggsave("Plots/celltype_groups/Poxn_sp_subtypes_dotplots.png", width=20, height=8)

## Re-annotate Poxn clusters ##

## Poxn_1 : 2,5,6, Poxn_2 : 0, Poxn_3 : 1, Poxn_4 : 3, Poxn_5 : 7, Poxn_6 : 9, Poxn_unannotated : 4,8 

Poxn_1_Int <- WhichCells(list_subcluster_celltypes_Poxn[[5]], idents=c(2,5,6))
Idents(list_subcluster_celltypes_Poxn[[5]], cells = Poxn_1_Int) <- 'Poxn_1'
Poxn_2_Int <- WhichCells(list_subcluster_celltypes_Poxn[[5]], idents=0)
Idents(list_subcluster_celltypes_Poxn[[5]], cells = Poxn_2_Int) <- 'Poxn_2'
Poxn_3_Int <- WhichCells(list_subcluster_celltypes_Poxn[[5]], idents=1)
Idents(list_subcluster_celltypes_Poxn[[5]], cells = Poxn_3_Int) <- 'Poxn_3'
Poxn_4_Int <- WhichCells(list_subcluster_celltypes_Poxn[[5]], idents=3)
Idents(list_subcluster_celltypes_Poxn[[5]], cells = Poxn_4_Int) <- 'Poxn_4'
Poxn_5_Int <- WhichCells(list_subcluster_celltypes_Poxn[[5]], idents=7)
Idents(list_subcluster_celltypes_Poxn[[5]], cells = Poxn_5_Int) <- 'Poxn_5'
Poxn_6_Int <- WhichCells(list_subcluster_celltypes_Poxn[[5]], idents=9)
Idents(list_subcluster_celltypes_Poxn[[5]], cells = Poxn_6_Int) <- 'Poxn_6'
Unanno_Poxn_Int <- WhichCells(list_subcluster_celltypes_Poxn[[5]], idents=c(4,8))
Idents(list_subcluster_celltypes_Poxn[[5]], cells = Unanno_Poxn_Int) <- 'Unanno_Poxn'

Idents(list_subcluster_celltypes_Poxn[[1]], cells = Unanno_Poxn_Int) <- 'Unanno_Poxn'
Idents(list_subcluster_celltypes_Poxn[[2]], cells = Unanno_Poxn_Int) <- 'Unanno_Poxn'
Idents(list_subcluster_celltypes_Poxn[[3]], cells = Unanno_Poxn_Int) <- 'Unanno_Poxn'
Idents(list_subcluster_celltypes_Poxn[[4]], cells = Unanno_Poxn_Int) <- 'Unanno_Poxn'

Poxn_anno_plot_merge <- plot_celltype_summary(plotlist=list_subcluster_celltypes_Poxn, marker_list=Poxn_markers, 
                                               ncols=3, pt_size=1, labelsize=5,
                                               exclude = 'Unanno_Poxn', r_width=c(1,2)) 

cowplot::plot_grid(Poxn_anno_plot_merge[[5]],Poxn_anno_plot_merge[[1]],Poxn_anno_plot_merge[[2]],Poxn_anno_plot_merge[[3]],
                   ncol=1)
ggsave("Plots/celltype_groups/Poxn_sp_subtypes_summary.png", width=12, height=15)


### fru cell types ###

fru_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                    individual=Trio_ISub_DF_list, celltype='fru',
                                    color='red')
plot_grid_merge(fru_plots)
ggsave("Plots/celltype_groups/fru_plots.png", width=10, height=11)

list_subcluster_celltypes_fru <- Get_subcluster_celltypes_from_Int(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                                    celltype='fru', res=0.1)

dimplot_list_fru <- plot_sp_dimplot_merge(list_subcluster_celltypes_fru)
plot_grid_merge(dimplot_list_fru)
ggsave("Plots/celltype_groups/fru_sp_subtypes.png", width=10, height=11)

for (i in 1:5) {
  
  DefaultAssay(list_subcluster_celltypes_fru[[i]]) <- "SCT"
  list_subcluster_celltypes_fru[[i]] <- PrepSCTFindMarkers(list_subcluster_celltypes_fru[[i]])
  
}

fru_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_fru[[5]], only.pos = T)
fru_Dmel_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_fru[[1]], only.pos = T)

fru_Int_subtype_markers_15_22 <- FindMarkers(list_subcluster_celltypes_fru[[5]], ident.1 = 15, ident.2 = 22, only.pos = T)


FeaturePlot(list_subcluster_celltypes_fru[[5]], feature=c('fru', 'Ms', 'Dh44', 'Gad1', 'acj6', 'VGlut', 'Delta',
                                                          'Awh','spab', 'pnt','CG12239', 'Hug', 'Crz', 'Ilp2', 'ITP', 'Mip'), order=T)
FeaturePlot(list_subcluster_celltypes_fru[[5]], feature=c('twz'), order=T)

FeaturePlot(list_subcluster_celltypes_fru[[2]], feature=c('br', 'VGlut', 'Gad1'), order=T)

DimPlot(list_subcluster_celltypes_fru[[5]], cells.highlight = WhichCells(list_subcluster_celltypes_fru[[5]], idents = c(4,9,12)))

FeaturePlot(list_subcluster_celltypes_fru[[2]], feature=c('fru', 'Ms', 'Gad1', 'acj6', 'VGlut','Awh','spab', 'pnt','CG12239', 'Hug', 'Crz', 'Ilp2', 'ITP', 'Mip'), 
            order=T,pt.size=0.2)

fru_Int_subtype_markers_1522 <- fru_Int_subtype_markers %>%
  dplyr::filter(cluster %in% c(15,22)) %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(n()==2)

fru_Int_subtype_markers_7etal <- fru_Int_subtype_markers %>%
  dplyr::filter(cluster %in% c(7,2,6,14,23,25)) %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(n()==1) %>%
  dplyr::filter(cluster==7)

fru_markers <- c('fru', 'VGlut', 'Gad1', 'Awh', 'acj6', 'Delta', 'spab','Ms', 'pnt', 'CG12239', 'Mip', 'Hug', 'Crz', 'Ilp2', 'ITP')

dotplot_list_fru <- plot_sp_dotplot_merge(list_subcluster_celltypes_fru, fru_markers)
plot_grid_merge(dotplot_list_fru)
ggsave("Plots/celltype_groups/fru_sp_subtypes_dotplots.png", width=24, height=15)

## Re-annotate fru clusters ##

## Fru_1 : 0,1,19,21,26, Fru_2 : 2,6,14,23,25, Fru_3 : 4,9,12, Fru_4 : 3, Fru_5 : 5, Fru_6 : 7, Fru_7 : 8, Fru_8 : 10, Fru_9 : 11, Fru_10 : 13, 
#  Fru_11 : 15,22,  Fru_12 : 16,24, Fru_13 : 17, Fru_14 : 18, Fru_15 : 20

fru_1_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=c(0,1,19,21,26))
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_1_Int) <- 'fru_1'
fru_2_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=c(2,6,14,23,25))
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_2_Int) <- 'fru_2'
fru_3_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=c(4,9,12))
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_3_Int) <- 'fru_3'
fru_4_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=3)
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_4_Int) <- 'fru_4'
fru_5_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=5)
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_5_Int) <- 'fru_5'
fru_6_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=7)
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_6_Int) <- 'fru_6'
fru_7_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=8)
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_7_Int) <- 'fru_7'
fru_8_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=10)
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_8_Int) <- 'fru_8'
fru_9_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=11)
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_9_Int) <- 'fru_9'
fru_10_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=13)
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_10_Int) <- 'fru_10'
fru_11_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=c(15,22))
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_11_Int) <- 'fru_11'
fru_12_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=c(16,24))
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_12_Int) <- 'fru_12'
fru_13_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=17)
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_13_Int) <- 'fru_13'
fru_14_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=18)
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_14_Int) <- 'fru_14'
fru_15_Int <- WhichCells(list_subcluster_celltypes_fru[[5]], idents=20)
Idents(list_subcluster_celltypes_fru[[5]], cells = fru_15_Int) <- 'fru_15'

fru_anno_plot_merge <- plot_celltype_summary(plotlist=list_subcluster_celltypes_fru, marker_list=fru_markers, 
                                              ncols=5, pt_size=1, labelsize=5,
                                              exclude = NA, r_width=c(1,2)) 

cowplot::plot_grid(fru_anno_plot_merge[[5]],fru_anno_plot_merge[[1]],fru_anno_plot_merge[[2]],fru_anno_plot_merge[[3]],
                   ncol=1)
ggsave("Plots/celltype_groups/fru_sp_subtypes_summary.png", width=20, height=39)


## neuropeptide based grouping & annotation ##

### NP cell types ###

load("Processed_Data/NP_celltypes.RData")

NP_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                   individual=Trio_ISub_DF_list, celltype=NP_celltypes,
                                   color='red')
plot_grid_merge(NP_plots)
ggsave("Plots/celltype_groups/NP_plots.png", width=10, height=11)

list_subcluster_celltypes_NP <- Get_subcluster_celltypes_from_Int(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                                   celltype=NP_celltypes, res=1)

dimplot_list_NP <- plot_sp_dimplot_merge(list_subcluster_celltypes_NP)
plot_grid_merge(dimplot_list_NP)
ggsave("Plots/celltype_groups/NP_sp_subtypes.png", width=10, height=11)

for (i in 1:5) {
  
  DefaultAssay(list_subcluster_celltypes_NP[[i]]) <- "SCT"
  list_subcluster_celltypes_NP[[i]] <- PrepSCTFindMarkers(list_subcluster_celltypes_NP[[i]])
  
}

NP_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_NP[[5]], only.pos = T)

NP_Int_subtype_markers_filtered <- NP_Int_subtype_markers %>%
  dplyr::filter(gene %in% df_neuropep$gene)

NP_np_markers <- unique(NP_Int_subtype_markers_filtered$gene)

dotplot_list_NP <- plot_sp_dotplot_merge(list_subcluster_celltypes_NP, c('ChAT', 'VGlut', 'Gad1', 'pros','Imp', NP_np_markers))
plot_grid_merge(dotplot_list_NP)
ggsave("Plots/celltype_groups/NP_sp_subtypes_dotplots.png", width=26, height=20)

FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("sNPF"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("CCHa2"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("Mip"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("Dh44"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("Ms"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("Proc"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("Tk"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("ITP"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("spab"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("AstC"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("NPF"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("AstA"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("Dh31"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("FMRFa"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("RYa"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("Ilp7"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("CCHa1"), order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_NP[[5]], feature=c("Nplp1"), order=T) + NoAxes()

## sNPF : 0, CCHa2 : 1, Mip : 3, Tk_1 : 10, 27, Tk_2 : 12, TK_3 : 29, sNPF_CCHa2_AstC : 13, NPF_Dh31 : 14, Dh31 : 18, 42, NPF_AstA : 24, 
## AstC_FMRFa : 25, Tk_Mip_RYa : 26, RYa : 28, 44, Tk_sNPF : 32, FMRFa_Dh44 : 39, Ms : 40, FMRFa : 41, sNPF_Dh44 : 49

NP_markers_final <- c("sNPF","CCHa2","Mip","Tk","AstC","NPF","Dh31","AstA","FMRFa","RYa","Ms")

NP_NA <- setdiff(0:58, c(0,1,3,10,27,12,29,13,14,18,42,24,25,26,28,44,32,39,40,41,49))

sNPF_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=0)
Idents(list_subcluster_celltypes_NP[[5]], cells = sNPF_Int) <- 'sNPF'
CCHa2_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=1)
Idents(list_subcluster_celltypes_NP[[5]], cells = CCHa2_Int) <- 'CCHa2'
Mip_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=3)
Idents(list_subcluster_celltypes_NP[[5]], cells = Mip_Int) <- 'Mip'
Tk_1_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=c(10,27))
Idents(list_subcluster_celltypes_NP[[5]], cells = Tk_1_Int) <- 'Tk_1'
Tk_2_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=12)
Idents(list_subcluster_celltypes_NP[[5]], cells = Tk_2_Int) <- 'Tk_2'
Tk_3_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=29)
Idents(list_subcluster_celltypes_NP[[5]], cells = Tk_3_Int) <- 'Tk_3'
sNPF_CCHa2_AstC_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=13)
Idents(list_subcluster_celltypes_NP[[5]], cells = sNPF_CCHa2_AstC_Int) <- 'sNPF_CCHa2_AstC'
NPF_Dh31_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=14)
Idents(list_subcluster_celltypes_NP[[5]], cells = NPF_Dh31_Int) <- 'NPF_Dh31'
Dh31_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=c(18,42))
Idents(list_subcluster_celltypes_NP[[5]], cells = Dh31_Int) <- 'Dh31'
NPF_AstA_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=24)
Idents(list_subcluster_celltypes_NP[[5]], cells = NPF_AstA_Int) <- 'NPF_AstA'
AstC_FMRFa_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=25)
Idents(list_subcluster_celltypes_NP[[5]], cells = AstC_FMRFa_Int) <- 'AstC_FMRFa'
Tk_Mip_RYa_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=26)
Idents(list_subcluster_celltypes_NP[[5]], cells = Tk_Mip_RYa_Int) <- 'Tk_Mip_RYa'
RYa_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=c(28,44))
Idents(list_subcluster_celltypes_NP[[5]], cells = RYa_Int) <- 'RYa'
Tk_sNPF_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=32)
Idents(list_subcluster_celltypes_NP[[5]], cells = Tk_sNPF_Int) <- 'Tk_sNPF'
FMRFa_Dh44_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=39)
Idents(list_subcluster_celltypes_NP[[5]], cells = FMRFa_Dh44_Int) <- 'FMRFa_Dh44'
Ms_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=40)
Idents(list_subcluster_celltypes_NP[[5]], cells = Ms_Int) <- 'MS'
FMRFa_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=41)
Idents(list_subcluster_celltypes_NP[[5]], cells = FMRFa_Int) <- 'FMRFa'
sNPF_Dh44_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=49)
Idents(list_subcluster_celltypes_NP[[5]], cells = sNPF_Dh44_Int) <- 'sNPF_Dh44'
Unanno_NP_Int <- WhichCells(list_subcluster_celltypes_NP[[5]], idents=NP_NA)
Idents(list_subcluster_celltypes_NP[[5]], cells = Unanno_NP_Int) <- 'Unanno_NP'

DimPlot(list_subcluster_celltypes_NP[[4]], label=T) + NoLegend()

Idents(list_subcluster_celltypes_NP[[1]], cells = Unanno_NP_Int) <- 'Unanno_NP'
Idents(list_subcluster_celltypes_NP[[2]], cells = Unanno_NP_Int) <- 'Unanno_NP'
Idents(list_subcluster_celltypes_NP[[3]], cells = Unanno_NP_Int) <- 'Unanno_NP'
Idents(list_subcluster_celltypes_NP[[4]], cells = Unanno_NP_Int) <- 'Unanno_NP'

NP_anno_plot_merge <- plot_celltype_summary(plotlist=list_subcluster_celltypes_NP, marker_list=NP_markers_final, 
                                             ncols=4, pt_size=1, labelsize=5,
                                             exclude = 'Unanno_NP', r_width=c(2,3)) 

cowplot::plot_grid(NP_anno_plot_merge[[5]],NP_anno_plot_merge[[1]],NP_anno_plot_merge[[2]],NP_anno_plot_merge[[3]],
                   ncol=1)
ggsave("Plots/celltype_groups/NP_sp_subtypes_summary.png", width=20, height=39)


## Neurotransmitter based grouping & annotation ##

### Ach cell types ###

Ach_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                    individual=Trio_ISub_DF_list, celltype='Ach',
                                    color='red')
plot_grid_merge(Ach_plots)
ggsave("Plots/celltype_groups/Ach_plots.png", width=10, height=11)

list_subcluster_celltypes_Ach <- Get_subcluster_celltypes_from_Int(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                                    celltype='Ach', res=0.5)
dimplot_list_Ach <- plot_sp_dimplot_merge(list_subcluster_celltypes_Ach, pts=0.2)
plot_grid_merge(dimplot_list_Ach)
ggsave("Plots/celltype_groups/Ach_sp_subtypes.png", width=10, height=11)

for (i in 1:5) {
  
  DefaultAssay(list_subcluster_celltypes_Ach[[i]]) <- "SCT"
  list_subcluster_celltypes_Ach[[i]] <- PrepSCTFindMarkers(list_subcluster_celltypes_Ach[[i]])
  
}

Ach_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Ach[[5]], only.pos = T)
Ach_Dmel_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Ach[[1]], only.pos = T)
Ach_Dsec_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Ach[[3]], only.pos = T)

FeaturePlot(list_subcluster_celltypes_Ach[[5]], feature=c('pros', 'Imp'), order=T)
FeaturePlot(list_subcluster_celltypes_Ach[[1]], feature=c('pros', 'Imp'), order=T)
FeaturePlot(list_subcluster_celltypes_Ach[[2]], feature=c('pros', 'Imp'), order=T)
FeaturePlot(list_subcluster_celltypes_Ach[[3]], feature=c('pros', 'Imp'), order=T)

FeaturePlot(list_subcluster_celltypes_Ach[[5]], feature=c('heph','rsh'), order=F)
FeaturePlot(list_subcluster_celltypes_Ach[[5]], feature=c('ChAT'), order=T)

DotPlot(list_subcluster_celltypes_Ach[[5]], feature=c('pros', 'Imp', 'heph','rsh', 'CG14459',
                                                      'Lim1','acj6', 'ct','br', 
                                                      'Lim3','fkh','salm','Dh31',
                                                      'Sp1','Octbeta1R','sNPF','bab1', 
                                                      'Ret', 'Dh44', 'Sox102F','Ac78C',
                                                      'sv','SiaT','dac','Wnt4'))

Ach_markers <- c('ChAT', 'VAChT', 'pros', 'Imp', 'heph','rsh', 'CG14459',
                 'Lim1','acj6', 'ct','br', 
                 'Lim3','fkh','salm','Dh31',
                 'Sp1','Octbeta1R','sNPF','bab1', 
                 'Ret', 'Dh44', 'Sox102F','Ac78C',
                 'sv','SiaT','dac','Wnt4')
  
dotplot_list_Ach <- plot_sp_dotplot_merge(list_subcluster_celltypes_Ach, Ach_markers)
plot_grid_merge(dotplot_list_Ach)
ggsave("Plots/celltype_groups/Ach_sp_subtypes_dotplots.png", width=38, height=16)

#Ach_1 : 0,6,19,32,33, Ach_2 : 1,10,12,13,31, Ach_3 : 3,5,16, Ach_4 : 4,8,9,15, Ach_5 : 7, Ach_6 : 11, Ach_7 : 17, Ach_8 : 20, Ach_9 : 21, 
#Ach_10 : 22, Ach_11 : 23, Ach_12 : 24, Ach_13 : 26, Ach_14 : 27, Ach_15 : 28, Ach_16 : 29, Ach_17 : 30, Ach_unannotated : 2,14,18,25,34,35,36

Ach_1_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=c(0,6,19,32,33))
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_1_Int) <- 'Ach_1'
Ach_2_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=c(1,10,12,13,31))
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_2_Int) <- 'Ach_2'
Ach_3_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=c(3,5,16))
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_3_Int) <- 'Ach_3'
Ach_4_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=c(4,8,9,15))
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_4_Int) <- 'Ach_4'
Ach_5_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=7)
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_5_Int) <- 'Ach_5'
Ach_6_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=11)
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_6_Int) <- 'Ach_6'
Ach_7_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=17)
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_7_Int) <- 'Ach_7'
Ach_8_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=20)
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_8_Int) <- 'Ach_8'
Ach_9_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=21)
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_9_Int) <- 'Ach_9'
Ach_10_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=22)
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_10_Int) <- 'Ach_10'
Ach_11_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=23)
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_11_Int) <- 'Ach_11'
Ach_12_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=24)
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_12_Int) <- 'Ach_12'
Ach_13_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=26)
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_13_Int) <- 'Ach_13'
Ach_14_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=27)
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_14_Int) <- 'Ach_14'
Ach_15_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=28)
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_15_Int) <- 'Ach_15'
Ach_16_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=29)
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_16_Int) <- 'Ach_16'
Ach_17_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=30)
Idents(list_subcluster_celltypes_Ach[[5]], cells = Ach_17_Int) <- 'Ach_17'
Unanno_Ach_Int <- WhichCells(list_subcluster_celltypes_Ach[[5]], idents=c(2,14,18,25,34,35,36))
Idents(list_subcluster_celltypes_Ach[[5]], cells = Unanno_Ach_Int) <- 'Unanno_Ach'

Idents(list_subcluster_celltypes_Ach[[1]], cells = Unanno_Ach_Int) <- 'Unanno_Ach'
Idents(list_subcluster_celltypes_Ach[[2]], cells = Unanno_Ach_Int) <- 'Unanno_Ach'
Idents(list_subcluster_celltypes_Ach[[3]], cells = Unanno_Ach_Int) <- 'Unanno_Ach'
Idents(list_subcluster_celltypes_Ach[[4]], cells = Unanno_Ach_Int) <- 'Unanno_Ach'

DimPlot(list_subcluster_celltypes_Ach[[5]], label=T) + NoLegend()

Ach_anno_plot_merge <- plot_celltype_summary(plotlist=list_subcluster_celltypes_Ach, marker_list=Ach_markers, 
                                             ncols=5, pt_size=0.1, labelsize=5,
                                             exclude = 'Unanno_Ach', r_width=c(1,2)) 

cowplot::plot_grid(Ach_anno_plot_merge[[5]],Ach_anno_plot_merge[[1]],Ach_anno_plot_merge[[2]],Ach_anno_plot_merge[[3]],
                   ncol=1)
ggsave("Plots/celltype_groups/Ach_sp_subtypes_summary.png", width=20, height=39)


### Glu cell types ###

Glu_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                   individual=Trio_ISub_DF_list, celltype='Glu',
                                   color='red')
plot_grid_merge(Glu_plots)
ggsave("Plots/celltype_groups/Glu_plots.png", width=10, height=11)

list_subcluster_celltypes_Glu<- Get_subcluster_celltypes_from_Int(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                                   celltype='Glu', res=0.5)

DimPlot(list_subcluster_celltypes_Glu[[5]], label=T) + NoLegend()

dimplot_list_Glu <- plot_sp_dimplot_merge(list_subcluster_celltypes_Glu, pts=0.2)
plot_grid_merge(dimplot_list_Glu)
ggsave("Plots/celltype_groups/Glu_sp_subtypes.png", width=10, height=11)

for (i in 1:5) {
  
  DefaultAssay(list_subcluster_celltypes_Glu[[i]]) <- "SCT"
  list_subcluster_celltypes_Glu[[i]] <- PrepSCTFindMarkers(list_subcluster_celltypes_Glu[[i]])
  
}

Glu_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_Glu[[5]], only.pos = T)

FeaturePlot(list_subcluster_celltypes_Glu[[5]], feature=c('vvl','CG14459'), order=F) * NoAxes()
FeaturePlot(list_subcluster_celltypes_Glu[[5]], feature='bsh', order=F) + NoAxes()
FeaturePlot(list_subcluster_celltypes_Glu[[5]], feature='VGlut', order=T) + NoAxes()

DotPlot(list_subcluster_celltypes_Glu[[5]], feature=c('VGlut', 'pros', 'Imp','rsh','bi','CG14459','Dh31','br',
                                                      'Dh44','grn','otp','Sp1','DIP-lambda', 'nAChRalpha6','acj6',
                                                      'bsh','vvl','Lmx1a','SiaT', 'dve','Tk','lz','tey','Sox21b', 'disco-r', 'witty'))


Glu_markers <- c('VGlut', 'Imp', 'rsh', 'CG14459', 'Dh31', 'br', 'grn', 'Sp1', 'DIP-lambda', 'bsh', 'vvl', 'Lmx1a', 'SiaT', 'disco-r', 'lz')

dotplot_list_Glu <- plot_sp_dotplot_merge(list_subcluster_celltypes_Glu, Glu_markers)
plot_grid_merge(dotplot_list_Glu)
ggsave("Plots/celltype_groups/Glu_sp_subtypes_dotplots.png", width=25, height=20)

# Glu_1 : 1,7,23, Glu_2 : 0, Glu_3 : 2,4,10, Glu_4 : 5, Glu_5 : 6, Glu_6 : 8,13, Glu_7 : 14, Glu_8 : 15, 
# Glu_9 : 17, Glu_10 : 18, Glu_11 : 19, Glu_12 : 20, Glu_13 : 21, Glu_14 : 26, Glu_15 : 28

Glu_1_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=c(1,7,23))
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_1_Int) <- 'Glu_1'
Glu_2_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=0)
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_2_Int) <- 'Glu_2'
Glu_3_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=c(2,4,10))
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_3_Int) <- 'Glu_3'
Glu_4_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=5)
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_4_Int) <- 'Glu_4'
Glu_5_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=6)
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_5_Int) <- 'Glu_5'
Glu_6_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=c(8,13))
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_6_Int) <- 'Glu_6'
Glu_7_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=14)
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_7_Int) <- 'Glu_7'
Glu_8_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=15)
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_8_Int) <- 'Glu_8'
Glu_9_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=17)
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_9_Int) <- 'Glu_9'
Glu_10_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=18)
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_10_Int) <- 'Glu_10'
Glu_11_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=19)
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_11_Int) <- 'Glu_11'
Glu_12_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=20)
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_12_Int) <- 'Glu_12'
Glu_13_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=21)
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_13_Int) <- 'Glu_13'
Glu_14_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=26)
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_14_Int) <- 'Glu_14'
Glu_15_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=28)
Idents(list_subcluster_celltypes_Glu[[5]], cells = Glu_15_Int) <- 'Glu_15'

Glu_NA <- setdiff(0:42, c(1,7,23,0,2,4,10,5,6,8,13,14,15,17,18,19,20,21,26,28))
Unanno_Glu_Int <- WhichCells(list_subcluster_celltypes_Glu[[5]], idents=Glu_NA)
Idents(list_subcluster_celltypes_Glu[[5]], cells = Unanno_Glu_Int) <- 'Unanno_Glu'

Idents(list_subcluster_celltypes_Glu[[1]], cells = Unanno_Glu_Int) <- 'Unanno_Glu'
Idents(list_subcluster_celltypes_Glu[[2]], cells = Unanno_Glu_Int) <- 'Unanno_Glu'
Idents(list_subcluster_celltypes_Glu[[3]], cells = Unanno_Glu_Int) <- 'Unanno_Glu'
Idents(list_subcluster_celltypes_Glu[[4]], cells = Unanno_Glu_Int) <- 'Unanno_Glu'

DimPlot(list_subcluster_celltypes_Glu[[5]], label=T) + NoLegend()

Glu_anno_plot_merge <- plot_celltype_summary(plotlist=list_subcluster_celltypes_Glu, marker_list=Glu_markers, 
                                             ncols=5, pt_size=0.1, labelsize=5,
                                             exclude = 'Unanno_Glu', r_width=c(1,2)) 

cowplot::plot_grid(Glu_anno_plot_merge[[5]],Glu_anno_plot_merge[[1]],Glu_anno_plot_merge[[2]],Glu_anno_plot_merge[[3]],
                   ncol=1)
ggsave("Plots/celltype_groups/Glu_sp_subtypes_summary.png", width=20, height=39)



### GABA cell types ###

GABA_plots <- plot_cell_type_groups(integrated=TrioBrain.integrated_slim_ISub_DF_labeled_species, 
                                   individual=Trio_ISub_DF_list, celltype='GABA',
                                   color='red')
plot_grid_merge(GABA_plots)
ggsave("Plots/celltype_groups/GABA_plots.png", width=10, height=11)

list_subcluster_celltypes_GABA <- Get_subcluster_celltypes_from_Int(integrated=TrioBrain.integrated_slim_ISub_DF_labeled, 
                                                                    celltype='GABA', res=0.5)
dimplot_list_GABA <- plot_sp_dimplot_merge(list_subcluster_celltypes_GABA, pts=0.2)
plot_grid_merge(dimplot_list_GABA)
ggsave("Plots/celltype_groups/GABA_sp_subtypes.png", width=10, height=11)

for (i in 1:5) {
  
  DefaultAssay(list_subcluster_celltypes_GABA[[i]]) <- "SCT"
  list_subcluster_celltypes_GABA[[i]] <- PrepSCTFindMarkers(list_subcluster_celltypes_GABA[[i]])
  
}

GABA_Int_subtype_markers <- FindAllMarkers(list_subcluster_celltypes_GABA[[5]], only.pos = T)

FeaturePlot(list_subcluster_celltypes_GABA[[5]], feature='DIP-alpha', order=T) + NoAxes()
FeaturePlot(list_subcluster_celltypes_GABA[[5]], feature='Gad1', order=T) + NoAxes()

DotPlot(list_subcluster_celltypes_GABA[[5]], feature=c('Gad1','heph','rsh', 'grn','Scr','trh','pb','tey','fkh','ome', 'acj6', 'inv', 'Lim3', 'Nep1','rk'))

GABA_markers <- c('Gad1', 'heph', 'rsh', 'Scr', 'pb', 'fkh', 'ome', 'acj6', 'inv', 'Lim3', 'trh', 'ey', 'Ets65A', 'chas', 'rk')

dotplot_list_GABA <- plot_sp_dotplot_merge(list_subcluster_celltypes_GABA, GABA_markers)
plot_grid_merge(dotplot_list_GABA)
ggsave("Plots/celltype_groups/GABA_sp_subtypes_dotplots.png", width=20, height=15)

# GABA_1 : 0,22, GABA_2 : 2,6, GABA_3 : 3, GABA_4 : 4,8,12, GABA_5 : 7, GABA_6 : 9, GABA_7 : 10, 
# GABA_8 : 13,18, GABA_9 : 14, GABA_10 : 15, GABA_11 : 16, GABA_12 : 17, GABA_13 : 19, GABA_14 : 20

GABA_1_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=c(0,22))
Idents(list_subcluster_celltypes_GABA[[5]], cells = GABA_1_Int) <- 'GABA_1'
GABA_2_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=c(2,6))
Idents(list_subcluster_celltypes_GABA[[5]], cells = GABA_2_Int) <- 'GABA_2'
GABA_3_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=3)
Idents(list_subcluster_celltypes_GABA[[5]], cells = GABA_3_Int) <- 'GABA_3'
GABA_4_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=c(4,8,12))
Idents(list_subcluster_celltypes_GABA[[5]], cells = GABA_4_Int) <- 'GABA_4'
GABA_5_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=7)
Idents(list_subcluster_celltypes_GABA[[5]], cells = GABA_5_Int) <- 'GABA_5'
GABA_6_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=9)
Idents(list_subcluster_celltypes_GABA[[5]], cells = GABA_6_Int) <- 'GABA_6'
GABA_7_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=10)
Idents(list_subcluster_celltypes_GABA[[5]], cells = GABA_7_Int) <- 'GABA_7'
GABA_8_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=c(13,18))
Idents(list_subcluster_celltypes_GABA[[5]], cells = GABA_8_Int) <- 'GABA_8'
GABA_9_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=14)
Idents(list_subcluster_celltypes_GABA[[5]], cells = GABA_9_Int) <- 'GABA_9'
GABA_10_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=15)
Idents(list_subcluster_celltypes_GABA[[5]], cells = GABA_10_Int) <- 'GABA_10'
GABA_11_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=16)
Idents(list_subcluster_celltypes_GABA[[5]], cells = GABA_11_Int) <- 'GABA_11'
GABA_12_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=17)
Idents(list_subcluster_celltypes_GABA[[5]], cells = GABA_12_Int) <- 'GABA_12'
GABA_13_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=19)
Idents(list_subcluster_celltypes_GABA[[5]], cells = GABA_13_Int) <- 'GABA_13'
GABA_14_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=20)
Idents(list_subcluster_celltypes_GABA[[5]], cells = GABA_14_Int) <- 'GABA_14'

GABA_NA <- setdiff(0:34, c(0,22,2,6,3,4,8,12,7,9,10,13,18,14,15,16,17,19,20))
Unanno_GABA_Int <- WhichCells(list_subcluster_celltypes_GABA[[5]], idents=GABA_NA)
Idents(list_subcluster_celltypes_GABA[[5]], cells = Unanno_GABA_Int) <- 'Unanno_GABA'

Idents(list_subcluster_celltypes_GABA[[1]], cells = Unanno_GABA_Int) <- 'Unanno_GABA'
Idents(list_subcluster_celltypes_GABA[[2]], cells = Unanno_GABA_Int) <- 'Unanno_GABA'
Idents(list_subcluster_celltypes_GABA[[3]], cells = Unanno_GABA_Int) <- 'Unanno_GABA'
Idents(list_subcluster_celltypes_GABA[[4]], cells = Unanno_GABA_Int) <- 'Unanno_GABA'

DimPlot(list_subcluster_celltypes_GABA[[5]], label=T) + NoLegend()

GABA_anno_plot_merge <- plot_celltype_summary(plotlist=list_subcluster_celltypes_GABA, marker_list=GABA_markers, 
                                             ncols=5, pt_size=0.5, labelsize=5,
                                             exclude = 'Unanno_GABA', r_width=c(1,2)) 

cowplot::plot_grid(GABA_anno_plot_merge[[5]],GABA_anno_plot_merge[[1]],GABA_anno_plot_merge[[2]],GABA_anno_plot_merge[[3]],
                   ncol=1)
ggsave("Plots/celltype_groups/GABA_sp_subtypes_summary.png", width=20, height=39)


## Final summary - total 107 annotated cell types

save(list_subcluster_celltypes_glia, list_subcluster_celltypes_KC, list_subcluster_celltypes_MA, 
     list_subcluster_celltypes_Clock, list_subcluster_celltypes_Poxn, 
     list_subcluster_celltypes_OPN,
     list_subcluster_celltypes_fru, list_subcluster_celltypes_NP, list_subcluster_celltypes_Ach, 
     list_subcluster_celltypes_Glu, list_subcluster_celltypes_GABA,
     file="Processed_Data/cell_type_lists.RData")


