library(ggplot2)
library(dplyr)
library(ggsankey) #devtools::install_github("davidsjoberg/ggsankey")
library(readr)
library(ggplot2)
library(ggrepel)
library(sctransform)
library(cowplot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/cell_info_join_all.RData")

TrioBrain.integrated_slim_ISub_DF_labeled <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_ISub_DF_labeled.rds")
Trio_ISub_DF_list_labeled <- readRDS(file = "Processed_Data/Trio_ISub_DF_list_labeled.rds")
secsim_DmelRef_ISub_DF_list_labeled <- readRDS(file = "Processed_Data/secsim_DmelRef_ISub_DF_list_labeled.rds")

DimPlot(subset(TrioBrain.integrated_slim_ISub_DF_labeled, ClusterLabel != 'Doublets'), 
        label=T, reduction = 'tsne', repel=T, group.by='ClusterLabel', raster=F) + NoLegend() + NoAxes() +
  ggtitle("")

ggsave("Plots/TrioBrain.integrated_slim_ISub_DF_labeled.png", width=14, height = 12)

species_vector <- c("Dmel", "Dsim", "Dsec", "DsecNoni", "Dsim_to_DmelRef", "Dsec_to_DmelRef", "DsecNoni_to_DmelRef")

for (i in 1:length(species_vector)) {
  
  if (i<=4) {
    
    DimPlot(subset(Trio_ISub_DF_list_labeled[[i]], ClusterLabel != 'Doublets'), 
            label=T, reduction = 'tsne', repel=T, group.by='ClusterLabel', raster=F) + NoLegend() + NoAxes() +
      ggtitle("")
    
    ggsave(glue::glue("Plots/{species_vector[i]}_ISub_DF_list_labeled.png"), width=9.8, height = 8.4)
    
  }
  
  else if (i>4) {
    
    DimPlot(subset(secsim_DmelRef_ISub_DF_list_labeled[[i-4]], ClusterLabel != 'Doublets'), 
            label=T, reduction = 'tsne', repel=T, group.by='ClusterLabel', raster=F) + NoLegend() + NoAxes() +
      ggtitle("")
    
    ggsave(glue::glue("Plots/{species_vector[i]}_ISub_DF_list_labeled.png"), width=9.8, height = 8.4)
    
  }

  
}

## cluster comparison ##

df_Dmel_celltypes <- as.character(unique(filter(cell_info_join_Trio, orig=='Dmel')$Cluster.y))
df_Dsim_celltypes <- as.character(unique(filter(cell_info_join_Trio, orig=='Dsim')$Cluster.y))
df_Dsec_celltypes <- as.character(unique(filter(cell_info_join_Trio, orig=='Dsec')$Cluster.y))
df_DsecNoni_celltypes <- as.character(unique(filter(cell_info_join_Trio, orig=='DsecNoni')$Cluster.y))
df_integrated_celltypes <- as.character(unique(cell_info_join_Trio$Cluster.x))

cell_info_join_Trio_summary <- cell_info_join_Trio %>%
  dplyr::distinct(Cluster.x, orig, Cluster.y, xy_x, xy_y)

cell_info_join_Trio_high_con <- cell_info_join_Trio %>%
  dplyr::filter(xy_x >= 0.8 & xy_y >= 0.1, orig != "DsecNoni") %>%
  dplyr::distinct(Cluster.x, orig, Cluster.y, xy_x, xy_y) %>%
  dplyr::group_by(Cluster.x) %>%
  dplyr::mutate(n_species = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(-n_species, Cluster.x)

cell_info_join_Trio_high_con_selected <- cell_info_join_Trio_high_con %>%
  dplyr::filter(n_species>=2)

conserved_celltype <- as.character(unique(cell_info_join_Trio_high_con_selected$Cluster.x))
conserved_celltype

count_high_con_selected_trio <- cell_info_join_Trio_high_con %>%
  dplyr::filter(orig %in% c("Dmel","Dsim", "Dsec")) %>%
  dplyr::filter(n_species==3) %>%
  dplyr::distinct(Cluster.x) %>%
  pull(Cluster.x) %>%
  length()

cell_info_join_Trio_high_con %>%
  dplyr::filter(orig %in% c("Dmel","Dsim", "Dsec")) %>%
  dplyr::filter(n_species==2) %>%
  dplyr::group_by(Cluster.x) %>%
  dplyr::summarise(pair = paste(orig, collapse = "")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(pair) %>%
  dplyr::summarise(n_pair=n()) %>%
  dplyr::ungroup()

## type-specific sankey plots ##

generate_sankey_celltype()

target_ct=KC_celltypes

generate_sankey_celltype <- function(df=cell_info_join_Trio, target_ct=glial_celltypes, target_sp=c("Dmel","Dsec"), filter_threshold =0.05) {
  
  df_target_prep <- cell_info_join_Trio %>%
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
  
  df_target <- cell_info_join_Trio %>%
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
         nrow = length(glial_celltypes), 
         byrow = TRUE)[t(sapply(target_ct, function(ct) grepl(ct, df_target_cluster_sp1)))])
  
  filter_cluster_Integrated <- as.character(matrix(df_target_cluster_Integrated, 
                                            ncol = length(df_target_cluster_Integrated), 
                                            nrow = length(glial_celltypes), 
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

# glia sankeyplot ##

glial_celltypes <- c("ENS","AST","PRN","SUB", "CTX")
pl_glia_melsec <- generate_sankey_celltype(df=cell_info_join_Trio, target_ct=glial_celltypes, target_sp=c("Dmel","Dsec"))
pl_glia_melsec

ggsave("Plots/sankey/pl_glia_melsec.png", width=7, height=19)

pair_list <- list(c("Dmel","Dsec"),c("Dmel","Dsim"), c("Dsim","Dsec"), c("Dsec", "DsecNoni"))

for (i in 1:4) {
  
  pair <- pair_list[[i]]
  
  sp1 <- pair[1]
  sp2 <- pair[2]
  
  generate_sankey_celltype(df=cell_info_join_Trio, target_ct=glial_celltypes, target_sp=pair)
  ggsave(glue::glue("Plots/sankey/pl_glia_{sp1}{sp2}.png"), width=12, height=12)
  
}

# KC sankeyplot ##

KC_celltypes <- c("αβ-KC","γ-KC","α'β'-KC")

for (i in 1:4) {
  
  pair <- pair_list[[i]]
  
  sp1 <- pair[1]
  sp2 <- pair[2]
  
  generate_sankey_celltype(df=cell_info_join_Trio, target_ct=KC_celltypes, target_sp=pair)
  ggsave(glue::glue("Plots/sankey/pl_KC_{sp1}{sp2}.png"), width=12, height=12)
  
}

## MA sankeyplot ##

MA_celltypes <- c("MON","TY","SER","DOP", "OCTY")

for (i in 1:4) {
  
  pair <- pair_list[[i]]
  
  sp1 <- pair[1]
  sp2 <- pair[2]
  
  generate_sankey_celltype(df=cell_info_join_Trio, target_ct=MA_celltypes, target_sp=pair)
  ggsave(glue::glue("Plots/sankey/pl_MA_{sp1}{sp2}.png"), width=12, height=12)
  
}

## KD sankeyplot ##

KD_celltypes <- c("Fru","Poxn","Clock","OPN")

for (ct in KD_celltypes) {
  
  for (i in 1:4) {
    
    pair <- pair_list[[i]]
    
    sp1 <- pair[1]
    sp2 <- pair[2]
    
    generate_sankey_celltype(df=cell_info_join_Trio, target_ct=ct, target_sp=pair)
    ggsave(glue::glue("Plots/sankey/pl_{ct}_{sp1}{sp2}.png"), width=12, height=12)
    
  }

}

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

for (i in 1:4) {
  
  pair <- pair_list[[i]]
  
  sp1 <- pair[1]
  sp2 <- pair[2]
  
  generate_sankey_celltype(df=cell_info_join_Trio, target_ct=NP_celltypes, target_sp=pair)
  ggsave(glue::glue("Plots/sankey/pl_NP_{sp1}{sp2}.png"), width=12, height=12)
  
}


## individual nps ## - need to fix the code

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

list_np_markers <- NULL

for (i in 1:4) {
  
  np_markers_temp <- get_marker_nps(seurat_object=Trio_ISub_DF_list_labeled[[i]], df_marker_sigs=Trio_ISub_DF_sig_marker_list[[i]])
  
  list_np_markers <- unique(c(list_np_markers, np_markers_temp))
  
}

pair_list_id <- list(c(1,3), c(1,2), c(2,3), c(3,4))

for (i in 1:4) {
  
  pair <- pair_list[[i]]
  
  sp1 <- pair[1]
  sp2 <- pair[2]
  
  sp1_id <- pair_list_id[[i]][1]
  sp2_id <- pair_list_id[[i]][2]
  
  np_markers_sp1 <- get_marker_nps(seurat_object=Trio_ISub_DF_list_labeled[[sp1_id]], df_marker_sigs=Trio_ISub_DF_sig_marker_list[[sp1_id]])
  np_markers_sp2 <- get_marker_nps(seurat_object=Trio_ISub_DF_list_labeled[[sp2_id]], df_marker_sigs=Trio_ISub_DF_sig_marker_list[[sp2_id]])
  
  np_shared_markers <- intersect(np_markers_sp1, np_markers_sp2)
  
  for (ct in np_shared_markers) {
    
    grepl(cell_info_join_Trio$Cluster.x)
    
    generate_sankey_celltype(df=cell_info_join_Trio, target_ct=ct, target_sp=pair)
    ggsave(glue::glue("Plots/sankey/NP/pl_{ct}_{sp1}{sp2}.png"), width=12, height=12)
    
  }
  
}

## highly conserved cell types - robust mapping across datasets ##



generate_sankey_matching <- function(df=cell_info_join_Trio, matching_frac1=0.9, matching_frac2=0.1, target_sp=c("Dmel","Dsec"), rm_doublets=T) {
  
  df_target1 <- cell_info_join_Trio %>%
    filter(orig %in% target_sp[1], xy_x >= matching_frac1 & xy_y >= matching_frac2) %>%
    distinct(Cluster.x, Cluster.y, .keep_all = T)
  
  Cluster_x1 <- df_target1$Cluster.x
  Cluster_y1 <- df_target1$Cluster.y
  
  df_target2 <- cell_info_join_Trio %>%
    filter(orig %in% target_sp[2], xy_x >= matching_frac1 & xy_y >= matching_frac2) %>%
    distinct(Cluster.x, Cluster.y, .keep_all = T)
  
  Cluster_x2 <- df_target2$Cluster.x
  Cluster_y2 <- df_target2$Cluster.y  
  
  df_target <- cell_info_join_Trio %>%
    filter(orig %in% target_sp) %>%
    dplyr::filter(ifelse(orig==target_sp[1], Cluster.x %in% c(Cluster_x1,Cluster_x2) | Cluster.y %in% Cluster_y1, 
                         Cluster.x %in% c(Cluster_x1,Cluster_x2) | Cluster.y %in% Cluster_y2)) %>%
    dplyr::select(Cell, Cluster.x, Cluster.y, orig) %>%
    tidyr::spread(key='orig', value='Cluster.y') %>%
    dplyr::select(sp1=target_sp[1], Integrated=Cluster.x, sp2=target_sp[2])
  
  if (rm_doublets==T) {
    
    df_target <- df_target %>%
      dplyr::filter((sp1 != "Doublets" & is.na(sp2) & Integrated != "Doublets") | (sp2 != "Doublets" & is.na(sp1) & Integrated != "Doublets"))
    
  }
  
  df_plot_target <- df_target %>%
    make_long(sp1,Integrated, 
              sp2) %>%
    dplyr::filter(if_any(c(node, next_x, next_node), ~ !is.na(.))) %>%
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
                na.rm=T)  +
    geom_sankey_label(aes(label = node), 
                      space = 0,
                      size = 3,
                      show.legend = F,
                      na.rm=T) +  # Adjust size and vjust as needed
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          panel.grid = element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_text(size=12, color='black'))
  
  pl_target
  
  return(pl_target)
  
}

generate_sankey_matching()

sankey_matching_095 <- generate_sankey_matching(df=cell_info_join_Trio, matching_frac1=0.95, matching_frac2=0.95, target_sp=c("Dmel","Dsec"), rm_doublets=F)
sankey_matching_095

ggsave("Plots/sankey_matching_095.pdf", width=10, height=20)

KD_celltypes <- c("Fru","OPN","Poxn","Clock")

pl_clu9_3_melsec <- generate_sankey_celltype(df=cell_info_join_Trio, target_ct='cluster9_3', target_sp=c("Dmel","Dsec"))
pl_clu9_3_melsec

pl_Ilp6_melsec <- generate_sankey_celltype(df=cell_info_join_Trio, target_ct='Ilp6', target_sp=c("Dmel","Dsec"))
pl_Ilp6_melsec

pl_Doublets_melsec <- generate_sankey_celltype(df=cell_info_join_Trio, target_ct='Doublets', target_sp=c("Dmel","Dsec"))
pl_Doublets_melsec
generate_sankey_celltype(df=cell_info_join_Trio, target_ct='Doublets', target_sp=c("Dsec","Dmel"))

melsim_ct_intersect <- intersect(df_Dmel_celltypes, df_Dsim_celltypes)[!grepl("Glu|GABA|Ach|cluster", intersect(df_Dmel_celltypes, df_Dsim_celltypes))]
pl_melsim_intersect <- generate_sankey_celltype(df=cell_info_join_Trio, target_ct=melsim_ct_intersect, target_sp=c("Dmel","Dsec"))
pl_melsim_intersect





generate_sankey_matching(df=cell_info_join_Trio, matching_frac1=0.8, matching_frac2=0.8, target_sp=c("Dmel","Dsec"), rm_doublets=T)

cell_info_join_Trio_high_con <- cell_info_join_Trio %>%
  dplyr::filter(xy_x >= 1 & xy_y >= 1)

as.character(unique(cell_info_join_Trio_high_con$Cluster.x))
as.character(unique(cell_info_join_Trio_high_con$Cluster.y))

high_con_ct <- unique(as.character(unique(cell_info_join_Trio_high_con$Cluster.y))) ### need to fix somehow

#high_con_ct[grepl('cluster',high_con_ct)]

pl_high_con_melsec <- generate_sankey_celltype(df=cell_info_join_Trio_high_con, target_ct=high_con_ct, target_sp=c("Dmel","Dsec"))
pl_high_con_melsec

