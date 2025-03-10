library(ggplot2)
library(dplyr)
#devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(Seurat)

load("Processed_Data/cell_info_join_Trio.RData")

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

## cluster comparison ##

df_Dmel_celltypes <- as.character(unique(filter(cell_info_join_Trio, orig=='Dmel')$Cluster.y))
df_Dsim_celltypes <- as.character(unique(filter(cell_info_join_Trio, orig=='Dsim')$Cluster.y))
df_Dsec_celltypes <- as.character(unique(filter(cell_info_join_Trio, orig=='Dsec')$Cluster.y))
df_DsecNoni_celltypes <- as.character(unique(filter(cell_info_join_Trio, orig=='DsecNoni')$Cluster.y))
df_integrated_celltypes <- as.character(unique(cell_info_join_Trio$Cluster.x))

## all cells ##

df_test1 <- cell_info_join_Trio %>%
  dplyr::filter(orig %in% c("Dmel","Dsec") & Cluster.x != "Doublets" & Cluster.y != "Doublets") %>%
  #dplyr::filter(xy_x > 0.75) %>%
  dplyr::select(Cell, Cluster.x, Cluster.y, orig) %>%
  tidyr::spread(key='orig', value='Cluster.y') %>%
  dplyr::select(Dmel, Integrated=Cluster.x, Dsec)

df_plot <- df_test1 %>%
  make_long(Dsec,Integrated, 
            Dmel) %>%
  dplyr::filter(if_any(c(node, next_x, next_node), ~ !is.na(.))) %>%
  dplyr::filter(x!='Dsec' | !is.na(node))

#df_plot$node <- factor(df_plot$node, levels=order)
#df_plot$next_node <- factor(df_plot$next_node, levels=c(order_nextnode,NA))

pl <- ggplot(df_plot, aes(x = x,                        
                          next_x = next_x,                                     
                          node = node,
                          next_node = next_node,        
                          fill = factor(node))) +
  geom_sankey(flow.alpha = 0.5,          #This Creates the transparency of your node 
              node.color = 'black',     # This is your node color  - set NA for transparency       
              space = 0,
              show.legend = F,
              na.rm=T)  +
  geom_sankey_text(aes(label = node), 
                   space = 0,
                   size = 3,
                   na.rm=T) +  # Adjust size and vjust as needed
  theme_minimal()

pl

ggsave("Plots/melsec_sankey.pdf", width=20, height=45)


## glial cells only ##



df_glia <- cell_info_join_Trio %>%
  dplyr::filter(orig %in% c("Dmel","Dsec") & (Cluster.x %in% glial_celltypes | Cluster.y %in% glial_celltypes )) %>%
  #dplyr::filter(xy_x > 0.75) %>%
  dplyr::select(Cell, Cluster.x, Cluster.y, orig) %>%
  tidyr::spread(key='orig', value='Cluster.y') %>%
  dplyr::select(Dmel, Integrated=Cluster.x, Dsec)

df_plot_glia <- df_glia %>%
  make_long(Dsec,Integrated, 
            Dmel) %>%
  dplyr::filter(if_any(c(node, next_x, next_node), ~ !is.na(.))) %>%
  dplyr::filter(x!='Dsec' | !is.na(node))

pl_glia <- ggplot(df_plot_glia, aes(x = x,                        
                          next_x = next_x,                                     
                          node = node,
                          next_node = next_node,        
                          fill = factor(node))) +
  geom_sankey(flow.alpha = 0.5,          #This Creates the transparency of your node 
              node.color = 'black',     # This is your node color  - set NA for transparency       
              space = 0,
              show.legend = F,
              na.rm=T)  +
  geom_sankey_text(aes(label = node), 
                   space = 0,
                   size = 3,
                   na.rm=T) +  # Adjust size and vjust as needed
  theme_minimal()

pl_glia

ggsave("Plots/melsec_glia_sankey.pdf", width=10, height=20)
