library(tidyverse)
library(Seurat)
library(ggrepel)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_sum_exp_SCT.RData")
load("Processed_Data/shared_genes.RData")
load("Processed_Data/celltype_order.RData")
TrioBrain.integrated_slim_labeled <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled.rds")

cc_genes <- c("VAChT", "ChAT", "mAChR-A", "mAChR-B", "mAChR-C", ## Ach signaling
              "nAChRalpha1", "nAChRalpha2", "nAChRalpha3","nAChRalpha5",  "nAChRalpha6", "nAChRalpha7", 
              "nAChRbeta1","nAChRbeta2", "nAChRbeta3", 
              "VGlut", "mGluR", "GluRIA", "GluRIB", "GluRIIE", "GluClalpha", ## Glu signaling
              "Gad1", "VGAT","Rdl", "Lcch3", "Grd", "CG8916", "GABA-B-R1","GABA-B-R2","GABA-B-R3", ## GABA
              "Vmat", "ple", "DAT", "Dop1R1","Dop1R2", "Dop2R", "DopEcR",  ## MON-Dop signaling
              "Tdc2", "Tbh","Oct-TyrR", "TyrR", "TyrRII", "Octalpha2R", "Octbeta1R",  "Octbeta2R",  "Octbeta3R", "Oamb", ## Tyr-Oct signaling
              "Trh", "SerT",   "5-HT1A","5-HT1B",  "5-HT2A", "5-HT2B", "5-HT7",  ## Ser signaling
              "Hdc",  "ort", "HisCl1")

cluster_list <- CellType_order$CellType

### normalized expression level

df_scale_cc <- df_sum_exp_SCT %>%
  dplyr::filter(species != "DsecNoni", gene %in% cc_genes) %>%
  dplyr::group_by(gene, cluster) %>%
  dplyr::summarise(median_expression = median(expression)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(max_expression=max(median_expression)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(max_expression > 0.3) %>%
  dplyr::mutate(norm_med_exp = median_expression/max_expression) 

df_scale_cc %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  ggplot(.) +
  geom_tile(color='black', size=0.1) +
  aes(y=factor(gene, levels=rev(cc_genes)), x=factor(cluster, levels=anno_cluster), fill=norm_med_exp) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y= element_text(color='black', face='italic', size=9),
        axis.text.x = element_text(color='black', size=9, angle=45, hjust = 1, vjust = 1),
        legend.title = element_text(size=9.5, color='black'),
        legend.text = element_text(size=9, color = 'black'),
        axis.title=element_blank(),
        legend.position='bottom',
        strip.text.x = element_text(margin = margin(0.02,0,0.02,0, "in"))) +
  scale_fill_gradient(low="white", high="darkblue", breaks = c(0, 0.2,0.4,0.6,0.8,1)) +
  #scale_fill_manual(values=c("white","darkblue")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(fill="expression")

ggsave(file="Plots/Manuscript/Fig2_heatmap_ccgene_cluster_exp.pdf", width=10, height= 8)

