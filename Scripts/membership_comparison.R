library(Seurat)
library(tidyverse)
library(readr)
library(ggplot2)
library(ggrepel)
library(Polychrome)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/celltype_order.RData")
TrioBrain.integrated_slim_labeled <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled.rds")

## set col palette ##

set.seed(723451) # for reproducibility
col_pal_55 <- createPalette(55, c("#ff0000"), M=1000)
names(col_pal_55) <- CellType_order$CellType

## Dmel ##

Dmel_decon_combined_DF <- readRDS(file = "Processed_Data/Dmel_combined_full_intron_decon_DF.rds")

DimPlot(Dmel_decon_combined_DF, reduction = "umap",  label = TRUE, pt.size = 0.01)

DotPlot(Dmel_decon_combined_DF, 
        dot.scale = 4.5,
        feature=c("VAChT","VGlut", "Gad1", "Vmat", "DAT", "Trh", "Tdc2", "Tbh", 
                  "nAChRalpha6","nAChRalpha5","nAChRalpha7", "ort",
                  "axo","Gs2","Tret1-1", "baz", "CG40470",
                  "crb","CG32204","Pka-C1", "mamo",
                  "sNPF", "CCHa2",  "Proc", "Mip", "AstA", "RYa", "Hug", "Ilp2","Ilp3","Ilp5","Tk", "SIFa", "Crz", 
                  "Dh31", "Dh44",
                  "Oaz", "SiaT", "fru", "tim", "Clk", "Poxn", "fit", "Hml")) +
  coord_flip() + labs(x=NULL, y="Clusters") +
  theme(legend.position='bottom',
        axis.title = element_blank())

## Dsim ##

Dsim_decon_combined_DF <- readRDS(file = "Processed_Data/Dsim_combined_full_intron_decon_DF.rds")

DimPlot(Dsim_decon_combined_DF, reduction = "umap",  label = TRUE, pt.size = 0.01)

DotPlot(Dsim_decon_combined_DF, 
        dot.scale = 4.5,
        feature=c("VAChT","VGlut", "Gad1", "Vmat", "DAT", "Trh", "Tdc2", "Tbh", 
                  "nAChRalpha6","nAChRalpha5","nAChRalpha7", "ort",
                  "axo","Gs2","Tret1-1", "baz", "CG40470",
                  "crb","CG32204","Pka-C1", "mamo",
                  "sNPF", "CCHa2",  "Proc", "Mip", "AstA", "RYa", "Hug", "Ilp2","Ilp3","Ilp5","Tk", "SIFa", "Crz", 
                  "Dh31", "Dh44",
                  "Oaz", "SiaT", "fru", "tim", "Clk", "Poxn", "fit", "Hml")) +
  coord_flip() + labs(x=NULL, y="Clusters") +
  theme(legend.position='bottom',
        axis.title = element_blank())

## Dsec ##

Dsec_decon_combined_DF <- readRDS(file = "Processed_Data/Dsec_combined_full_intron_decon_DF.rds")

DimPlot(Dsec_decon_combined_DF, reduction = "umap",  label = TRUE, pt.size = 0.01)

DotPlot(Dsec_decon_combined_DF, 
        dot.scale = 4.5,
        feature=c("VAChT","VGlut", "Gad1", "Vmat", "DAT", "Trh", "Tdc2", "Tbh", 
                  "nAChRalpha6","nAChRalpha5","nAChRalpha7", "ort",
                  "axo","Gs2","Tret1-1", "baz", "CG40470",
                  "crb","CG32204","Pka-C1", "mamo",
                  "sNPF", "CCHa2",  "Proc", "Mip", "AstA", "RYa", "Hug", "Ilp2","Ilp3","Ilp5","Tk", "SIFa", "Crz", 
                  "Dh31", "Dh44",
                  "Oaz", "SiaT", "fru", "tim", "Clk", "Poxn", "fit", "Hml")) +
  coord_flip() + labs(x=NULL, y="Clusters") +
  theme(legend.position='bottom',
        axis.title = element_blank())

### cross-comparison ###

mel_ID <- Idents(Dmel_decon_combined_DF)
sim_ID <- Idents(Dsim_decon_combined_DF)
sec_ID <- Idents(Dsec_decon_combined_DF)

Integrated_ID <- Idents(TrioBrain.integrated_slim_labeled)

Integrated_ID_mel <- Integrated_ID[grep("1_1", names(Integrated_ID))]
Integrated_ID_sim <- Integrated_ID[grep("1_2", names(Integrated_ID))]
Integrated_ID_sec <- Integrated_ID[grep("1_3", names(Integrated_ID))]

head(mel_ID)
head(Integrated_ID_mel)
head(sim_ID)
head(Integrated_ID_sim)
head(sec_ID)
head(Integrated_ID_sec)

names(Integrated_ID_mel) <- gsub("1_1", "1", names(Integrated_ID_mel))
names(Integrated_ID_sim) <- gsub("1_2", "1", names(Integrated_ID_sim))
names(Integrated_ID_sec) <- gsub("1_3", "1", names(Integrated_ID_sec))


df_Idents_Dmel_raw <- as.data.frame(mel_ID) %>%
  tibble::rownames_to_column('ID')  

df_Idents_Dmel_int <- as.data.frame(Integrated_ID_mel) %>%
  tibble::rownames_to_column('ID') %>%
  dplyr::left_join(., df_Idents_Dmel_raw, by="ID")

df_Idents_Dmel_int_summary <- df_Idents_Dmel_int %>%
  dplyr::group_by(mel_ID) %>%
  dplyr::mutate(total_n_cluster=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(mel_ID, Integrated_ID_mel) %>%
  dplyr::summarise(freq=n()/total_n_cluster) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

df_Idents_Dmel_int_summary %>%
  ggplot(.) +
  geom_col() +
  aes(x=mel_ID, y=freq, fill=factor(Integrated_ID_mel, levels=CellType_order$CellType)) +
  theme_bw() +
  scale_fill_manual(values=col_pal_55)

df_Idents_Dmel_int_summary %>%
  ggplot(.) +
  geom_tile(color='black') +
  aes(x=mel_ID, y=factor(Integrated_ID_mel, levels=rev(CellType_order$CellType)),fill=freq) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y= element_text(color='black', face='italic', size=9),
        axis.text.x = element_text(color='black', size=9, angle=90, hjust = 1, vjust = 0.5),
        legend.title = element_text(size=9.5, color='black'),
        legend.text = element_text(size=9, color = 'black'),
        axis.title=element_blank(),
        legend.position='bottom',
        strip.text.x = element_text(margin = margin(0.02,0,0.02,0, "in"))) +
  scale_fill_gradient(low="white", high="blue") +
  labs(x="Dmel cluster", y="Integrated cluster")

df_Idents_Dmel_int_summary %>%
  dplyr::filter(Integrated_ID_mel %in% anno_cluster) %>%
  ggplot(.) +
  geom_tile(color='black') +
  aes(x=mel_ID, y=factor(Integrated_ID_mel, levels=rev(CellType_order$CellType)),fill=freq) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y= element_text(color='black', face='italic', size=9),
        axis.text.x = element_text(color='black', size=9, angle=90, hjust = 1, vjust = 0.5),
        legend.title = element_text(size=9.5, color='black'),
        legend.text = element_text(size=9, color = 'black'),
        axis.title=element_blank(),
        legend.position='bottom',
        strip.text.x = element_text(margin = margin(0.02,0,0.02,0, "in"))) +
  scale_fill_gradient(low="white", high="blue") +
  labs(x="Dmel cluster", y="Integrated cluster")

df_Idents_Dmel_int_summary %>%
  dplyr::filter(Integrated_ID_mel %in% anno_cluster, freq > 0.9) %>%
  ggplot(.) +
  geom_tile(color='black') +
  aes(x=mel_ID, y=factor(Integrated_ID_mel, levels=rev(CellType_order$CellType)),fill=freq) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y= element_text(color='black', face='italic', size=9),
        axis.text.x = element_text(color='black', size=9, angle=90, hjust = 1, vjust = 0.5),
        legend.title = element_text(size=9.5, color='black'),
        legend.text = element_text(size=9, color = 'black'),
        axis.title=element_blank(),
        legend.position='bottom',
        strip.text.x = element_text(margin = margin(0.02,0,0.02,0, "in"))) +
  scale_fill_gradient(low="red", high="blue") +
  labs(x="Dmel cluster", y="Integrated cluster")



df_Idents_Dmel_int_summary_maxfreq <- df_Idents_Dmel_int_summary %>%
  dplyr::group_by(mel_ID) %>%
  dplyr::mutate(max_freq=max(freq)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(freq==max_freq)

df_Idents_Dmel_int_summary_maxfreq %>%
  ggplot(.) +
  geom_histogram(binwidth=0.01) +
  aes(x=max_freq) +
  theme_bw()

df_Idents_Dmel_int_summary_maxfreq %>%
  ggplot(.) +
  geom_col() +
  aes(x=mel_ID, y=max_freq, fill=factor(Integrated_ID_mel, levels=CellType_order$CellType)) +
  theme_bw() +
  scale_fill_manual(values=col_pal_55)
