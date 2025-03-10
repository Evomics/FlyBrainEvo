library(tidyverse)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/../.."))

load("Processed_Data/cluster_specific_DEGs_ref_self.RData")
load("Processed_Data/cluster_specific_DEGs_refmel.RData")

df_deg_clusters_ref_merge1 <- df_deg_clusters_ref_mel %>%
  dplyr::mutate(ref="MEL")

df_deg_clusters_ref_merge2 <- df_deg_clusters_ref_self %>%
  dplyr::mutate(ref="OWN")

df_deg_clusters_ref_merge <- rbind(df_deg_clusters_ref_merge1, df_deg_clusters_ref_merge2)

df_per_expressed_ref <- df_deg_clusters_ref_merge %>%
  dplyr::distinct(gene, cluster, pct.1, pct.2, pair, ref) %>%
  dplyr::mutate(Dmel_refMEL = ifelse(pair =="melsec" & ref=="MEL", pct.1, 
                                     ifelse(pair =="melsim" & ref=="MEL", pct.1, NA)),
                Dsim_refMEL = ifelse(pair =="melsim" & ref=="MEL", pct.2, 
                                     ifelse(pair =="simsec" & ref=="MEL", pct.1, NA)),
                Dsec_refMEL = ifelse(pair =="melsec" & ref=="MEL", pct.2, 
                                     ifelse(pair =="simsec" & ref=="MEL", pct.2, NA)),
                Dmel_refOWN = ifelse(pair =="melsec" & ref=="OWN", pct.1, 
                                     ifelse(pair =="melsim" & ref=="OWN", pct.1, NA)),
                Dsim_refOWN = ifelse(pair =="melsim" & ref=="OWN", pct.2, 
                                     ifelse(pair =="simsec" & ref=="OWN", pct.1, NA)),
                Dsec_refOWN = ifelse(pair =="melsec" & ref=="OWN", pct.2, 
                                     ifelse(pair =="simsec" & ref=="OWN", pct.2, NA)),
                
  ) %>%
  dplyr::distinct(gene, cluster, Dmel_refMEL, Dsim_refMEL, Dsec_refMEL, Dmel_refOWN, Dsim_refOWN, Dsec_refOWN) %>%
  tidyr::gather(-gene, -cluster, key="species_ref", value="pct_exp") %>%
  tidyr::separate(species_ref, into = c('species','ref'), sep='_ref') %>%
  na.omit() %>%
  dplyr::distinct()

save(df_per_expressed_ref, file="Processed_Data/df_per_expressed_ref.RData")


