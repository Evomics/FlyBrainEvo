library(tidyverse)
library(Seurat)
library(ggrepel)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_sum_exp_SCT.RData")
load("Processed_Data/DEG_cluster_ref_sigs_genelist.RData")
load("Processed_Data/shared_genes.RData")
load("Processed_Data/celltype_order.RData")
load("Processed_Data/df_cc_cor_gather_pct_join.RData")
load("Processed_Data/df_per_ref_global_filtered.RData")


TrioBrain.integrated_slim_labeled <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled.rds")

## replicates, pseudobulk ##

### dotplots - across replicates ###

sum_exp_reps <- AverageExpression(TrioBrain.integrated_slim_labeled, features = shared_genes, slot = 'data', add.ident = 'orig.ident')

df_sum_exp_reps_SCT <- as.data.frame(sum_exp_reps$SCT)  %>%
  tibble::rownames_to_column('gene') %>%
  tidyr::gather(-gene, key='cluster', value='expression') %>%
  dplyr::mutate(species=ifelse(grepl("Dmel", cluster)==T, "Dmel",
                               ifelse(grepl("Dsim", cluster)==T, "Dsim",
                                      ifelse(grepl("DsecNoni", cluster)==T, "DsecNoni","Dsec")))) %>%
  dplyr::mutate(reps = ifelse(grepl("rep1", cluster)==T, 1,
                              ifelse(grepl("rep2", cluster)==T, 2,
                                     ifelse(grepl("rep3", cluster)==T, 3,
                                            ifelse(grepl("rep4", cluster)==T, 4,
                                                   ifelse(grepl("rep5", cluster)==T, 5,6))))))

df_sum_exp_reps_SCT$cluster <- gsub("_Dmel_rep.", "", df_sum_exp_reps_SCT$cluster)
df_sum_exp_reps_SCT$cluster <- gsub("_Dsim_rep.", "", df_sum_exp_reps_SCT$cluster)
df_sum_exp_reps_SCT$cluster <- gsub("_DsecNoni_rep.", "", df_sum_exp_reps_SCT$cluster)
df_sum_exp_reps_SCT$cluster <- gsub("_Dsec_rep.", "", df_sum_exp_reps_SCT$cluster)

unique(df_sum_exp_reps_SCT$cluster)

df_sum_exp_reps_SCT$species <- factor(df_sum_exp_reps_SCT$species, levels = c("Dmel", "Dsim", "Dsec", "DsecNoni"))

#save(anno_cluster, df_sum_exp_reps_SCT, file="Processed_Data/df_sum_exp_reps_SCT.RData")
load("Processed_Data/df_sum_exp_reps_SCT.RData")

## wirins ##

Dmel_gene_name <- read.table(file="Dmel/Dmel_ID_genename.tsv") %>% dplyr::distinct(V1,V2) %>%
  dplyr::rename(FBID=V1, gene=V2)

df_wirin <- data.table::fread("genelist/list_wirins_FBgg0001381.txt", header = F) %>%
  dplyr::mutate(class="Wirin") %>%
  dplyr::select(FBID=V1, class) %>%
  dplyr::left_join(., Dmel_gene_name, by="FBID") %>% ### few genes without matching gene_name
  na.omit()

wirin_genes <- df_wirin$gene

wirin_genes_clean <- intersect(wirin_genes, clean_gene_set) ### clean gene list by filtering genome effects ###

## correlation ##

df_cc_cor_gather_pct_join %>%
  dplyr::filter(pair=="melsec", pct_exp > 0.1, gene %in% wirin_genes_clean, value=='r') %>%
  ggplot(.) +
  geom_smooth(span=2, size=0.5) +
  geom_point(alpha=0.5, size =1.5) +
  geom_text_repel(aes(label=gene), fontface = "italic", size=3) +
  aes(x=pct_exp*100, y=cor) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=12, color='black'),
        axis.text= element_text(size=11, color='black'),
        strip.text = element_text(size=11, color='black'),
        legend.position = 'none') +
  facet_grid(~factor(pair, levels=c('melsec','melsim','simsec'), labels = c("Dmel vs Dsec", "Dmel vs Dsim", "Dsim vs Dsec"))) +
  scale_x_continuous(expand=c(0.01,0.01)) +
  #scale_y_continuous(expand=c(0.001,0.001)) +
  labs(x="Percent expressed (%)", y="Pearson's r")

ggsave(file="Plots/Manuscript/Fig6_wirin_melsec_r.pdf", width=4, height= 4)

## dpr6 ##

gen='dpr6'

df_wirin_join <- df_sum_exp_SCT %>%
  dplyr::filter(gene %in% wirin_genes_clean, species != "DsecNoni")

df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
  dplyr::filter(gene==gen) %>%
  tidyr::spread(key='species', value='expression')

## heatplot ##

df_sum_exp_SCT %>%
  dplyr::filter(gene ==gen, species != "DsecNoni", cluster %in% anno_cluster) %>%
  ggplot(.) +
  geom_tile(color='black') +
  aes(x=factor(cluster, levels=anno_cluster), y=factor(species, levels = c("Dsec", "Dsim", "Dmel")), fill=expression) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y= element_text(color='black', face='italic', size=9),
        axis.text.x = element_text(color='black', size=9, angle=45, hjust = 1, vjust = 1),
        legend.title = element_text(size=9.5, color='black'),
        legend.text = element_text(size=9, color = 'black'),
        axis.title=element_blank(),
        legend.position='bottom',
        strip.text.x = element_text(margin = margin(0.02,0,0.02,0, "in"))) +
  scale_fill_gradient(low="white", high="navy") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0), position = "right")

ggsave(file="Plots/Manuscript/Fig6_heatmap_dpr6.pdf", width=10, height= 2)

## sample dotplots ##

select_cluster <- c("Tbh", "SIFa", "fru(Glu)",  "Da6(I)","α'β'-KC",
                    "RYa","Da6(E)","fru(Ms)",
                    "αβ-KC","γ-KC")

df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster %in% select_cluster, gene ==gen, species != "DsecNoni") %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.4, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=8),
        axis.text.y= element_text(color='black', size=8),
        axis.title.y= element_text(color='black', size=9),
        axis.title.x=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=8),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=8),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8"))+
  facet_wrap(~factor(cluster, levels=select_cluster), scales = 'free', nrow=2) 

ggsave(glue::glue("Plots/Manuscript/Fig6_dotplot_{gen}.pdf"), width=6, height =3)

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_smooth(size= 0.3, method="glm") +
  geom_point(alpha = 0.6, size=1.2) +
  geom_text_repel(data=dplyr::filter(df_sum_exp_SCT_select, 
                  #                   cluster %in% anno_cluster),
                                      cluster %in% select_cluster),
                  aes(label=cluster), label.size=0.1, segment.size = 0.1, box.padding = 0.1,nudge_x = 0.1,
                  label.padding = 0.1, size = 3) +
  aes(x=Dmel, y=Dsec) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title= element_text(color='black', size=11, face='italic')) +
  ggtitle(paste(gen, "| r =", round(cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method='pearson'),3)))

## sr - acitity regulated gene ##

gen='sr'

gen %in% clean_gene_set

df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
  dplyr::filter(gene==gen) %>%
  tidyr::spread(key='species', value='expression')

## heatplot ##

df_sum_exp_SCT %>%
  dplyr::filter(gene ==gen, cluster %in% anno_cluster, cluster != "Blood") %>%
  ggplot(.) +
  geom_tile(color='black') +
  aes(x=factor(cluster, levels=anno_cluster), y=factor(species, levels = c("DsecNoni", "Dsec", "Dsim", "Dmel")), fill=expression) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y= element_text(color='black', face='italic', size=9),
        axis.text.x = element_text(color='black', size=9, angle=45, hjust = 1, vjust = 1),
        legend.title = element_text(size=9.5, color='black'),
        legend.text = element_text(size=9, color = 'black'),
        axis.title=element_blank(),
        legend.position='bottom',
        strip.text.x = element_text(margin = margin(0.02,0,0.02,0, "in"))) +
  scale_fill_gradient(low="white", high="navy") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0), position = "right")

ggsave(file="Plots/Manuscript/Fig6_heatmap_dpr6.pdf", width=10, height= 2)

## sample dotplots ##

select_cluster <- c("Tbh", "SIFa", "fru(Glu)",  "Da6(I)","α'β'-KC",
                    "RYa","Da6(E)","fru(Ms)",
                    "αβ-KC","γ-KC")

df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster %in% select_cluster, gene ==gen, species != "DsecNoni") %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 1.4, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=8),
        axis.text.y= element_text(color='black', size=8),
        axis.title.y= element_text(color='black', size=9),
        axis.title.x=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=8),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=8),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8"))+
  facet_wrap(~factor(cluster, levels=select_cluster), scales = 'free', nrow=2) 

ggsave(glue::glue("Plots/Manuscript/Fig7_dotplot_{gen}.pdf"), width=6, height =3)

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_smooth(size= 0.3, method="glm") +
  geom_point(alpha = 0.6, size=1.2) +
  geom_text_repel(data=dplyr::filter(df_sum_exp_SCT_select, 
                                     #                   cluster %in% anno_cluster),
                                     cluster %in% select_cluster),
                  aes(label=cluster), label.size=0.1, segment.size = 0.1, box.padding = 0.1,nudge_x = 0.1,
                  label.padding = 0.1, size = 3) +
  aes(x=Dmel, y=Dsec) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title= element_text(color='black', size=11, face='italic')) +
  ggtitle(paste(gen, "| r =", round(cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method='pearson'),3)))







 #### chemoconnectome ###

cc_genes <- c("VAChT", "ChAT", "mAChR-A", "mAChR-B", "mAChR-C", ## Ach signaling
              "nAChRalpha1", "nAChRalpha2", "nAChRalpha3","nAChRalpha5",  "nAChRalpha6", "nAChRalpha7", 
              "nAChRbeta1","nAChRbeta2", "nAChRbeta3", 
              "VGlut", "mGluR", "GluRIA", "GluRIB", "GluRIIE", "GluClalpha", ## Glu signaling
              "Gad1", "VGAT","Rdl", "Lcch3", "Grd", "CG8916", "GABA-B-R1","GABA-B-R2","GABA-B-R3", ## GABA
              "Vmat", "ple", "DAT", "Dop1R1","Dop1R2", "Dop2R", "DopEcR",  ## MON-Dop signaling
              "Tdc2", "Tbh","Oct-TyrR", "TyrR", "TyrRII", "Octalpha2R", "Octbeta1R",  "Octbeta2R",  "Octbeta3R", "Oamb", ## Tyr-Oct signaling
              "Trh", "SerT",   "5-HT1A","5-HT1B",  "5-HT2A", "5-HT2B", "5-HT7",  ## Ser signaling
              "Hdc",  "ort", "HisCl1")

df_cc_cor_gather_pct_join %>%
  dplyr::filter(pair=="melsec", pct_exp > 0.1, value=='rho', gene %in% cc_genes) %>%
  ggplot(.) +
  geom_point(alpha=0.5, size =1.5) +
  geom_smooth(span=2) +
  geom_text_repel(aes(label=gene)) +
  aes(x=pct_exp*100, y=cor) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=12, color='black'),
        axis.text= element_text(size=11, color='black'),
        strip.text = element_text(size=11.5, color='black'),
        legend.position = 'none') +
  facet_grid(~factor(pair, levels=c('melsec','melsim','simsec'), labels = c("Dmel vs Dsec", "Dmel vs Dsim", "Dsim vs Dsec"))) +
  #scale_x_continuous(expand=c(0.001,0.001)) +
  #scale_y_continuous(expand=c(0.001,0.001)) +
  labs(x="Percent expressed (%)", y="Spearman's rho")

## GABA-B-R1 ##

gen='GABA-B-R1'

df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
  dplyr::filter(gene==gen) %>%
  tidyr::spread(key='species', value='expression')

select_cluster <- c("Da6(I)", "OPN", "Hug", "ort", "Da5/7", "sNPF(I)", "Tk", "Crz", "Proc", "RYa")

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_smooth(size= 0.3, method="glm") +
  geom_point(alpha = 0.6, size=1.2) +
  geom_text_repel(data=dplyr::filter(df_sum_exp_SCT_select, 
                                      cluster %in% anno_cluster),
#                                    cluster %in% select_cluster),
             aes(label=cluster), label.size=0.1, segment.size = 0.1, box.padding = 0.1,nudge_x = 0.1,
            label.padding = 0.1, size = 3) +
  aes(x=Dmel, y=Dsec) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title= element_text(color='black', size=11, face='italic')) +
  ggtitle(paste(gen, "| rho =", round(cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method='spearman'),3)))

ggsave(glue::glue("Plots/Manuscript/Fig5_{gen}_allcluster_scatter_simsec.pdf"), width=3.5, height =3.5)


df_sum_exp_reps_SCT %>%
  dplyr::filter(cluster %in% select_cluster, gene ==gen, species != "DsecNoni") %>%
  ggplot(.) +
  #geom_point() +
  geom_dotplot(dotsize = 0.9, binaxis='y', stackdir='center') +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="darkorange", alpha =0.8, size = 0.2) +
  aes(x=species, y=expression, fill=species) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x= element_text(color='black', face='italic', size=8),
        axis.text.y= element_text(color='black', size=8),
        axis.title.y= element_text(color='black', size=9),
        axis.title.x=element_blank(), 
        strip.text.x = element_text(margin = margin(0.04,0,0.04,0, "in"), size=8),
        strip.text.y = element_text(margin = margin(0,0.04,0,0.04, "in"), face='italic', size=8),
        legend.position='none') +
  scale_fill_manual(values= c("#E41A1C", "#4DAF4A",  "#377EB8"))+
  facet_wrap(~factor(cluster, levels=select_cluster), scales = 'free', nrow=1) 

ggsave(file="Plots/Manuscript/Fig4_dotplot_5HT2A.pdf", width=9, height= 1.5)


