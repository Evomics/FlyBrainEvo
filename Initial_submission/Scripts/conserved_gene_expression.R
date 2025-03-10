library(Seurat)
library(tidyverse)
library(ggrepel)
library(cowplot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/df_per_expressed_ref.RData")
load("Processed_Data/df_sum_exp_SCT.RData")
load("Processed_Data/global_DEGs.RData")
load("Processed_Data/df_deg_global_ref_merge.RData")
load("Processed_Data/df_deg_sigs.RData")


### number of expressed clusters ###

threshold = 0.1

df_exp_clusters <- df_per_expressed_ref %>%
  dplyr::filter(pct_exp >= threshold)

df_exp_clusters_summary <- df_exp_clusters %>%
  dplyr::group_by(gene, species, ref) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup()

df_exp_clusters_summary %>%
  ggplot(.) +
  geom_histogram() +
  aes(x=n) +
  theme_bw() +
  facet_grid(ref~species)

### clean gene list by filtering genome effects - global ###

df_per_expressed_ref_global <- df_deg_global_ref_merge %>%
  dplyr::distinct(gene, pct.1, pct.2, pair, ref) %>%
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
  dplyr::distinct(gene, Dmel_refMEL, Dsim_refMEL, Dsec_refMEL, Dmel_refOWN, Dsim_refOWN, Dsec_refOWN) %>%
  tidyr::gather(-gene, key="species_ref", value="pct_exp") %>%
  tidyr::separate(species_ref, into = c('species','ref'), sep='_ref') %>%
  na.omit() %>%
  dplyr::distinct()

df_per_expressed_ref_global_diff <- df_per_expressed_ref_global %>%
  tidyr::spread(key="ref", value='pct_exp') %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(RefMEL_rank_exp=min_rank(-MEL), RefOWN_rank_exp=min_rank(-OWN)) %>%
  dplyr::mutate(rank_diff = abs(RefMEL_rank_exp-RefOWN_rank_exp), rank_diff_ratio=rank_diff/n())

df_per_expressed_ref_global_diff_filter <- df_per_expressed_ref_global_diff %>%
  dplyr::filter(rank_diff_ratio < 0.05) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(n_species = n_distinct(species)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_species ==3)

length(unique(df_per_expressed_ref_global_diff_filter$gene))

clean_gene_set <- unique(df_per_expressed_ref_global_diff_filter$gene)

#save(clean_gene_set, df_per_expressed_ref_global_diff_filter, file="Processed_Data/df_per_ref_global_filtered.RData")
load("Processed_Data/df_per_ref_global_filtered.RData")

### correlation ###

## scatter plots ##

gene_list <- intersect(unique(df_per_expressed_ref_global_diff_filter$gene), intersect(unique(df_exp_clusters$gene), unique(df_sum_exp_SCT$gene)))

df_cc_cor_rho <- NULL
df_cc_cor_r <- NULL

for (gen in gene_list) {
  
  df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
    dplyr::filter(gene==gen) %>%
    tidyr::spread(key='species', value='expression')
  
  cc_cor_melsec_rho <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method='spearman')
  cc_cor_melsim_rho <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsim, method='spearman')
  cc_cor_simsec_rho <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$Dsim, method='spearman')
  cc_cor_secNoni_rho <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$DsecNoni, method='spearman')
  
  cc_cor_melsec_r <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method='pearson')
  cc_cor_melsim_r <- cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsim, method='pearson')
  cc_cor_simsec_r <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$Dsim, method='pearson')
  cc_cor_secNoni_r <- cor(df_sum_exp_SCT_select$Dsec, df_sum_exp_SCT_select$DsecNoni, method='pearson')
  
  df_cc_cor_r_temp <- data.frame(gene=gen, 
                               melsec=cc_cor_melsec_r,
                               melsim=cc_cor_melsim_r, 
                               simsec=cc_cor_simsec_r, 
                               secNoni=cc_cor_secNoni_r,
                               value='r')
  
  df_cc_cor_r <- rbind(df_cc_cor_r, df_cc_cor_r_temp)
  
  df_cc_cor_rho_temp <- data.frame(gene=gen, 
                               melsec=cc_cor_melsec_rho,
                               melsim=cc_cor_melsim_rho, 
                               simsec=cc_cor_simsec_rho, 
                               secNoni=cc_cor_secNoni_rho,
                               value='rho')
  
  df_cc_cor_rho <- rbind(df_cc_cor_rho, df_cc_cor_rho_temp)


}

df_cc_cor_rho_gather <- df_cc_cor_rho %>%
  tidyr::gather(-gene, -value, key='pair', value='cor')

df_cc_cor_r_gather <- df_cc_cor_r %>%
  tidyr::gather(-gene,-value, key='pair', value='cor')

df_cc_cor_gather <- rbind(df_cc_cor_rho_gather, df_cc_cor_r_gather)

df_cc_cor <- df_cc_cor_gather %>%
  tidyr::spread(key='value', value='cor')


#save(df_cc_cor_rho, df_cc_cor_r, df_cc_cor_gather,  df_cc_cor, file="Processed_Data/exp_gene_cor.RData")
load("Processed_Data/exp_gene_cor.RData")

df_cc_cor %>%
  dplyr::filter(pair == 'melsec') %>%
  ggplot(.) +
  geom_point(size=0.5, alpha=0.6) +
  aes(x=r, y=rho) +
  theme_bw()

## Examples for correlated/uncorrelated expression pattern ##

gen='bru3'
gen='unc-13'
gen='gw'
gen='MED26'

df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
  dplyr::filter(gene==gen) %>%
  tidyr::spread(key='species', value='expression')

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_smooth(size= 0.3, method="glm") +
  geom_point(alpha = 0.6, size=1.2) +
 # geom_text_repel(data=dplyr::filter(df_sum_exp_SCT_select, 
  #                                    cluster %in% c("ENS", "AST","PRN", "SUB", "Blood", "CTX", "Fat",
   #                                                  "FMRFa(1)","fru(Ms)", "Clock","Mip(2)", "IPC","sNPF(E)", "Hug", "CCHa2(2)", "Dh44", "Crz",
    #                                                 "fru(E)", "fru(Glu)", "SIFa", "FMRFa(2)", 
     #                                                "Proc", "OCTY", "DOP", "Poxn", "Da5/7", "Da6(E)",
      #                                               "αβ-KC", "γ-KC" ,    "α'β'-KC")),
       #           aes(label=cluster), label.size=0.1, segment.size = 0.1, box.padding = 0.1,nudge_y = 0.5,
        #          label.padding = 0.1, size = 3) +
  aes(x=Dmel, y=Dsec) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title= element_text(color='black', size=11, face='italic')) +
  ggtitle(paste(gen, "| rho =", round(cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method='spearman'),3)))

ggsave(glue::glue("Plots/Manuscript/Fig5_{gen}_allcluster_scatter_simsec.pdf"), width=3.5, height =3.5)


plot_bru3 <- df_sum_exp_SCT %>%
  dplyr::filter(gene=='bru3') %>%
  tidyr::spread(key='species', value='expression') %>%
  ggplot(.) +
  geom_smooth(size= 0.3, method="glm") +
  geom_point(alpha = 0.6, size=1.2) +
  aes(x=Dmel, y=Dsec) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title.y= element_text(color='black', size=11, face='italic'),
        axis.title.x = element_blank()) +
  ggtitle(paste('bru3', "| rho =", round(cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method='spearman'),3)))

plot_unc13 <- df_sum_exp_SCT %>%
  dplyr::filter(gene=='unc-13') %>%
  tidyr::spread(key='species', value='expression') %>%
  ggplot(.) +
  geom_smooth(size= 0.3, method="glm") +
  geom_point(alpha = 0.6, size=1.2) +
  aes(x=Dmel, y=Dsec) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title.y= element_text(color='black', size=11, face='italic'),
        axis.title.x = element_blank()) +
  ggtitle(paste('unc-13', "| rho =", round(cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method='spearman'),3)))

plot_gw <- df_sum_exp_SCT %>%
  dplyr::filter(gene=='gw') %>%
  tidyr::spread(key='species', value='expression') %>%
  ggplot(.) +
  geom_smooth(size= 0.3, method="glm") +
  geom_point(alpha = 0.6, size=1.2) +
  aes(x=Dmel, y=Dsec) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title= element_text(color='black', size=11, face='italic')) +
  ggtitle(paste('gw', "| rho =", round(cor(df_sum_exp_SCT_select$Dmel, df_sum_exp_SCT_select$Dsec, method='spearman'),3)))
  
plot_rho_gene_examples <- cowplot::plot_grid(plot_bru3,plot_unc13, plot_gw, ncol=1, align='v')

plot_rho_gene_examples

ggsave("Plots/Manuscript/Fig5_exp_cor_examples.pdf", width=2.8, height = 8)


 ### pct exp x correlation ###

df_per_select_Dmel <- df_per_expressed_ref_global_diff_filter %>%
  dplyr::filter(OWN > 0.05, species == 'Dmel') %>%
  dplyr::select(gene, pct_exp=OWN) ## 1686 genes are selected for the analysis

cor_gene_set <- unique(df_per_select_Dmel$gene)

## specificity of these 1686 genes ##

df_per_expressed_ref_cor <- df_per_expressed_ref %>%
  dplyr::filter(gene %in% cor_gene_set) %>%
  dplyr::group_by(gene, cluster, species) %>%
  dplyr::summarise(pct_exp_max=max(pct_exp)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(pct_exp_max >= 0.05) %>%
  dplyr::group_by(species, gene) %>%
  dplyr::mutate(n_celltype=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(species == "Dmel") %>%
  dplyr::distinct(gene, n_celltype)

mean(df_per_expressed_ref_cor$n_celltype)

df_deg_clusters_sigs_select <- df_deg_clusters_ref_merge_sigs %>%
  dplyr::distinct(gene, pair) %>%
  dplyr::mutate(DEG=T) 

df_cc_cor_gather_pct_join <- df_cc_cor_gather %>%
  dplyr::left_join(., df_per_select_Dmel, by='gene') %>%
  na.omit() %>%
  dplyr::left_join(., df_deg_clusters_sigs_select, by=c('gene', 'pair')) %>%
  dplyr::mutate(DEG=ifelse(is.na(DEG), F, T))

save(df_cc_cor_gather_pct_join, file="Processed_Data/df_cc_cor_gather_pct_join.RData")
load("Processed_Data/df_cc_cor_gather_pct_join.RData")

df_cc_cor_gather_pct_join %>%
  dplyr::filter(pair=="melsec") %>%
  ggplot(.) +
  geom_point(alpha=0.8, size=1, shape=21) +
  #geom_smooth() +
  aes(x=pct_exp, y=cor, fill=DEG) +
  theme_bw() +
  facet_grid(~value)



### density plot ###

df_cc_cor_gather_pct_join %>%
  dplyr::filter(pair!="secNoni", pct_exp > 0.1) %>%
  ggplot(.) +
  geom_density_2d_filled(alpha=0.9) +
  geom_point(alpha=0.5, size =0.1) +
  aes(x=pct_exp, y=cor) +
  theme_bw() +
  theme() +
  facet_grid(value~pair) +
  scale_x_continuous(expand=c(0.01,0.01)) +
  scale_y_continuous(expand=c(0.01,0.01))

plot_2d_density_rho <- df_cc_cor_gather_pct_join %>%
  dplyr::filter(pair!="secNoni", pct_exp > 0.1, value=='rho') %>%
  ggplot(.) +
  geom_density_2d_filled(alpha=0.85) +
  geom_point(alpha=0.5, size =0.01) +
  geom_hline(yintercept=0.7, linetype=2, color='red', alpha=0.8, size=0.5) +
  aes(x=pct_exp*100, y=cor) +
  theme_bw() +
  theme(axis.title.y = element_text(size=12, color='black'),
        axis.text.y = element_text(size=11, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size=11.5, color='black'),
        legend.position = 'none') +
  facet_grid(~factor(pair, levels=c('melsec','melsim','simsec'), labels = c("Dmel vs Dsec", "Dmel vs Dsim", "Dsim vs Dsec"))) +
  scale_x_continuous(expand=c(0.001,0.001)) +
  scale_y_continuous(expand=c(0.001,0.001)) +
  labs(x="Percent expressed (%)", y="Spearman's rho")

plot_histo_rho_pct <- df_cc_cor_gather_pct_join %>%
  dplyr::filter(pair!="secNoni", cor > 0.7, pct_exp > 0.1, value=='rho') %>%
  ggplot(.) +
  geom_histogram(binwidth=3, boundary=10) +
  aes(x=pct_exp*100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=12, color='black'),
        axis.text = element_text(size=11, color='black'),
        strip.text.x = element_blank()) +
  facet_grid(~pair) +
  scale_x_continuous(expand=c(0.01,0.01)) +
  labs(x="Percent expressed (%)", y="Gene counts (rho > 0.7)")

plot_rho_pct_exp <- cowplot::plot_grid(plot_2d_density_rho, plot_histo_rho_pct, ncol=1, align='v', rel_heights = c(2,1))

plot_rho_pct_exp

sggsave("Plots/Manuscript/Fig5_2D_density_rho_pct.pdf", width=10, height = 5)

df_cc_cor_gather_pct_join_strong_rho <- df_cc_cor_gather_pct_join %>%
  dplyr::filter(value=='rho', pair != "secNoni") %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(min(cor) > 0.7) %>%
  dplyr::ungroup()

length(unique(df_cc_cor_gather_pct_join_strong_rho$gene))

conser_genes <- unique(df_cc_cor_gather_pct_join_strong_rho$gene)
write.csv(conser_genes, file="conserved_genes.csv")
conser_genes_control <- unique(df_cc_cor_gather_pct_join$gene)
write.csv(conser_genes_control, file="conserved_genes_control.csv")

## FlyEnrichR ###

FER_cons_BP <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved368/GO_Biological_Process_2018_table.txt") %>%
  dplyr::mutate(class="BP", group='conserved', n_gene=368) %>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_cons_MF <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved368/GO_Molecular_Function_2018_table.txt") %>%
  dplyr::mutate(class="MF", group='conserved', n_gene=368)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_cons_CC <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved368/GO_Cellular_Component_2018_table.txt") %>%
  dplyr::mutate(class="CC", group='conserved', n_gene=368)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')

FER_control_BP <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved_control1686/GO_Biological_Process_2018_table.txt") %>%
  dplyr::mutate(class="BP", group='control', n_gene=1686)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_control_MF <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved_control1686/GO_Molecular_Function_2018_table.txt") %>%
  dplyr::mutate(class="MF", group='control', n_gene=1686)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_control_CC <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved_control1686/GO_Cellular_Component_2018_table.txt") %>%
  dplyr::mutate(class="CC", group='control', n_gene=1686)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')

df_FER_cons_int <- rbind(FER_cons_BP, FER_cons_MF, FER_cons_CC, FER_control_BP, FER_control_MF, FER_control_CC) %>%
  tidyr::separate(Overlap, c('n_overlap','n_gene_class')) %>%
  dplyr::filter(as.numeric(n_overlap) >=5) %>%
  dplyr::mutate(frac_to_group=as.numeric(n_overlap)/n_gene, frac_to_class=as.numeric(n_overlap)/as.numeric(n_gene_class)) %>%
  dplyr::mutate(frac_to_group_norm_class=as.numeric(n_overlap)/n_gene/as.numeric(n_gene_class)) %>%
  dplyr::group_by(Term) %>%
  dplyr::mutate(n_Term = n(),
                n_overlap_conserved=min(as.numeric(n_overlap)), n_overlap_control=max(as.numeric(n_overlap)), 
                fraction=n_overlap_conserved/n_overlap_control) %>%
  dplyr::mutate(fold_enrich = ifelse(n_Term==2, 1686/368*n_overlap_conserved/n_overlap_control, NA)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Genes) %>%
  dplyr::filter(p_adj==min(p_adj) & Z_score==min(Z_score)) %>%
  dplyr::ungroup()
  
## add fold enrichment

df_FER_cons_int_top10 <- df_FER_cons_int %>%
  dplyr::filter(group == "conserved") %>%
  dplyr::group_by(class) %>%
  dplyr::top_n(10, fold_enrich)

top10_term <- unique(df_FER_cons_int_top10$Term)

df_FER_cons_int %>%
  dplyr::filter(Term %in% top10_term, group == "conserved", fold_enrich > 1) %>%
  ggplot(.) +
  geom_point(shape=21) +
  aes(x=fold_enrich, y=reorder(Term, fold_enrich), fill=-log10(p_adj)) +
  theme_bw() +
  facet_wrap(~class, nrow=1, scales='free')

plot_FER_cons_BP <- df_FER_cons_int %>%
  dplyr::filter(Term %in% top10_term, group == "conserved", fold_enrich > 1, class=="BP") %>%
  ggplot(.) +
  geom_point(shape=21, size=3) +
  aes(x=fold_enrich, y=reorder(Term, fold_enrich), fill=-log10(p_adj)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title.x= element_text(color='black', size=11),
        legend.key.height = unit(0.1,"in"),
        legend.position=c(0.25,0.85),
        legend.direction='horizontal',
        legend.title = element_text(color='black', size=9),
        legend.text = element_text(color='black', size=8.5)) +
  facet_wrap(~class, ncol=1, scales='free') +
  labs(fill='-log10(adjusted P-value)', x='Fold Enrichment')

plot_FER_cons_BP

plot_FER_cons_MF <- df_FER_cons_int %>%
  dplyr::filter(Term %in% top10_term, group == "conserved", fold_enrich > 1, class=="MF") %>%
  ggplot(.) +
  geom_point(shape=21, size=3) +
  aes(x=fold_enrich, y=reorder(Term, fold_enrich), fill=-log10(p_adj)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title.x= element_text(color='black', size=11),
        legend.key.height = unit(0.1,"in"),
        legend.position=c(0.25,0.85),
        legend.direction='horizontal',
        legend.title = element_text(color='black', size=9),
        legend.text = element_text(color='black', size=8.5)) +
  facet_wrap(~class, ncol=1, scales='free') +
  labs(fill='-log10(adjusted P-value)', x='Fold Enrichment')

plot_FER_cons_MF

plot_FER_cons_CC <- df_FER_cons_int %>%
  dplyr::filter(Term %in% top10_term, group == "conserved", fold_enrich > 1, class=="CC") %>%
  ggplot(.) +
  geom_point(shape=21, size=3) +
  aes(x=fold_enrich, y=reorder(Term, fold_enrich), fill=-log10(p_adj)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title.x= element_text(color='black', size=11),
        legend.key.height = unit(0.1,"in"),
        legend.position=c(0.25,0.7),
        legend.direction='horizontal',
        legend.title = element_text(color='black', size=9),
        legend.text = element_text(color='black', size=8.5)) +
  facet_wrap(~class, ncol=1, scales='free') +
  labs(fill='-log10(adjusted P-value)', x='Fold Enrichment')

plot_FER_cons_CC

plot_FER_cons_GO <- cowplot::plot_grid(plot_FER_cons_BP, plot_FER_cons_MF, plot_FER_cons_CC, ncol=1, align='v', rel_heights = c(1.04,1,0.67))

plot_FER_cons_GO

ggsave(plot_FER_cons_GO, file="Plots/Manuscript/working/230403_Extended_Data_Fig2c_ver4.pdf", width=12, height = 6.3)

#### TFs ###

Dmel_gene_name <- read.table(file="Dmel/Dmel_ID_genename.tsv") %>% dplyr::distinct(V1,V2) %>%
  dplyr::rename(FBID=V1, gene=V2)

df_TF <- data.table::fread("genelist/list_TF_GO0000981.txt", header = F) %>%
  dplyr::mutate(class="TF") %>%
  dplyr::select(FBID=V1, class) %>%
  dplyr::left_join(., Dmel_gene_name, by="FBID") %>% ### few genes without matching gene_name
  na.omit()

TF_genes <- df_TF$gene

df_cc_cor_gather_pct_join_strong_rho_TF <- df_cc_cor_gather_pct_join_strong_rho %>%
  dplyr::filter(gene %in% TF_genes) 

TF_conser <- unique(df_cc_cor_gather_pct_join_strong_rho_TF$gene)

##

df_sum_exp_SCT_TF <- df_sum_exp_SCT %>%
  dplyr::filter(gene %in% TF_conser, species != "DsecNoni") 

df_sum_exp_SCT_TF_summary <- df_sum_exp_SCT_TF %>%
  dplyr::filter(expression>0) %>%
  dplyr::group_by(cluster, species) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() 

df_sum_exp_SCT_TF_summary %>%
  dplyr::summarise(mean=mean(n))

### normalized expression level

df_scale_ts <- df_sum_exp_SCT %>%
  dplyr::filter(species == "Dmel", gene %in% TF_conser, cluster %in% anno_cluster) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(max_expression=max(expression)) %>%
  dplyr::ungroup() %>%
  #dplyr::filter(max_expression > 0.3) %>%
  dplyr::mutate(norm_exp = expression/max_expression) 

df_scale_ts %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  ggplot(.) +
  geom_tile(color='black', size=0.1) +
  aes(y=factor(gene, levels=rev(TF_conser)), x=factor(cluster, levels=anno_cluster), fill=norm_exp) +
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

#ggsave(file="Plots/Manuscript/Fig2_heatmap_ts_gene_cluster_exp.pdf", width=10, height= 8)

df_scale_ts$cluster <- factor(df_scale_ts$cluster, levels= anno_cluster)

mt_scale_ts <- df_scale_ts %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  dplyr::select(gene, cluster, norm_exp) %>%
  tidyr::spread(key='gene', value='norm_exp') %>%
  dplyr::select(-cluster) %>%
  as.matrix()

rownames(mt_scale_ts) <- anno_cluster

plot_ts <- pheatmap::pheatmap(mt_scale_ts, angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                              annotation_names_col = T, show_colnames = T, cluster_rows = T,
                              treeheight_col =20, treeheight_row = 20, border_color = "black", na_col = "grey90")



### With pearson's R

df_cc_cor_gather_pct_join_strong_r <- df_cc_cor_gather_pct_join %>%
  dplyr::filter(value=='r', pair != "secNoni") %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(min(cor) > 0.7) %>%
  dplyr::ungroup()

length(unique(df_cc_cor_gather_pct_join_strong_r$gene))

conser_genes_r <- unique(df_cc_cor_gather_pct_join_strong_r$gene)
write.csv(conser_genes_r, file="conserved_genes_r.csv")


## FlyEnrichR ###

FER_cons_BP <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved368/GO_Biological_Process_2018_table.txt") %>%
  dplyr::mutate(class="BP", group='conserved', n_gene=368) %>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_cons_MF <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved368/GO_Molecular_Function_2018_table.txt") %>%
  dplyr::mutate(class="MF", group='conserved', n_gene=368)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_cons_CC <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved368/GO_Cellular_Component_2018_table.txt") %>%
  dplyr::mutate(class="CC", group='conserved', n_gene=368)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')

FER_control_BP <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved_control1686/GO_Biological_Process_2018_table.txt") %>%
  dplyr::mutate(class="BP", group='control', n_gene=1686)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_control_MF <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved_control1686/GO_Molecular_Function_2018_table.txt") %>%
  dplyr::mutate(class="MF", group='control', n_gene=1686)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_control_CC <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved_control1686/GO_Cellular_Component_2018_table.txt") %>%
  dplyr::mutate(class="CC", group='control', n_gene=1686)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')

df_FER_cons_int <- rbind(FER_cons_BP, FER_cons_MF, FER_cons_CC, FER_control_BP, FER_control_MF, FER_control_CC) %>%
  tidyr::separate(Overlap, c('n_overlap','n_gene_class')) %>%
  dplyr::filter(as.numeric(n_overlap) >=5) %>%
  dplyr::mutate(frac_to_group=as.numeric(n_overlap)/n_gene, frac_to_class=as.numeric(n_overlap)/as.numeric(n_gene_class)) %>%
  dplyr::mutate(frac_to_group_norm_class=as.numeric(n_overlap)/n_gene/as.numeric(n_gene_class)) %>%
  dplyr::group_by(Term) %>%
  dplyr::mutate(n_Term = n(),
                n_overlap_conserved=min(as.numeric(n_overlap)), n_overlap_control=max(as.numeric(n_overlap)), 
                fraction=n_overlap_conserved/n_overlap_control) %>%
  dplyr::mutate(fold_enrich = ifelse(n_Term==2, 1686/368*n_overlap_conserved/n_overlap_control, NA)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Genes) %>%
  dplyr::filter(p_adj==min(p_adj) & Z_score==min(Z_score)) %>%
  dplyr::ungroup()

## add fold enrichment

df_FER_cons_int_top10 <- df_FER_cons_int %>%
  dplyr::filter(group == "conserved") %>%
  dplyr::group_by(class) %>%
  dplyr::top_n(10, fold_enrich)

top10_term <- unique(df_FER_cons_int_top10$Term)

df_FER_cons_int %>%
  dplyr::filter(Term %in% top10_term, group == "conserved", fold_enrich > 1) %>%
  ggplot(.) +
  geom_point(shape=21) +
  aes(x=fold_enrich, y=reorder(Term, fold_enrich), fill=-log10(p_adj)) +
  theme_bw() +
  facet_wrap(~class, ncol=1, scales='free')


#### TFs ###

Dmel_gene_name <- read.table(file="Dmel/Dmel_ID_genename.tsv") %>% dplyr::distinct(V1,V2) %>%
  dplyr::rename(FBID=V1, gene=V2)

df_TF <- data.table::fread("genelist/list_TF_GO0000981.txt", header = F) %>%
  dplyr::mutate(class="TF") %>%
  dplyr::select(FBID=V1, class) %>%
  dplyr::left_join(., Dmel_gene_name, by="FBID") %>% ### few genes without matching gene_name
  na.omit()

TF_genes <- df_TF$gene

df_cc_cor_gather_pct_join_strong_r_TF <- df_cc_cor_gather_pct_join_strong_r %>%
  dplyr::filter(gene %in% TF_genes) 

TF_conser <- unique(df_cc_cor_gather_pct_join_strong_rho_TF$gene)

##

df_sum_exp_SCT_TF <- df_sum_exp_SCT %>%
  dplyr::filter(gene %in% TF_conser, species != "DsecNoni") 

df_sum_exp_SCT_TF_summary <- df_sum_exp_SCT_TF %>%
  dplyr::filter(expression>0) %>%
  dplyr::group_by(cluster, species) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() 

df_sum_exp_SCT_TF_summary %>%
  dplyr::summarise(mean=mean(n))

### normalized expression level

df_scale_ts <- df_sum_exp_SCT %>%
  dplyr::filter(species == "Dmel", gene %in% TF_conser, cluster %in% anno_cluster) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(max_expression=max(expression)) %>%
  dplyr::ungroup() %>%
  #dplyr::filter(max_expression > 0.3) %>%
  dplyr::mutate(norm_exp = expression/max_expression) 

df_scale_ts %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  ggplot(.) +
  geom_tile(color='black', size=0.1) +
  aes(y=factor(gene, levels=rev(TF_conser)), x=factor(cluster, levels=anno_cluster), fill=norm_exp) +
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

#ggsave(file="Plots/Manuscript/Fig2_heatmap_ts_gene_cluster_exp.pdf", width=10, height= 8)

df_scale_ts$cluster <- factor(df_scale_ts$cluster, levels= anno_cluster)

mt_scale_ts <- df_scale_ts %>%
  dplyr::filter(cluster %in% anno_cluster) %>%
  dplyr::select(gene, cluster, norm_exp) %>%
  tidyr::spread(key='gene', value='norm_exp') %>%
  dplyr::select(-cluster) %>%
  as.matrix()

rownames(mt_scale_ts) <- anno_cluster

plot_ts <- pheatmap::pheatmap(mt_scale_ts, angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                              annotation_names_col = T, show_colnames = T, cluster_rows = T,
                              treeheight_col =20, treeheight_row = 20, border_color = "black", na_col = "grey90")







#############





df_cc_cor_gather_pct_join %>%
  dplyr::filter(pair!="secNoni", pct_exp > 0.1) %>%
  ggplot(.) +
  geom_boxplot(alpha=0.9) +
  aes(x=pair, y=cor) +
  theme_bw() +
  facet_grid(~value)





## GO term plot for conserved genes ##

bg_go_BP_select <- bg_go_BP@result %>%
  dplyr::select(ID, Description, GeneRatio_bg=GeneRatio,BgRatio, p.adjust.bg = p.adjust) %>%
  dplyr::mutate(group="bg")

ego_cons_rho_select <- ego_conser_rho_BP@result %>%
  dplyr::select(ID, Description, GeneRatio_cons=GeneRatio,BgRatio, p.adjust.cons = p.adjust) %>%
  dplyr::mutate(group="conserved(rho)")

go_comp_cons <- ego_cons_rho_select %>%
  dplyr::left_join(., bg_go_BP_select, by=c('ID','Description','BgRatio')) %>%
  tidyr::separate(col=GeneRatio_cons, into=c("A1","A2"),sep="/") %>%
  dplyr::mutate(GeneRatio_cons=as.numeric(A1)/as.numeric(A2))%>%
  tidyr::separate(col=GeneRatio_bg, into=c("B1","B2"),sep="/") %>%
  dplyr::mutate(GeneRatio_bg=as.numeric(B1)/as.numeric(B2)) %>%
  dplyr::mutate(ratio_diff = GeneRatio_cons/GeneRatio_bg)

go_comp_cons_filtered <- go_comp_cons %>%
  dplyr::filter(p.adjust.cons < 0.05 & as.numeric(A1) >= 5) %>%
  dplyr::top_n(n=10, wt=ratio_diff)

df_GO_enrich_BP_cons <- rbind(dplyr::mutate(ego_conser_rho_BP@result, group="conserved(rho)"), 
                                dplyr::mutate(bg_go_BP@result, group="bg")) %>%
  tidyr::separate(GeneRatio, c("N_div_class", "N_div_all")) %>%
  tidyr::separate(BgRatio, c("Total_class", "Total_genome")) %>%
  dplyr::mutate(N_div_class = as.numeric(N_div_class), N_div_all = as.numeric(N_div_all), Total_class = as.numeric(Total_class), Total_genome = as.numeric(Total_genome)) %>%
  dplyr::mutate(GeneRatio = N_div_class/N_div_all, BgRatio = Total_class/Total_genome) %>%
  dplyr::mutate(EnrichRatio = GeneRatio/BgRatio) %>%
  dplyr::arrange(-EnrichRatio) %>%
  dplyr::filter(ID %in% go_comp_cons_filtered$ID)

GO_list_BP_cons <- unique(df_GO_enrich_BP_cons$Description)

GO_list_BP_cons_plotpoint <- data.frame(Description=GO_list_BP_cons, plotpoint=length(GO_list_BP_cons):1)

df_GO_enrich_BP_cons_sum <- df_GO_enrich_BP_cons %>%
  dplyr::filter(Description %in% GO_list_BP_cons) %>%
  dplyr::left_join(., GO_list_BP_cons_plotpoint, by = "Description") %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(class_gene_total = max(Total_class)) %>%
  dplyr::ungroup()

df_class_total_BP_cons <- df_GO_enrich_BP_cons_sum %>%
  dplyr::distinct(Description, class_gene_total, plotpoint) %>%
  dplyr::arrange(-plotpoint)

plot_GO_BP_cons <- df_GO_enrich_BP_cons_sum %>%
  ggplot(.) +
  geom_vline(xintercept = 1, color='blue', size=0.4) +
  geom_point(alpha=0.9, size=2) +
  aes(x=EnrichRatio, y=plotpoint, shape=group, fill=-log10(p.adjust), alpha=0.9) +
  theme_bw() +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10), expand = c(0.2,0.2)) +
  scale_y_continuous(breaks = length(GO_list_BP_cons):1, labels = GO_list_BP_cons, name = "", expand = c(0.05,0.05)) +
  #scale_y_continuous(breaks = length(GO_list_BP_cons):1, labels = GO_list_BP_cons, name = "",cons.axis = dup_axis(labels = df_class_total_BP$class_gene_total, name = "Total gene counts")) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  scale_shape_manual(values=c(21, 22,23,24)) +
  theme(axis.text = element_text(size=9, color='black'), axis.title = element_text(size=10, color='black'), 
        legend.title = element_text(size=9, color='black'), legend.text = element_text(size=8, color='black'), 
        legend.position = 'bottom', panel.grid = element_blank(),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.05, 'in')) +
  guides(shape= guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3)), 
         size= guide_legend(nrow=1, order = 2, title.position = "top"), 
         fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2)) +
  labs(x="Fold enrichment", shape = "Comparison", fill = "-log10(Corrected P-value)") 

plot_GO_BP_cons

ggsave("Plots/Manuscript/Fig5_cons_GO_BP.pdf", width=7, height=2.2)


## GO - MF ##

bg_go_MF_select <- bg_go_MF@result %>%
  dplyr::select(ID, Description, GeneRatio_bg=GeneRatio,BgRatio, p.adjust.bg = p.adjust) %>%
  dplyr::mutate(group="bg")

ego_cons_rho_select <- ego_conser_rho_MF@result %>%
  dplyr::select(ID, Description, GeneRatio_cons=GeneRatio,BgRatio, p.adjust.cons = p.adjust) %>%
  dplyr::mutate(group="conserved(rho)")

go_comp_cons <- ego_cons_rho_select %>%
  dplyr::left_join(., bg_go_MF_select, by=c('ID','Description','BgRatio')) %>%
  tidyr::separate(col=GeneRatio_cons, into=c("A1","A2"),sep="/") %>%
  dplyr::mutate(GeneRatio_cons=as.numeric(A1)/as.numeric(A2))%>%
  tidyr::separate(col=GeneRatio_bg, into=c("B1","B2"),sep="/") %>%
  dplyr::mutate(GeneRatio_bg=as.numeric(B1)/as.numeric(B2)) %>%
  dplyr::mutate(ratio_diff = GeneRatio_cons/GeneRatio_bg)

go_comp_cons_filtered <- go_comp_cons %>%
  dplyr::filter(p.adjust.cons < 0.05 & as.numeric(A1) >= 5) %>%
  dplyr::top_n(n=10, wt=ratio_diff)

df_GO_enrich_MF_cons <- rbind(dplyr::mutate(ego_conser_rho_MF@result, group="conserved(rho)"), 
                              dplyr::mutate(bg_go_MF@result, group="bg")) %>%
  tidyr::separate(GeneRatio, c("N_div_class", "N_div_all")) %>%
  tidyr::separate(BgRatio, c("Total_class", "Total_genome")) %>%
  dplyr::mutate(N_div_class = as.numeric(N_div_class), N_div_all = as.numeric(N_div_all), Total_class = as.numeric(Total_class), Total_genome = as.numeric(Total_genome)) %>%
  dplyr::mutate(GeneRatio = N_div_class/N_div_all, BgRatio = Total_class/Total_genome) %>%
  dplyr::mutate(EnrichRatio = GeneRatio/BgRatio) %>%
  dplyr::arrange(-EnrichRatio) %>%
  dplyr::filter(ID %in% go_comp_cons_filtered$ID)

GO_list_MF_cons <- unique(df_GO_enrich_MF_cons$Description)

GO_list_MF_cons_plotpoint <- data.frame(Description=GO_list_MF_cons, plotpoint=length(GO_list_MF_cons):1)

df_GO_enrich_MF_cons_sum <- df_GO_enrich_MF_cons %>%
  dplyr::filter(Description %in% GO_list_MF_cons) %>%
  dplyr::left_join(., GO_list_MF_cons_plotpoint, by = "Description") %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(class_gene_total = max(Total_class)) %>%
  dplyr::ungroup()

df_class_total_MF_cons <- df_GO_enrich_MF_cons_sum %>%
  dplyr::distinct(Description, class_gene_total, plotpoint) %>%
  dplyr::arrange(-plotpoint)

plot_GO_MF_cons <- df_GO_enrich_MF_cons_sum %>%
  ggplot(.) +
  geom_vline(xintercept = 1, color='blue', size=0.4) +
  geom_point(alpha=0.9, size=2) +
  aes(x=EnrichRatio, y=plotpoint, shape=group, fill=-log10(p.adjust), alpha=0.9) +
  theme_bw() +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10), expand = c(0.2,0.2)) +
  scale_y_continuous(breaks = length(GO_list_MF_cons):1, labels = GO_list_MF_cons, name = "", expand = c(0.05,0.05)) +
  #scale_y_continuous(breaks = length(GO_list_MF_cons):1, labels = GO_list_MF_cons, name = "",cons.axis = dup_axis(labels = df_class_total_MF$class_gene_total, name = "Total gene counts")) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  scale_shape_manual(values=c(21, 22,23,24)) +
  theme(axis.text = element_text(size=9, color='black'), axis.title = element_text(size=10, color='black'), 
        legend.title = element_text(size=9, color='black'), legend.text = element_text(size=8, color='black'), 
        legend.position = 'bottom', panel.grid = element_blank(),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.05, 'in')) +
  guides(shape= guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3)), 
         size= guide_legend(nrow=1, order = 2, title.position = "top"), 
         fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2)) +
  labs(x="Fold enrichment", shape = "Comparison", fill = "-log10(Corrected P-value)") 

plot_GO_MF_cons

ggsave("Plots/Manuscript/Fig5_cons_GO_MF.pdf", width=7, height=2.2)

####### archive

## average expression ##

df_per_exp_mel <- test_degs_melsec %>%
  dplyr::select(gene, pct_exp=pct.1)

df_cc_cor %>%
  dplyr::filter(pair == 'melsec') %>%
  dplyr::left_join(., df_per_exp_mel, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_point(size=1.5, alpha=0.8, shape=21) +
  aes(x=r, y=rho, fill=pct_exp) +
  theme_bw() +
  scale_fill_gradient(low = "blue", high = "red", na.value = NA)


df_cc_cor %>%
  dplyr::filter(pair == 'melsec') %>%
  dplyr::left_join(., df_per_exp_mel, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_point(size=1.5, alpha=0.8, shape=21) +
  aes(x=pct_exp, y=rho) +
  theme_bw()

df_cc_cor %>%
  dplyr::filter(pair == 'melsec') %>%
  dplyr::left_join(., df_per_exp_mel, by='gene') %>%
  na.omit() %>%
  ggplot(.) +
  geom_point(size=1.5, alpha=0.8, shape=21) +
  aes(x=pct_exp, y=r) +
  theme_bw()




## average expression ##

df_per_exp_mel <- test_degs_melsec %>%
  dplyr::select(gene, pct_exp=pct.1)

df_cc_cor_gather_av_exp_melsec <- df_cc_cor_gather %>%
  dplyr::filter(pair == "melsec") %>%
  dplyr::left_join(., df_per_exp_mel, by='gene')

df_cc_cor_gather_av_exp_melsec %>%
  ggplot(.) +
  geom_point(size=1, alpha=0.8) +
  geom_text_repel(data=dplyr::filter(df_cc_cor_gather_av_exp_melsec, gene=="unc-13"), aes(label=gene)) +
  #geom_smooth() +
  aes(y=pct_exp*100, x=cor) +
  theme_bw() +
  theme(panel.grid=element_blank())

## top 50 expression level, top 50 correlation, DEGs, comparison


## scatter plot - x axis: N of exp clusters, y axis : correlation ##

df_exp_clusters_average_refOWN <- df_exp_clusters %>%
  dplyr::filter(ref=="OWN") %>%
  dplyr::group_by(gene, )

