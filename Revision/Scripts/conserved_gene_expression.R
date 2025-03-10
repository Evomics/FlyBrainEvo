library(Seurat)
library(tidyverse)
library(ggrepel)
library(cowplot)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

TrioBrain.integrated_slim_labeled_final <- readRDS(file = "Processed_Data/TrioBrain.integrated_slim_labeled_final.rds")

load("Processed_Data/df_per_expressed_ref.RData")
load("Processed_Data/df_sum_exp_SCT.RData")
load("Processed_Data/df_deg_sigs.RData")
load("Processed_Data/Trio_Dmelref_DF_labeled.RData")

### select genes with >%5 expression in Dmel

Trio_Dmelref_DF_labeled_species <- Trio_Dmelref_DF_labeled
Trio_Dmelref_DF_labeled_species$orig.ident <- gsub("_rep[1-6]", "", Trio_Dmelref_DF_labeled_species$orig.ident)
Trio_Dmelref_DF_labeled_species@meta.data$orig.ident <- factor(Trio_Dmelref_DF_labeled_species@meta.data$orig.ident, 
                                                                               levels=c("Dmel", "Dsim_to_DmelRef", "Dsec_to_DmelRef"))

DotPlot_Trio_DmelRef <- DotPlot(Trio_Dmelref_DF_labeled_species, group.by='orig.ident', feature=rownames(Trio_Dmelref_DF_labeled_species[["SCT"]]))$data

TrioBrain.integrated_slim_labeled_final_species <- TrioBrain.integrated_slim_labeled_final
TrioBrain.integrated_slim_labeled_final_species$orig.ident <- gsub("_rep[1-6]", "", TrioBrain.integrated_slim_labeled_final_species$orig.ident)
TrioBrain.integrated_slim_labeled_final_species@meta.data$orig.ident <- factor(TrioBrain.integrated_slim_labeled_final_species@meta.data$orig.ident, 
                                                                               levels=c("Dmel", "Dsim", "Dsec", "DsecNoni"))

DotPlot_Trio_Self <- DotPlot(TrioBrain.integrated_slim_labeled_final_species, group.by='orig.ident', feature=rownames(Dmel_seurat[["SCT"]]))$data

Dsim_per_exp <- DotPlot_Trio_DmelRef %>%
  dplyr::filter(id=="Dsim_to_DmelRef") %>%
  dplyr::select(features.plot, pct.exp_DmelRef = pct.exp) %>%
  dplyr::left_join(., dplyr::filter(DotPlot_Trio_Self, id=="Dsim"), by='features.plot') %>%
  dplyr::select(features.plot, pct.exp, pct.exp_DmelRef) %>%
  dplyr::mutate(pct_diff =abs(pct.exp-pct.exp_DmelRef)) %>%
  dplyr::mutate(RefMEL_rank_exp=min_rank(-pct.exp_DmelRef), RefOWN_rank_exp=min_rank(-pct.exp)) %>%
  dplyr::mutate(rank_diff = abs(RefMEL_rank_exp-RefOWN_rank_exp), rank_diff_ratio=rank_diff/n())

sim_per_exp_05 <- Dsim_per_exp %>%
  dplyr::filter(rank_diff_ratio < 0.05) %>%
  dplyr::distinct(features.plot) %>%
  pull()

Dsec_per_exp <- DotPlot_Trio_DmelRef %>%
  dplyr::filter(id=="Dsec_to_DmelRef") %>%
  dplyr::select(features.plot, pct.exp_DmelRef = pct.exp) %>%
  dplyr::left_join(., dplyr::filter(DotPlot_Trio_Self, id=="Dsec"), by='features.plot') %>%
  dplyr::select(features.plot, pct.exp, pct.exp_DmelRef) %>%
  dplyr::mutate(pct_diff =abs(pct.exp-pct.exp_DmelRef)) %>%
  dplyr::mutate(RefMEL_rank_exp=min_rank(-pct.exp_DmelRef), RefOWN_rank_exp=min_rank(-pct.exp)) %>%
  dplyr::mutate(rank_diff = abs(RefMEL_rank_exp-RefOWN_rank_exp), rank_diff_ratio=rank_diff/n())

sec_per_exp_05 <- Dsec_per_exp %>%
  dplyr::filter(rank_diff_ratio < 0.05) %>%
  dplyr::distinct(features.plot) %>%
  pull()

simsec_clean_genes <- intersect(sim_per_exp_05, sec_per_exp_05)

Dmel_05_select <- DotPlot_Trio_Self %>%
  dplyr::filter(pct.exp >= 5, id=="Dmel") %>%
  dplyr::distinct(features.plot) %>%
  pull()

simsec_clean_Dmel05_genes <- intersect(Dmel_05_select, simsec_clean_genes) ## 2405 genes are selected for the analysis

save(DotPlot_Trio_DmelRef, DotPlot_Trio_Self, simsec_clean_Dmel05_genes, file = "Processed_Data/Trio_gene_pct_exp.RData")

load("Processed_Data/Trio_gene_pct_exp.RData")

## specificity of these 1686 genes ##

df_per_expressed_ref_cor <- df_per_expressed_ref %>%
  dplyr::filter(gene %in% simsec_clean_Dmel05_genes) %>%
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

### correlation ###

## scatter plots ##

df_cc_cor_rho <- NULL
df_cc_cor_r <- NULL

for (gen in simsec_clean_Dmel05_genes) {
  
  df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
    dplyr::filter(gene==gen, cluster %in% Anno_idents) %>%
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
gen='gw'

df_sum_exp_SCT_select <- df_sum_exp_SCT %>%
  dplyr::filter(gene==gen, cluster %in% Anno_idents) %>%
  tidyr::spread(key='species', value='expression')

df_sum_exp_SCT_select %>%
  ggplot(.) +
  geom_smooth(size= 0.3, method="glm") +
  geom_point(alpha = 0.6, size=1.2) +
  #geom_text_repel(aes(label=cluster)) +
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

df_per_select_Dmel <- DotPlot_Trio_Self %>%
  dplyr::filter(pct.exp >= 5, id=="Dmel") %>%
  dplyr::select(gene=features.plot, pct_exp=pct.exp)

df_cc_cor_gather_pct_join <- df_cc_cor_gather %>%
  dplyr::left_join(., df_per_select_Dmel, by='gene') %>%
  na.omit()

#save(df_cc_cor_gather_pct_join, file="Processed_Data/df_cc_cor_gather_pct_join.RData")
load("Processed_Data/df_cc_cor_gather_pct_join.RData")


### density plot ###

df_cc_cor_gather_pct_join %>%
  dplyr::filter(pair!="secNoni", pct_exp > 0.05) %>%
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
  dplyr::filter(pair!="secNoni", pct_exp > 0.05, value=='rho') %>%
  ggplot(.) +
  geom_density_2d_filled(alpha=0.85) +
  geom_point(alpha=0.5, size =0.01) +
  geom_hline(yintercept=0.7, linetype=2, color='red', alpha=0.8, size=0.5) +
  aes(x=pct_exp, y=cor) +
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
  dplyr::filter(pair!="secNoni", cor > 0.7, pct_exp > 0.05, value=='rho') %>%
  ggplot(.) +
  geom_histogram(binwidth=3, boundary=10) +
  aes(x=pct_exp) +
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

ggsave("Plots/Manuscript/Fig5_2D_density_rho_pct.pdf", width=10, height = 5)

df_cc_cor_gather_pct_join_strong_rho <- df_cc_cor_gather_pct_join %>%
  dplyr::filter(value=='rho', pair != "secNoni") %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(min(cor) > 0.7) %>%
  dplyr::ungroup()

length(unique(df_cc_cor_gather_pct_join_strong_rho$gene))

conser_genes <- unique(df_cc_cor_gather_pct_join_strong_rho$gene)
write.csv(conser_genes, file="Processed_Data/conserved_genes.csv", row.names=F, quote=F)
conser_genes_control <- unique(df_cc_cor_gather_pct_join$gene)
write.csv(conser_genes_control, file="Processed_Data/conserved_genes_control.csv", row.names=F, quote=F)

df_high_cor_summary <- df_cc_cor_gather_pct_join %>%
  dplyr::filter(gene %in% conser_genes, pair != "secNoni") %>%
  dplyr::select(gene, pair, correlation=value, value=cor, percent_expressed=pct_exp) %>%
  dplyr::arrange(gene, pair)

write.csv(df_high_cor_summary, file="Processed_Data/df_high_cor_summary.csv", quote=F, row.names=F)



## FlyEnrichR ###

### conserved 411 https://maayanlab.cloud/FlyEnrichr/enrich?dataset=8ae5d23c524053d3ccb827a1799dd772 ##

FER_cons_BP <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved413/GO_Biological_Process_2018_table.txt") %>%
  dplyr::mutate(class="BP", group='conserved', n_gene=411) %>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_cons_MF <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved413/GO_Molecular_Function_2018_table.txt") %>%
  dplyr::mutate(class="MF", group='conserved', n_gene=411)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_cons_CC <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved413/GO_Cellular_Component_2018_table.txt") %>%
  dplyr::mutate(class="CC", group='conserved', n_gene=411)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')

### control 2405 https://maayanlab.cloud/FlyEnrichr/enrich?dataset=be3d48b392bcc6187874b1c932cc630f 

FER_control_BP <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved_control2405/GO_Biological_Process_2018_table.txt") %>%
  dplyr::mutate(class="BP", group='control', n_gene=2405)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_control_MF <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved_control2405/GO_Molecular_Function_2018_table.txt") %>%
  dplyr::mutate(class="MF", group='control', n_gene=2405)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_control_CC <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved_control2405/GO_Cellular_Component_2018_table.txt") %>%
  dplyr::mutate(class="CC", group='control', n_gene=2405)%>%
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
  dplyr::mutate(fold_enrich = ifelse(n_Term==2, 2405/411*n_overlap_conserved/n_overlap_control, NA)) %>%
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

plot_FER_cons_GO <- cowplot::plot_grid(plot_FER_cons_BP, plot_FER_cons_MF, plot_FER_cons_CC, ncol=1, align='v', rel_heights = c(1,1,0.6))

plot_FER_cons_GO

ggsave(plot_FER_cons_GO, file="Plots/Manuscript/Extended_Data_Fig2c_ver4.pdf", width=12, height = 7.5)


### Different threshold GO term. rho = 0.6 or 0.8 ###

 # rho > 0.6 #

### conserved 656 https://maayanlab.cloud/FlyEnrichr/enrich?dataset=a93ec83b5decfaec023f050bd2c2c8d6 ## 

FER_cons_BP_60 <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved656/GO_Biological_Process_2018_table.txt") %>%
  dplyr::mutate(class="BP", group='conserved', n_gene=656) %>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_cons_MF_60 <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved656/GO_Molecular_Function_2018_table.txt") %>%
  dplyr::mutate(class="MF", group='conserved', n_gene=656)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_cons_CC_60 <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved656/GO_Cellular_Component_2018_table.txt") %>%
  dplyr::mutate(class="CC", group='conserved', n_gene=656)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')

df_FER_cons_60_int <- rbind(FER_cons_BP_60, FER_cons_MF_60, FER_cons_CC_60, FER_control_BP, FER_control_MF, FER_control_CC) %>%
  tidyr::separate(Overlap, c('n_overlap','n_gene_class')) %>%
  dplyr::filter(as.numeric(n_overlap) >=5) %>%
  dplyr::mutate(frac_to_group=as.numeric(n_overlap)/n_gene, frac_to_class=as.numeric(n_overlap)/as.numeric(n_gene_class)) %>%
  dplyr::mutate(frac_to_group_norm_class=as.numeric(n_overlap)/n_gene/as.numeric(n_gene_class)) %>%
  dplyr::group_by(Term) %>%
  dplyr::mutate(n_Term = n(),
                n_overlap_conserved=min(as.numeric(n_overlap)), n_overlap_control=max(as.numeric(n_overlap)), 
                fraction=n_overlap_conserved/n_overlap_control) %>%
  dplyr::mutate(fold_enrich = ifelse(n_Term==2, 2405/656*n_overlap_conserved/n_overlap_control, NA)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Genes) %>%
  dplyr::filter(p_adj==min(p_adj) & Z_score==min(Z_score)) %>%
  dplyr::ungroup()

## add fold enrichment

df_FER_cons_60_int_top10 <- df_FER_cons_60_int %>%
  dplyr::filter(group == "conserved") %>%
  dplyr::group_by(class) %>%
  dplyr::top_n(10, fold_enrich)

top10_term_60 <- unique(df_FER_cons_60_int_top10$Term)

plot_FER_cons_60_BP <- df_FER_cons_60_int %>%
  dplyr::filter(Term %in% top10_term_60, group == "conserved", fold_enrich > 1, class=="BP") %>%
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
  labs(fill='-log10(adjusted P-value)', x='Fold Enrichment') +
  lims(x=c(1,4))

plot_FER_cons_60_BP

plot_FER_cons_60_MF <- df_FER_cons_60_int %>%
  dplyr::filter(Term %in% top10_term_60, group == "conserved", fold_enrich > 1, class=="MF") %>%
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
  labs(fill='-log10(adjusted P-value)', x='Fold Enrichment') +
  lims(x=c(1,4))

plot_FER_cons_60_MF

plot_FER_cons_60_CC <- df_FER_cons_60_int %>%
  dplyr::filter(Term %in% top10_term_60, group == "conserved", fold_enrich > 1, class=="CC") %>%
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
  labs(fill='-log10(adjusted P-value)', x='Fold Enrichment') +
  lims(x=c(1,4))

plot_FER_cons_60_CC

plot_FER_cons_60_GO <- cowplot::plot_grid(plot_FER_cons_60_BP, plot_FER_cons_60_MF, plot_FER_cons_60_CC, ncol=1, align='v', rel_heights = c(1,1,0.9))

plot_FER_cons_60_GO

ggsave(plot_FER_cons_60_GO, file="Plots/Manuscript/Revision_fig_rho60.pdf", width=12, height = 9)


# rho > 0.8 #

### conserved 185 https://maayanlab.cloud/FlyEnrichr/enrich?dataset=1ce25fe9a49e99b9bb29984a24320515 ## 

FER_cons_BP_80 <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved185/GO_Biological_Process_2018_table.txt") %>%
  dplyr::mutate(class="BP", group='conserved', n_gene=185) %>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_cons_MF_80 <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved185/GO_Molecular_Function_2018_table.txt") %>%
  dplyr::mutate(class="MF", group='conserved', n_gene=185)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')
FER_cons_CC_80 <- data.table::fread(file="Processed_Data/FlyEnrichR/conserved185/GO_Cellular_Component_2018_table.txt") %>%
  dplyr::mutate(class="CC", group='conserved', n_gene=185)%>%
  dplyr::rename(p_adj = 'Adjusted P-value', Z_score = 'Z-score')

df_FER_cons_80_int <- rbind(FER_cons_BP_80, FER_cons_MF_80, FER_cons_CC_80, FER_control_BP, FER_control_MF, FER_control_CC) %>%
  tidyr::separate(Overlap, c('n_overlap','n_gene_class')) %>%
  dplyr::filter(as.numeric(n_overlap) >=5) %>%
  dplyr::mutate(frac_to_group=as.numeric(n_overlap)/n_gene, frac_to_class=as.numeric(n_overlap)/as.numeric(n_gene_class)) %>%
  dplyr::mutate(frac_to_group_norm_class=as.numeric(n_overlap)/n_gene/as.numeric(n_gene_class)) %>%
  dplyr::group_by(Term) %>%
  dplyr::mutate(n_Term = n(),
                n_overlap_conserved=min(as.numeric(n_overlap)), n_overlap_control=max(as.numeric(n_overlap)), 
                fraction=n_overlap_conserved/n_overlap_control) %>%
  dplyr::mutate(fold_enrich = ifelse(n_Term==2, 2405/185*n_overlap_conserved/n_overlap_control, NA)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Genes) %>%
  dplyr::filter(p_adj==min(p_adj) & Z_score==min(Z_score)) %>%
  dplyr::ungroup()

## add fold enrichment

df_FER_cons_80_int_top10 <- df_FER_cons_80_int %>%
  dplyr::filter(group == "conserved") %>%
  dplyr::group_by(class) %>%
  dplyr::top_n(10, fold_enrich)

top10_term_80 <- unique(df_FER_cons_80_int_top10$Term)

plot_FER_cons_80_BP <- df_FER_cons_80_int %>%
  dplyr::filter(Term %in% top10_term_80, group == "conserved", fold_enrich > 1, class=="BP") %>%
  ggplot(.) +
  geom_point(shape=21, size=3) +
  aes(x=fold_enrich, y=reorder(Term, fold_enrich), fill=-log10(p_adj)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title.x= element_text(color='black', size=11),
        legend.key.height = unit(0.1,"in"),
        legend.position=c(0.72,0.2),
        legend.direction='horizontal',
        legend.title = element_text(color='black', size=9),
        legend.text = element_text(color='black', size=8.5)) +
  facet_wrap(~class, ncol=1, scales='free') +
  labs(fill='-log10(adjusted P-value)', x='Fold Enrichment') 

plot_FER_cons_80_BP

plot_FER_cons_80_MF <- df_FER_cons_80_int %>%
  dplyr::filter(Term %in% top10_term_80, group == "conserved", fold_enrich > 1, class=="MF") %>%
  ggplot(.) +
  geom_point(shape=21, size=3) +
  aes(x=fold_enrich, y=reorder(Term, fold_enrich), fill=-log10(p_adj)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title.x= element_text(color='black', size=11),
        legend.key.height = unit(0.1,"in"),
        legend.position=c(0.72,0.2),
        legend.direction='horizontal',
        legend.title = element_text(color='black', size=9),
        legend.text = element_text(color='black', size=8.5)) +
  facet_wrap(~class, ncol=1, scales='free') +
  labs(fill='-log10(adjusted P-value)', x='Fold Enrichment') 

plot_FER_cons_80_MF

plot_FER_cons_80_CC <- df_FER_cons_80_int %>%
  dplyr::filter(Term %in% top10_term_80, group == "conserved", fold_enrich > 1, class=="CC") %>%
  ggplot(.) +
  geom_point(shape=21, size=3) +
  aes(x=fold_enrich, y=reorder(Term, fold_enrich), fill=-log10(p_adj)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title.x= element_text(color='black', size=11),
        legend.key.height = unit(0.1,"in"),
        legend.position=c(0.72,0.2),
        legend.direction='horizontal',
        legend.title = element_text(color='black', size=9),
        legend.text = element_text(color='black', size=8.5)) +
  facet_wrap(~class, ncol=1, scales='free') +
  labs(fill='-log10(adjusted P-value)', x='Fold Enrichment')

plot_FER_cons_80_CC

plot_FER_cons_80_GO <- cowplot::plot_grid(plot_FER_cons_80_BP, plot_FER_cons_80_MF, plot_FER_cons_80_CC, ncol=1, align='v', rel_heights = c(1,1,0.6))

plot_FER_cons_80_GO

ggsave(plot_FER_cons_80_GO, file="Plots/Manuscript/Revision_fig_rho80.pdf", width=12, height = 8.5)

