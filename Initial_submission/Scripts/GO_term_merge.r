library(clusterProfiler)
library("org.Dm.eg.db")
library(enrichplot)
library(ggupset)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("Processed_Data/shared_genes.RData")
load("Processed_Data/background_vari_features.RData")
load("Processed_Data/DEG_cluster_ref_sigs_genelist.RData")
load("Processed_Data/df_cc_cor_gather_pct_join.RData")


## GO term enrichment ##

## background ##

bg_go_BP <- enrichGO(gene          = backgroud_DEG,
                     universe      = shared_genes,
                     OrgDb         = org.Dm.eg.db,
                     keyType = "SYMBOL",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = F)

dotplot(bg_go_BP, showCategory = 20)
upsetplot(bg_go_BP, n = 15)
View(bg_go_BP@result)

bg_go_MF <- enrichGO(gene          = backgroud_DEG,
                     universe      = shared_genes,
                     OrgDb         = org.Dm.eg.db,
                     keyType = "SYMBOL",
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = F)

dotplot(bg_go_MF, showCategory = 20)
upsetplot(bg_go_MF, n = 15)

## DEGs ##

## speciation genes ##

deg_secspe_go_BP <- enrichGO(gene          = Dsec_spe_genes,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05,
                             readable      = F)

dotplot(deg_secspe_go_BP, showCategory = 20)

View(deg_secspe_go_BP@result)

deg_secspe_go_MF <- enrichGO(gene          = Dsec_spe_genes,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "MF",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05,
                             readable      = F)

dotplot(deg_secspe_go_MF, showCategory = 20)

View(deg_secspe_go_MF@result)

## speciation genes - glia ##


deg_secspe_glia_go_BP <- enrichGO(gene          = Dsec_spe_glia,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05,
                             readable      = F)

dotplot(deg_secspe_glia_go_BP, showCategory = 20)

View(deg_secspe_glia_go_BP@result)

deg_secspe_glia_go_MF <- enrichGO(gene          = Dsec_spe_glia,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "MF",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05,
                             readable      = F)

dotplot(deg_secspe_glia_go_MF, showCategory = 20)

View(deg_secspe_glia_go_MF@result)

deg_secspe_glia_go_All <- enrichGO(gene          = Dsec_spe_glia,
                                  universe      = shared_genes,
                                  OrgDb         = org.Dm.eg.db,
                                  keyType = "SYMBOL",
                                  ont           = "ALL",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.05,
                                  readable      = F)

View(deg_secspe_glia_go_All@result)

deg_secspe_nonglia_go_All <- enrichGO(gene          = Dsec_spe_non_glia,
                                   universe      = shared_genes,
                                   OrgDb         = org.Dm.eg.db,
                                   keyType = "SYMBOL",
                                   ont           = "ALL",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = F)

View(deg_secspe_nonglia_go_All@result)


## KEGG analysis ##

df_mel_genename <- read.table(file="Dmel/Dmel_ID_genename.tsv")

search_kegg_organism('dme', by='kegg_code')

kk <- enrichKEGG(gene = paste("Dmel_", Dsec_spe_glia, sep=""), organism = 'dme')

head(kk)

entrez_id_Dsec_spe_glia <- as.data.frame(org.Dm.egALIAS2EG) %>%
  dplyr::filter(alias_symbol %in% Dsec_spe_glia) %>%
  dplyr::filter(!gene_id %in% c(32627, 43541, 45233, 42249,43825, 42424,
                                40336, 37890, 45845))



KEGGdb <- download_KEGG('dme', keggType = "KEGG", keyType = "kegg")

KEGG_genenames <- KEGGdb$KEGGPATHID2EXTID$to

"Dmel_CG9470" %in% KEGG_genenames

entrez_id_KEGGname <- as.data.frame(org.Dm.egALIAS2EG) %>%
  dplyr::filter(alias_symbol %in% KEGG_genenames)






## sim ##

deg_simspe_go_BP <- enrichGO(gene          = Dsim_spe_genes,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05,
                             readable      = F)

dotplot(deg_simspe_go_BP, showCategory = 20)

View(deg_simspe_go_BP@result)

deg_simspe_go_MF <- enrichGO(gene          = Dsim_spe_genes,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "MF",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05,
                             readable      = F)

dotplot(deg_simspe_go_MF, showCategory = 20)

View(deg_simspe_go_MF@result)


## BP ##

deg_melsec_go_BP <- enrichGO(gene          = DEG_melsec_celltype,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = F)

dotplot(deg_melsec_go_BP, showCategory = 20)

deg_melsim_go_BP <- enrichGO(gene          = DEG_melsim_celltype,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = F)

dotplot(deg_melsim_go_BP, showCategory = 20)


deg_simsec_go_BP <- enrichGO(gene          = DEG_simsec_celltype,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = F)

dotplot(deg_simsec_go_BP, showCategory = 20)

## species specific ##

DEG_sec_specific <- intersect(DEG_simsec_celltype, DEG_melsec_celltype)

Dsec_specific_BP <- enrichGO(gene          = DEG_sec_specific,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = F)

View(Dsec_specific_BP@result)
dotplot(Dsec_specific_BP, showCategory=10)

DEG_sim_specific <- intersect(DEG_simsec_celltype, DEG_melsim_celltype)

Dsim_specific_BP <- enrichGO(gene          = DEG_sim_specific,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = F)

View(Dsim_specific_BP@result)
dotplot(Dsim_specific_BP, showCategory=10)

## MF ##

deg_melsec_go_MF <- enrichGO(gene          = DEG_melsec_celltype,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "MF",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = F)

dotplot(deg_melsec_go_MF, showCategory = 20)

deg_melsim_go_MF <- enrichGO(gene          = DEG_melsim_celltype,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "MF",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = F)

dotplot(deg_melsim_go_MF, showCategory = 20)


deg_simsec_go_MF <- enrichGO(gene          = DEG_simsec_celltype,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "MF",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = F)

dotplot(deg_simsec_go_MF, showCategory = 20)


## species specific ##

Dsec_specific_MF <- enrichGO(gene          = DEG_sec_specific,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "MF",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = F)

View(Dsec_specific_MF@result)
dotplot(Dsec_specific_MF, showCategory=10)

Dsim_specific_MF <- enrichGO(gene          = DEG_sim_specific,
                             universe      = shared_genes,
                             OrgDb         = org.Dm.eg.db,
                             keyType = "SYMBOL",
                             ont           = "MF",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = F)

View(Dsim_specific_MF@result)
dotplot(Dsim_specific_MF, showCategory=10)


## conserved genes ##

## conserved ##

cor_threshold=0.8

df_conserved_exp_gene <- df_cc_cor_gather_pct_join %>%
  dplyr::filter(DEG==F, cor > cor_threshold, pair != "secNoni") %>%
  dplyr::group_by(gene, value) %>%
  dplyr::filter(n()==3) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(gene)


ego_conser_BP <- enrichGO(gene          = df_conserved_exp_gene$gene,
                          universe      = shared_genes,
                          OrgDb         = org.Dm.eg.db,
                          keyType = "SYMBOL",
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = F)

dotplot(ego_conser_BP, showCategory = 20)
upsetplot(ego_conser_BP, n = 15)
View(ego_conser_BP@result)

ego_conser_MF <- enrichGO(gene          = df_conserved_exp_gene$gene,
                          universe      = shared_genes,
                          OrgDb         = org.Dm.eg.db,
                          keyType = "SYMBOL",
                          ont           = "MF",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = F)

dotplot(ego_conser_MF, showCategory = 20)
upsetplot(ego_conser_MF, n = 15)
View(ego_conser_MF@result)

## rho only ##

rho_threshold=0.7

df_conserved_rho_exp_gene <- df_cc_cor_gather_pct_join %>%
  dplyr::filter(DEG==F, cor > rho_threshold, pair != "secNoni", value=='rho') %>%
  dplyr::group_by(gene, value) %>%
  dplyr::filter(n()==3) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(gene)

ego_conser_rho_BP <- enrichGO(gene          = df_conserved_rho_exp_gene$gene,
                              universe      = shared_genes,
                              OrgDb         = org.Dm.eg.db,
                              keyType = "SYMBOL",
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05,
                              readable      = F)

dotplot(ego_conser_rho_BP, showCategory = 20)
upsetplot(ego_conser_rho_BP, n = 15)
View(ego_conser_rho_BP@result)

ego_conser_rho_MF <- enrichGO(gene          = df_conserved_rho_exp_gene$gene,
                              universe      = shared_genes,
                              OrgDb         = org.Dm.eg.db,
                              keyType = "SYMBOL",
                              ont           = "MF",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05,
                              readable      = F)

dotplot(ego_conser_rho_MF, showCategory = 20)
upsetplot(ego_conser_rho_MF, n = 15)
View(ego_conser_rho_MF@result)

## merge all - GO term plot ##

## BP ##

bg_go_BP_select <- bg_go_BP@result %>%
  dplyr::select(ID, Description, GeneRatio_bg=GeneRatio,BgRatio, p.adjust.bg = p.adjust) %>%
  dplyr::mutate(group="bg")

ego_BP_melsec_select <- deg_melsec_go_BP@result %>%
  dplyr::select(ID, Description, GeneRatio_target=GeneRatio,BgRatio, p.adjust.target = p.adjust) %>%
  dplyr::mutate(group="melsec")

ego_BP_melsim_select <- deg_melsim_go_BP@result %>%
  dplyr::select(ID, Description, GeneRatio_target=GeneRatio,BgRatio, p.adjust.target = p.adjust) %>%
  dplyr::mutate(group="melsim")

ego_BP_simsec_select <- deg_simsec_go_BP@result %>%
  dplyr::select(ID, Description, GeneRatio_target=GeneRatio,BgRatio, p.adjust.target = p.adjust) %>%
  dplyr::mutate(group="simsec")

ego_BP_sim_select <- Dsim_specific_BP@result %>%
  dplyr::select(ID, Description, GeneRatio_target=GeneRatio,BgRatio, p.adjust.target = p.adjust) %>%
  dplyr::mutate(group="Dsim")

ego_BP_sec_select <- Dsec_specific_BP@result %>%
  dplyr::select(ID, Description, GeneRatio_target=GeneRatio,BgRatio, p.adjust.target = p.adjust) %>%
  dplyr::mutate(group="Dsec")

ego_cons_rho_select <- ego_conser_rho_BP@result %>%
  dplyr::select(ID, Description, GeneRatio_target=GeneRatio,BgRatio, p.adjust.target = p.adjust) %>%
  dplyr::mutate(group="conserved(rho)")

ego_BP_pair_cons_select_merge <- rbind(ego_BP_melsec_select, ego_BP_melsim_select, ego_BP_simsec_select, 
                                       ego_cons_rho_select)

go_BP_comp_pair_cons <- ego_BP_pair_cons_select_merge %>%
  dplyr::left_join(., bg_go_BP_select, by=c('ID','Description','BgRatio')) %>%
  tidyr::separate(col=GeneRatio_target, into=c("A1","A2"),sep="/") %>%
  dplyr::mutate(GeneRatio_target=as.numeric(A1)/as.numeric(A2))%>%
  tidyr::separate(col=GeneRatio_bg, into=c("B1","B2"),sep="/") %>%
  dplyr::mutate(GeneRatio_bg=as.numeric(B1)/as.numeric(B2)) %>%
  dplyr::mutate(ratio_diff = GeneRatio_target/GeneRatio_bg)

go_BP_comp_pair_cons_filtered <- go_BP_comp_pair_cons %>%
  dplyr::filter(p.adjust.target < 0.01 & as.numeric(A1) >= 10) %>%
  dplyr::group_by(group.x) %>%
  dplyr::top_n(n=10, wt=ratio_diff)

df_GO_enrich_BP_pair_cons <- rbind(dplyr::mutate(deg_melsec_go_BP@result, group="melsec"), 
                              dplyr::mutate(deg_melsim_go_BP@result, group="melsim"), 
                              dplyr::mutate(deg_simsec_go_BP@result, group="simsec"), 
                              dplyr::mutate(ego_conser_rho_BP@result, group="conserved"),
                              dplyr::mutate(bg_go_BP@result, group="bg")) %>%
  tidyr::separate(GeneRatio, c("N_div_class", "N_div_all")) %>%
  tidyr::separate(BgRatio, c("Total_class", "Total_genome")) %>%
  dplyr::mutate(N_div_class = as.numeric(N_div_class), N_div_all = as.numeric(N_div_all), Total_class = as.numeric(Total_class), Total_genome = as.numeric(Total_genome)) %>%
  dplyr::mutate(GeneRatio = N_div_class/N_div_all, BgRatio = Total_class/Total_genome) %>%
  dplyr::mutate(EnrichRatio = GeneRatio/BgRatio) %>%
  dplyr::arrange(-EnrichRatio) %>%
  dplyr::filter(ID %in% go_BP_comp_pair_cons_filtered$ID)

GO_list_BP_pair_cons <- unique(df_GO_enrich_BP_pair_cons$Description)

GO_list_BP_pair_cons_plotpoint <- data.frame(Description=GO_list_BP_pair_cons, plotpoint=length(GO_list_BP_pair_cons):1)

df_GO_enrich_BP_pair_cons_sum <- df_GO_enrich_BP_pair_cons %>%
  dplyr::filter(Description %in% GO_list_BP_pair_cons) %>%
  dplyr::left_join(., GO_list_BP_pair_cons_plotpoint, by = "Description") %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(class_gene_total = max(Total_class)) %>%
  dplyr::ungroup()

df_GO_enrich_BP_pair_cons_sum %>%
  ggplot(.) +
  geom_vline(xintercept = 1, color='blue', size=0.4) +
  geom_point(alpha=0.9, size=2) +
  aes(x=EnrichRatio, y=plotpoint, shape=group, fill=-log10(p.adjust), alpha=0.9) +
  theme_bw() +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10), expand = c(0.2,0.2)) +
  scale_y_continuous(breaks = length(GO_list_BP_pair_cons):1, labels = GO_list_BP_pair_cons, name = "", expand = c(0.05,0.05)) +
  #scale_y_continuous(breaks = length(GO_list_BP_pair_cons):1, labels = GO_list_BP_pair_cons, name = "",pair.axis = dup_axis(labels = df_class_total_BP$class_gene_total, name = "Total gene counts")) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  scale_shape_manual(values=c(21, 22,23,24,25)) +
  theme(axis.text = element_text(size=9, color='black'), axis.title = element_text(size=10, color='black'), 
        legend.title = element_text(size=9, color='black'), legend.text = element_text(size=8, color='black'), 
        legend.position = 'bottom', panel.grid = element_blank(),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.05, 'in')) +
  guides(shape= guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3)), 
         size= guide_legend(nrow=1, order = 2, title.position = "top"), 
         fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2)) +
  labs(x="Fold enrichment", shape = "Comparison", fill = "-log10(Corrected P-value)") 

select_overlap_Des <- c("cell-cell adhesion via plasma-membrane adhesion molecules",
                        "regulation of membrane potential",
                        "adenylate cyclase-activating G protein-coupled receptor signaling pathway",
                        "motor neuron axon guidance",
                        "compound eye photoreceptor development",
                        "potassium ion transmembrane transport",
                        "modulation of chemical synaptic transmission",
                        "genital disc development" ,
                        "mating behavior",
                        "cell fate specification",
                        "regulation of neurotransmitter levels",
                        "neuromuscular junction development",
                        "neurotransmitter transport",
                        "organic anion transport",
                        "regulation of vesicle-mediated transport")

GO_list_BP_pair_cons_plotpoint_select <- data.frame(Description=select_overlap_Des, plotpoint=length(select_overlap_Des):1)

df_GO_enrich_BP_pair_cons_sum_select <- df_GO_enrich_BP_pair_cons %>%
  dplyr::filter(Description %in% select_overlap_Des) %>%
  dplyr::left_join(., GO_list_BP_pair_cons_plotpoint_select, by = "Description") %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(class_gene_total = max(Total_class)) %>%
  dplyr::ungroup()

plot_GO_BP_pair_cons <- df_GO_enrich_BP_pair_cons_sum_select %>%
  dplyr::filter(Description %in% select_overlap_Des) %>%
  ggplot(.) +
  geom_vline(xintercept = 1, color='blue', size=0.4) +
  geom_point(alpha=0.9, size=2) +
  aes(x=EnrichRatio, y=plotpoint, shape=group, fill=-log10(p.adjust), alpha=0.9) +
  theme_bw() +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10), expand = c(0.2,0.2)) +
  scale_y_continuous(breaks = length(select_overlap_Des):1, labels = select_overlap_Des, name = "", expand = c(0.05,0.05)) +
  #scale_y_continuous(breaks = length(GO_list_BP_pair_cons):1, labels = GO_list_BP_pair_cons, name = "",pair.axis = dup_axis(labels = df_class_total_BP$class_gene_total, name = "Total gene counts")) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  scale_shape_manual(values=c(21, 22,23,24,25)) +
  theme(axis.text = element_text(size=9, color='black'), axis.title = element_text(size=10, color='black'), 
        legend.title = element_text(size=9, color='black'), legend.text = element_text(size=8, color='black'), 
        legend.position = 'bottom', panel.grid = element_blank(),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.05, 'in')) +
  guides(shape= guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3)), 
         size= guide_legend(nrow=1, order = 2, title.position = "top"), 
         fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2)) +
  labs(x="Fold enrichment", shape = "Comparison", fill = "-log10(Corrected P-value)") 

plot_GO_BP_pair_cons

ggsave("Plots/Manuscript/Fig5_deg_cons_GO.pdf", width=8, height=4.3)


### heatmap -pheatmap ###

## extract GO enriched genes ##

topGO_BP_genes <- NULL
topGO_BP_genes_GO <- NULL
topGO_BP_genes_group <- NULL

for (i in 1:nrow(df_GO_enrich_BP_pair_cons)) {
  
  genes <- strsplit(df_GO_enrich_BP_pair_cons$geneID[i], "/")[[1]]
  group <- rep(df_GO_enrich_BP_pair_cons$group[[i]], length(genes))
  GOs <- rep(df_GO_enrich_BP_pair_cons$Description[[i]], length(genes))
  
  topGO_BP_genes <- c(topGO_BP_genes, genes)
  topGO_BP_genes_GO <- c(topGO_BP_genes_GO, GOs)
  topGO_BP_genes_group <- c(topGO_BP_genes_group, group)
  
}

df_topGO_BP_genes <- data.frame(gene=topGO_BP_genes, GO=topGO_BP_genes_GO, group=topGO_BP_genes_group) %>%
  dplyr::mutate(hit=1)

## adhesion ##

df_topGO_BP_genes_adhesion <- data.frame(gene=topGO_BP_genes, GO=topGO_BP_genes_GO, group=topGO_BP_genes_group) %>%
  dplyr::mutate(hit=1) %>%
  dplyr::filter(GO =="cell-cell adhesion via plasma-membrane adhesion molecules")

matrix_adhesion <- matrix(ncol=length(unique(df_topGO_BP_genes_adhesion$gene)), nrow=length(unique(df_topGO_BP_genes_adhesion$group)),
                      dimnames = list(unique(df_topGO_BP_genes_adhesion$group),unique(df_topGO_BP_genes_adhesion$gene)))

df_adhesion <- matrix_adhesion %>%
  as.data.frame() %>%
  tibble::rownames_to_column('group') %>%
  tidyr::gather(-group, key='gene', value='na') %>%
  dplyr::select(-na) %>%
  dplyr::left_join(., df_topGO_BP_genes_adhesion, by=c('gene', 'group')) %>%
  dplyr::mutate(hit=ifelse(is.na(hit), 0, 1)) %>%
  dplyr::select(-GO) %>%
  tidyr::spread(key='gene', value='hit')

df_adhesion$group <- factor(df_adhesion$group, 
                            levels = c("bg", "conserved", "melsec", "melsim", "simsec"), 
                            labels = c("Background", "Conserved", "Dmel vs Dsec", "Dmel vs Dsim", "Dsim vs Dsec"))

rownames(df_adhesion) <- df_adhesion$group

mt_adhesion <- df_adhesion %>%
  dplyr::select(-group) %>%
  as.matrix() 

plot_go_adhesion <- pheatmap::pheatmap(mt_adhesion, angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                                                annotation_names_col = T, show_colnames = T, cluster_rows = F,
                                                treeheight_col = 0, treeheight_row = 0,
                                                color=c("white","darkblue"),border_color = "black")

col_dend <- plot_go_adhesion[[2]]
col_dend <- dendextend::rotate(col_dend, order=c(2,1))

plot_go_adhesion_rotate <- pheatmap::pheatmap(mt_adhesion, cluster_cols=as.hclust(col_dend), 
                   angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                   annotation_names_col = 0, show_colnames = T, cluster_rows = F,
                   treeheight_col = 0, treeheight_row = 0,
                   color=c("white","darkblue"),border_color = "black")

ggsave(plot_go_adhesion_rotate, file="Plots/Manuscript/Fig5_pheatmap_GO_adhesion.pdf", width=12, height= 1.8)

## cAMP signaling ##

df_topGO_BP_genes_cAMP <- data.frame(gene=topGO_BP_genes, GO=topGO_BP_genes_GO, group=topGO_BP_genes_group) %>%
  dplyr::mutate(hit=1) %>%
  dplyr::filter(GO =="adenylate cyclase-activating G protein-coupled receptor signaling pathway")

matrix_cAMP <- matrix(ncol=length(unique(df_topGO_BP_genes_cAMP$gene)), nrow=length(unique(df_topGO_BP_genes_cAMP$group)),
                      dimnames = list(unique(df_topGO_BP_genes_cAMP$group),unique(df_topGO_BP_genes_cAMP$gene)))

df_cAMP <- matrix_cAMP %>%
  as.data.frame() %>%
  tibble::rownames_to_column('group') %>%
  tidyr::gather(-group, key='gene', value='na') %>%
  dplyr::select(-na) %>%
  dplyr::left_join(., df_topGO_BP_genes_cAMP, by=c('gene', 'group')) %>%
  dplyr::mutate(hit=ifelse(is.na(hit), 0, 1)) %>%
  dplyr::select(-GO) %>%
  tidyr::spread(key='gene', value='hit')

df_cAMP$group <- factor(df_cAMP$group, 
                        levels = c("bg", "conserved", "melsec", "melsim", "simsec"), 
                        labels = c("Background (3000)", "Conserved (173)", "Dmel vs Dsec (487)", "Dmel vs Dsim (369)", "Dsim vs Dsec (242)"))

rownames(df_cAMP) <- df_cAMP$group

mt_cAMP <- df_cAMP %>%
  dplyr::select(-group) %>%
  as.matrix() 

plot_go_cAMP <- pheatmap::pheatmap(mt_cAMP, angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                                   annotation_names_col = T, show_colnames = T, cluster_rows = F,
                                   treeheight_col = 0, treeheight_row = 0,
                                   color=c("white","darkblue"),border_color = "black")

col_dend <- plot_go_cAMP[[2]]
col_dend <- dendextend::rotate(col_dend, order=c(2,1))

plot_go_cAMP_rotate <- pheatmap::pheatmap(mt_cAMP, cluster_cols=as.hclust(col_dend), 
                                              angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                                              annotation_names_col = 0, show_colnames = T, cluster_rows = F,
                                              treeheight_col = 0, treeheight_row = 0,
                                              color=c("white","darkblue"),border_color = "black")

ggsave(plot_go_cAMP, file="Plots/Manuscript/Fig5_pheatmap_GO_cAMP.pdf", width=7.5, height= 1.8)


## cell fate specification ##

df_topGO_BP_genes_cellfate <- data.frame(gene=topGO_BP_genes, GO=topGO_BP_genes_GO, group=topGO_BP_genes_group) %>%
  dplyr::mutate(hit=1) %>%
  dplyr::filter(GO =="cell fate specification")

matrix_cellfate <- matrix(ncol=length(unique(df_topGO_BP_genes_cellfate$gene)), nrow=length(unique(df_topGO_BP_genes_cellfate$group)),
                      dimnames = list(unique(df_topGO_BP_genes_cellfate$group),unique(df_topGO_BP_genes_cellfate$gene)))

df_cellfate <- matrix_cellfate %>%
  as.data.frame() %>%
  tibble::rownames_to_column('group') %>%
  tidyr::gather(-group, key='gene', value='na') %>%
  dplyr::select(-na) %>%
  dplyr::left_join(., df_topGO_BP_genes_cellfate, by=c('gene', 'group')) %>%
  dplyr::mutate(hit=ifelse(is.na(hit), 0, 1)) %>%
  dplyr::select(-GO) %>%
  tidyr::spread(key='gene', value='hit')

df_cellfate$group <- factor(df_cellfate$group, 
                        levels = c("bg", "conserved", "melsec", "melsim", "simsec"), 
                        labels = c("Background (3000)", "Conserved (173)", "Dmel vs Dsec (487)", "Dmel vs Dsim (369)", "Dsim vs Dsec (242)"))

rownames(df_cellfate) <- df_cellfate$group

mt_cellfate <- df_cellfate %>%
  dplyr::select(-group) %>%
  as.matrix() 

plot_go_cellfate <- pheatmap::pheatmap(mt_cellfate, angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                                   annotation_names_col = T, show_colnames = T, cluster_rows = F,
                                   treeheight_col = 0, treeheight_row = 0,
                                   color=c("white","darkblue"),border_color = "black")

ggsave(plot_go_cellfate, file="Plots/Manuscript/Fig5_pheatmap_GO_cellfate.pdf", width=10.5, height= 1.8)



## membrane potential ##

df_topGO_BP_genes_memP <- data.frame(gene=topGO_BP_genes, GO=topGO_BP_genes_GO, group=topGO_BP_genes_group) %>%
  dplyr::mutate(hit=1) %>%
  dplyr::filter(GO =="regulation of membrane potential")

matrix_memP <- matrix(ncol=length(unique(df_topGO_BP_genes_memP$gene)), nrow=length(unique(df_topGO_BP_genes_memP$group)),
                          dimnames = list(unique(df_topGO_BP_genes_memP$group),unique(df_topGO_BP_genes_memP$gene)))

df_memP <- matrix_memP %>%
  as.data.frame() %>%
  tibble::rownames_to_column('group') %>%
  tidyr::gather(-group, key='gene', value='na') %>%
  dplyr::select(-na) %>%
  dplyr::left_join(., df_topGO_BP_genes_memP, by=c('gene', 'group')) %>%
  dplyr::mutate(hit=ifelse(is.na(hit), 0, 1)) %>%
  dplyr::select(-GO) %>%
  tidyr::spread(key='gene', value='hit')

df_memP$group <- factor(df_cellfate$group, 
                        levels = c("bg", "conserved", "melsec", "melsim", "simsec"), 
                        labels = c("Background (3000)", "Conserved (173)", "Dmel vs Dsec (487)", "Dmel vs Dsim (369)", "Dsim vs Dsec (242)"))

rownames(df_memP) <- df_memP$group

mt_memP <- df_memP %>%
  dplyr::select(-group) %>%
  as.matrix() 

plot_go_memP <- pheatmap::pheatmap(mt_memP, angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                                       annotation_names_col = T, show_colnames = T, cluster_rows = F,
                                       treeheight_col = 0, treeheight_row = 0,
                                       color=c("white","darkblue"),border_color = "black")

ggsave(plot_go_memP, file="Plots/Manuscript/Fig5_pheatmap_GO_memP.pdf", width=12, height= 1.8)


## synaptic transmission ##

df_topGO_BP_genes_synT <- data.frame(gene=topGO_BP_genes, GO=topGO_BP_genes_GO, group=topGO_BP_genes_group) %>%
  dplyr::mutate(hit=1) %>%
  dplyr::filter(GO =="modulation of chemical synaptic transmission")

matrix_synT <- matrix(ncol=length(unique(df_topGO_BP_genes_synT$gene)), nrow=length(unique(df_topGO_BP_genes_synT$group)),
                      dimnames = list(unique(df_topGO_BP_genes_synT$group),unique(df_topGO_BP_genes_synT$gene)))

df_synT <- matrix_synT %>%
  as.data.frame() %>%
  tibble::rownames_to_column('group') %>%
  tidyr::gather(-group, key='gene', value='na') %>%
  dplyr::select(-na) %>%
  dplyr::left_join(., df_topGO_BP_genes_synT, by=c('gene', 'group')) %>%
  dplyr::mutate(hit=ifelse(is.na(hit), 0, 1)) %>%
  dplyr::select(-GO) %>%
  tidyr::spread(key='gene', value='hit')

df_synT$group <- factor(df_synT$group, 
                        levels = c("bg", "conserved", "melsec", "melsim", "simsec"), 
                        labels = c("Background (3000)", "Conserved (173)", "Dmel vs Dsec (487)", "Dmel vs Dsim (369)", "Dsim vs Dsec (242)"))

rownames(df_synT) <- df_synT$group

mt_synT <- df_synT %>%
  dplyr::select(-group) %>%
  as.matrix() 

plot_go_synT <- pheatmap::pheatmap(mt_synT, angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                                   annotation_names_col = T, show_colnames = T, cluster_rows = F,
                                   treeheight_col = 0, treeheight_row = 0,
                                   color=c("white","darkblue"),border_color = "black")

col_dend <- plot_go_synT[[2]]
col_dend <- dendextend::rotate(col_dend, order=c(2,4))

plot_go_synT_rotate <- pheatmap::pheatmap(mt_synT, cluster_cols=as.hclust(col_dend), 
                                          angle_col = 45, border = NA, annotation_names_row = T, show_rownames = T, 
                                          annotation_names_col = 0, show_colnames = T, cluster_rows = F,
                                          treeheight_col = 0, treeheight_row = 0,
                                          color=c("white","darkblue"),border_color = "black")

ggsave(plot_go_synT_rotate, file="Plots/Manuscript/Fig5_pheatmap_GO_synT.pdf", width=17, height= 1.8)



## merge all - GO term plot ##

## MF ##

bg_go_MF_select <- bg_go_MF@result %>%
  dplyr::select(ID, Description, GeneRatio_bg=GeneRatio,BgRatio, p.adjust.bg = p.adjust) %>%
  dplyr::mutate(group="bg")

ego_MF_melsec_select <- deg_melsec_go_MF@result %>%
  dplyr::select(ID, Description, GeneRatio_target=GeneRatio,BgRatio, p.adjust.target = p.adjust) %>%
  dplyr::mutate(group="melsec")

ego_MF_melsim_select <- deg_melsim_go_MF@result %>%
  dplyr::select(ID, Description, GeneRatio_target=GeneRatio,BgRatio, p.adjust.target = p.adjust) %>%
  dplyr::mutate(group="melsim")

ego_MF_simsec_select <- deg_simsec_go_MF@result %>%
  dplyr::select(ID, Description, GeneRatio_target=GeneRatio,BgRatio, p.adjust.target = p.adjust) %>%
  dplyr::mutate(group="simsec")

ego_MF_sim_select <- Dsim_specific_MF@result %>%
  dplyr::select(ID, Description, GeneRatio_target=GeneRatio,BgRatio, p.adjust.target = p.adjust) %>%
  dplyr::mutate(group="Dsim")

ego_MF_sec_select <- Dsec_specific_MF@result %>%
  dplyr::select(ID, Description, GeneRatio_target=GeneRatio,BgRatio, p.adjust.target = p.adjust) %>%
  dplyr::mutate(group="Dsec")

ego_cons_rho_select <- ego_conser_rho_MF@result %>%
  dplyr::select(ID, Description, GeneRatio_target=GeneRatio,BgRatio, p.adjust.target = p.adjust) %>%
  dplyr::mutate(group="conserved(rho)")

ego_MF_pair_cons_select_merge <- rbind(ego_MF_melsec_select, ego_MF_melsim_select, ego_MF_simsec_select, 
                                       ego_cons_rho_select)

go_MF_comp_pair_cons <- ego_MF_pair_cons_select_merge %>%
  dplyr::left_join(., bg_go_MF_select, by=c('ID','Description','BgRatio')) %>%
  tidyr::separate(col=GeneRatio_target, into=c("A1","A2"),sep="/") %>%
  dplyr::mutate(GeneRatio_target=as.numeric(A1)/as.numeric(A2))%>%
  tidyr::separate(col=GeneRatio_bg, into=c("B1","B2"),sep="/") %>%
  dplyr::mutate(GeneRatio_bg=as.numeric(B1)/as.numeric(B2)) %>%
  dplyr::mutate(ratio_diff = GeneRatio_target/GeneRatio_bg)

go_MF_comp_pair_cons_filtered <- go_MF_comp_pair_cons %>%
  dplyr::filter(p.adjust.target < 0.01 & as.numeric(A1) >= 10) %>%
  dplyr::group_by(group.x) %>%
  dplyr::top_n(n=10, wt=ratio_diff)

df_GO_enrich_MF_pair_cons <- rbind(dplyr::mutate(deg_melsec_go_MF@result, group="melsec"), 
                                   dplyr::mutate(deg_melsim_go_MF@result, group="melsim"), 
                                   dplyr::mutate(deg_simsec_go_MF@result, group="simsec"), 
                                   dplyr::mutate(ego_conser_rho_MF@result, group="conserved"),
                                   dplyr::mutate(bg_go_MF@result, group="bg")) %>%
  tidyr::separate(GeneRatio, c("N_div_class", "N_div_all")) %>%
  tidyr::separate(BgRatio, c("Total_class", "Total_genome")) %>%
  dplyr::mutate(N_div_class = as.numeric(N_div_class), N_div_all = as.numeric(N_div_all), Total_class = as.numeric(Total_class), Total_genome = as.numeric(Total_genome)) %>%
  dplyr::mutate(GeneRatio = N_div_class/N_div_all, BgRatio = Total_class/Total_genome) %>%
  dplyr::mutate(EnrichRatio = GeneRatio/BgRatio) %>%
  dplyr::arrange(-EnrichRatio) %>%
  dplyr::filter(ID %in% go_MF_comp_pair_cons_filtered$ID)

GO_list_MF_pair_cons <- unique(df_GO_enrich_MF_pair_cons$Description)

GO_list_MF_pair_cons_plotpoint <- data.frame(Description=GO_list_MF_pair_cons, plotpoint=length(GO_list_MF_pair_cons):1)

df_GO_enrich_MF_pair_cons_sum <- df_GO_enrich_MF_pair_cons %>%
  dplyr::filter(Description %in% GO_list_MF_pair_cons) %>%
  dplyr::left_join(., GO_list_MF_pair_cons_plotpoint, by = "Description") %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(class_gene_total = max(Total_class)) %>%
  dplyr::ungroup()

df_GO_enrich_MF_pair_cons_sum %>%
  ggplot(.) +
  geom_vline(xintercept = 1, color='blue', size=0.4) +
  geom_point(alpha=0.9, size=2) +
  aes(x=EnrichRatio, y=plotpoint, shape=group, fill=-log10(p.adjust), alpha=0.9) +
  theme_bw() +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10), expand = c(0.2,0.2)) +
  scale_y_continuous(breaks = length(GO_list_MF_pair_cons):1, labels = GO_list_MF_pair_cons, name = "", expand = c(0.05,0.05)) +
  #scale_y_continuous(breaks = length(GO_list_MF_pair_cons):1, labels = GO_list_MF_pair_cons, name = "",pair.axis = dup_axis(labels = df_class_total_MF$class_gene_total, name = "Total gene counts")) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  scale_shape_manual(values=c(21, 22,23,24,25)) +
  theme(axis.text = element_text(size=9, color='black'), axis.title = element_text(size=10, color='black'), 
        legend.title = element_text(size=9, color='black'), legend.text = element_text(size=8, color='black'), 
        legend.position = 'bottom', panel.grid = element_blank(),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.05, 'in')) +
  guides(shape= guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3)), 
         size= guide_legend(nrow=1, order = 2, title.position = "top"), 
         fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2)) +
  labs(x="Fold enrichment", shape = "Comparison", fill = "-log10(Corrected P-value)") 

select_overlap_Des <- c("cell-cell adhesion via plasma-membrane adhesion molecules",
                        "regulation of membrane potential",
                        "adenylate cyclase-activating G protein-coupled receptor signaling pathway",
                        "motor neuron axon guidance",
                        "compound eye photoreceptor development",
                        "potassium ion transmembrane transport",
                        "modulation of chemical synaptic transmission",
                        "genital disc development" ,
                        "mating behavior",
                        "cell fate specification",
                        "regulation of neurotransmitter levels",
                        "neuromuscular junction development",
                        "neurotransmitter transport",
                        "organic anion transport",
                        "regulation of vesicle-mediated transport")

GO_list_MF_pair_cons_plotpoint_select <- data.frame(Description=select_overlap_Des, plotpoint=length(select_overlap_Des):1)

df_GO_enrich_MF_pair_cons_sum_select <- df_GO_enrich_MF_pair_cons %>%
  dplyr::filter(Description %in% select_overlap_Des) %>%
  dplyr::left_join(., GO_list_MF_pair_cons_plotpoint_select, by = "Description") %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(class_gene_total = max(Total_class)) %>%
  dplyr::ungroup()

plot_GO_MF_pair_cons <- df_GO_enrich_MF_pair_cons_sum_select %>%
  dplyr::filter(Description %in% select_overlap_Des) %>%
  ggplot(.) +
  geom_vline(xintercept = 1, color='blue', size=0.4) +
  geom_point(alpha=0.9, size=2) +
  aes(x=EnrichRatio, y=plotpoint, shape=group, fill=-log10(p.adjust), alpha=0.9) +
  theme_bw() +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10), expand = c(0.2,0.2)) +
  scale_y_continuous(breaks = length(select_overlap_Des):1, labels = select_overlap_Des, name = "", expand = c(0.05,0.05)) +
  #scale_y_continuous(breaks = length(GO_list_MF_pair_cons):1, labels = GO_list_MF_pair_cons, name = "",pair.axis = dup_axis(labels = df_class_total_MF$class_gene_total, name = "Total gene counts")) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  scale_shape_manual(values=c(21, 22,23,24,25)) +
  theme(axis.text = element_text(size=9, color='black'), axis.title = element_text(size=10, color='black'), 
        legend.title = element_text(size=9, color='black'), legend.text = element_text(size=8, color='black'), 
        legend.position = 'bottom', panel.grid = element_blank(),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.05, 'in')) +
  guides(shape= guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3)), 
         size= guide_legend(nrow=1, order = 2, title.position = "top"), 
         fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2)) +
  labs(x="Fold enrichment", shape = "Comparison", fill = "-log10(Corrected P-value)") 

plot_GO_MF_pair_cons

ggsave("Plots/Manuscript/Fig5_deg_cons_GO.pdf", width=8, height=4.3)

