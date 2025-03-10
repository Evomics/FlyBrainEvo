library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))


## ortholog gene name ##

### one-to-one ortholog ###

single_ortho <- read.table("orthofinder/Dmel_Dsim/primary_transcripts/Results_Dec29_1/SingleCopyOrthogroups.txt") %>%
  dplyr::rename(orthogroup=V1)

ortho <- data.table::fread("orthofinder/Dmel_Dsim/primary_transcripts/Results_Dec29_1/Orthogroups.csv", header=T)

colnames(ortho) <- c("orthogroup", "Dmel", "Dsim")

ortho_onetoone <- ortho %>%
  dplyr::filter(orthogroup %in% single_ortho$orthogroup) %>%
  tidyr::gather(-orthogroup, key="species", value="gene_ID")

df_mel_genename <- read.table(file="Dmel/Dmel_ID_genename.tsv") %>%
  dplyr::rename(gene_ID=V1, gene_name=V2, transcript_ID=V3, transcript_name=V4)

df_mel_ortho_gene <- ortho_onetoone %>%
  dplyr::filter(species == "Dmel") %>%
  dplyr::left_join(., df_mel_genename, by='gene_ID') %>%
  dplyr::distinct(gene_name, orthogroup) %>%
  na.omit()  ## 21 proteins are withdrawn

df_sim_genename <- read.table(file="Dsim/Dsim_ID_genename.tsv") %>%
  dplyr::rename(species=V1, gene_ID=V2, gene_name=V3, transcript_ID=V4, transcript_name=V5) %>%
  dplyr::distinct(gene_ID, gene_name)

df_sim_ortho_gene <- ortho_onetoone %>%
  dplyr::filter(species == "Dsim") %>%
  dplyr::left_join(., df_sim_genename, by='gene_ID') %>%
  dplyr::distinct(gene_ID, orthogroup) %>%
  dplyr::left_join(., df_mel_ortho_gene, by='orthogroup')%>%
  na.omit()

#system("gunzip Dsim/S1_re/outs/filtered_feature_bc_matrix/archive/features.tsv.gz")
feature <- read_tsv(file="Dsim/S1_re/outs/filtered_feature_bc_matrix/archive/features.tsv", col_names=F)

### one-to-one orthologs from the blast ###

seq_id <- data.table::fread("orthofinder/Dmel_Dsim/primary_transcripts/Results_Dec29_1/WorkingDirectory/SequenceIDs.txt", header=F)
seq_id$V1 <- gsub(":","",seq_id$V1)

df_blast_melsim <- data.table::fread("orthofinder/Dmel_Dsim/primary_transcripts/Results_Dec29_1/WorkingDirectory/Blast0_1.txt", header=F) %>%
  dplyr::rename(qseqid=V1, sseqid=V2, pident=V3, length=V4, mismatch=V5, gapopen=V6, qstart=V7,
                qend=V8, sstart=V9, send=V10, evalue=V11, bitscore=V12) %>%
  dplyr::left_join(., dplyr::rename(seq_id, qseqid=V1, qseq_geneid=V2), by="qseqid") %>%
  dplyr::left_join(., dplyr::rename(seq_id, sseqid=V1, sseq_geneid=V2), by="sseqid") %>%
  dplyr::mutate(class="DmeltoDsim", pair = paste(qseq_geneid,sseq_geneid, sep="_"))

df_blast_melsim_top <- df_blast_melsim %>%
  dplyr::group_by(qseqid) %>%
  top_n(wt=bitscore, n=1) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(., dplyr::select(df_sim_genename, sseq_geneid=gene_ID, gene_name), by="sseq_geneid") %>%
  dplyr::distinct(sseqid, qseqid, pident, bitscore, qend, .keep_all=T)

df_blast_simmel <- data.table::fread("orthofinder/Dmel_Dsim/primary_transcripts/Results_Dec29_1/WorkingDirectory/Blast1_0.txt", header=F) %>%
  dplyr::rename(qseqid=V1, sseqid=V2, pident=V3, length=V4, mismatch=V5, gapopen=V6, qstart=V7,
                qend=V8, sstart=V9, send=V10, evalue=V11, bitscore=V12) %>%
  dplyr::left_join(., dplyr::rename(seq_id, qseqid=V1, qseq_geneid=V2), by="qseqid") %>%
  dplyr::left_join(., dplyr::rename(seq_id, sseqid=V1, sseq_geneid=V2), by="sseqid") %>%
  dplyr::mutate(class="DsimtoDmel", pair = paste(sseq_geneid, qseq_geneid,sep="_"))

df_blast_simmel_top <- df_blast_simmel %>%
  dplyr::group_by(qseqid) %>%
  top_n(wt=bitscore, n=1) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(., dplyr::select(df_mel_genename, sseq_geneid=gene_ID, gene_name), by="sseq_geneid") %>%
  dplyr::distinct(sseqid, qseqid, pident, bitscore, qend, .keep_all=T)

df_blast_top <- rbind(df_blast_melsim_top, df_blast_simmel_top) %>%
  dplyr::group_by(pair) %>%
  dplyr::mutate(recip = n_distinct(class)) %>%
  dplyr::ungroup()

df_blast_simmel_top_onetoone <- df_blast_top %>%
  dplyr::group_by(qseq_geneid) %>%
  dplyr::mutate(n1=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sseq_geneid) %>%
  dplyr::mutate(n2=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(recip==2, n1==1, n2==1, class=="DsimtoDmel")  

save(df_blast_simmel_top_onetoone, file="Processed_Data/df_blast_simmel_top_onetoone.RData")

feature_gene <- feature %>%
  dplyr::rename(gene_ID=X1) %>%
  dplyr::left_join(., dplyr::select(df_sim_genename, gene_ID, gene_name_sim = gene_name), by="gene_ID") %>%
  dplyr::distinct(gene_ID, .keep_all=T) %>%
  dplyr::left_join(., df_sim_ortho_gene, by="gene_ID") %>%
  dplyr::left_join(., dplyr::select(df_blast_simmel_top_onetoone, gene_ID=qseq_geneid, gene_name_tophit=gene_name, bitscore, pident, length), by='gene_ID') %>%
  dplyr::mutate(gene=ifelse(is.na(orthogroup) & bitscore > 10, gene_name_tophit, gene_name)) %>%
  dplyr::mutate(gene=ifelse(is.na(gene), gene_name_sim, gene)) %>%
  dplyr::mutate(gene=ifelse(gene_ID=="FBgn0192805", "GD21369",
                            ifelse(gene_ID=="FBgn0269846", "GD28556", gene))) %>% ## fix duplicating gene names
  dplyr::select(gene_ID,gene,X3) %>%
  na.omit()
  


  
