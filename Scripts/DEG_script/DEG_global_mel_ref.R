library(Seurat)
library(tidyverse)
library(MAST)

load("/work/FAC/FBM/CIG/rbenton/neuroflies/Daehan/FlyBrainEvo/R/ref_DEG/bin/Trio_Dmelref_DF_labeled.RData")

test <- Trio_Dmelref_DF_labeled

test@meta.data$orig.ident <- gsub("_rep1", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep2", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep3", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep4", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep5", "", test@meta.data$orig.ident)
test@meta.data$orig.ident <- gsub("_rep6", "", test@meta.data$orig.ident)

test@meta.data$orig.ident <- factor(test@meta.data$orig.ident, 
                                    levels=c("Dmel", "Dsim_to_Dmel", "Dsec_to_Dmel"))

Idents(test) <- test$orig.ident

test_degs_melsim_refmel <- FindMarkers(test, ident.1 = "Dmel", ident.2 = "Dsim_to_Dmel", verbose = T, test.use = "MAST", min.pct = 0.05, logfc.threshold = 0) %>%
  tibble::rownames_to_column('gene')
test_degs_melsec_refmel <- FindMarkers(test, ident.1 = "Dmel", ident.2 = "Dsec_to_Dmel", verbose = T, test.use = "MAST", min.pct = 0.05, logfc.threshold = 0) %>%
  tibble::rownames_to_column('gene')
test_degs_simsec_refmel <- FindMarkers(test, ident.1 = "Dsim_to_Dmel", ident.2 = "Dsec_to_Dmel", verbose = T, test.use = "MAST", min.pct = 0.05, logfc.threshold = 0) %>%
  tibble::rownames_to_column('gene')

save(test_degs_melsim_refmel,test_degs_melsec_refmel,test_degs_simsec_refmel,
     file="/work/FAC/FBM/CIG/rbenton/neuroflies/Daehan/FlyBrainEvo/R/ref_DEG/out/DEG_global_refmel.RData")




