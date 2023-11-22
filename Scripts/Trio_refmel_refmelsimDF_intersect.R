library(Seurat)
library(tidyverse)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

### load required dataset ##

load("Processed_Data/shared_genes.RData")
load("Processed_Data/celltype_order.RData")
load("Processed_Data/Decon_unfiltered_seurat.RData")

TrioBrain.integrated_slim_labeled <- readRDS(file = "Processed_Data/TrioBrain_DF.integrated_slim_labeled.rds")

## Dmel ref only ##

### change Idents ###

Idents(Dmel_decon_combined) <- "Dmel_to_Dmel"
Idents(Dsim_decon_combined) <- "Dsim_to_Dsim"
Idents(Dsim_to_Dmel_decon_combined) <- "Dsim_to_Dmel"
Idents(Dsec_decon_combined) <- "Dsec_to_Dsim"
Idents(Dsec_to_Dmel_decon_combined) <- "Dsec_to_Dmel"

Dmelref_merge_melsim <- merge(Dmel_decon_combined, Dsim_to_Dmel_decon_combined)
Dmelref_merge_melsimsec <- merge(Dmelref_merge_melsim, Dsec_to_Dmel_decon_combined)

## intersect ID ##

Integrated_ID <- Idents(TrioBrain.integrated_slim_labeled)

intersect_ID1 <- intersect(names(Idents(Dmelref_merge_melsimsec)), gsub("-1_3", "-1", names(Integrated_ID)))

melref_merge_ID2 <- Idents(Dmelref_merge_melsimsec)
names(melref_merge_ID2) <- gsub("-1", "-1_3", names(melref_merge_ID2))
names(melref_merge_ID2) <- gsub("-1_3_1", "-1_1", names(melref_merge_ID2))
names(melref_merge_ID2) <- gsub("-1_3_2", "-1_2", names(melref_merge_ID2))

intersect_ID2 <- intersect(names(melref_merge_ID2), names(Integrated_ID))

#subset(Dmelref_merge_melsimsec, cells = intersect_ID)

Trio_Dmelref_DF_labeled <- subset(Dmelref_merge_melsimsec, cells = intersect_ID1, feature=shared_genes)
Trio_Dmelref_DF_labeled <- SCTransform(Trio_Dmelref_DF_labeled, method = "glmGamPoi", verbose = T)

TrioBrain.integrated_slim_labeled_intersect <- subset(TrioBrain.integrated_slim_labeled, cells = intersect_ID2)

Idents(Trio_Dmelref_DF_labeled) <- Idents(TrioBrain.integrated_slim_labeled_intersect)

save(Trio_Dmelref_DF_labeled, file="Processed_Data/Trio_Dmelref_DF_labeled.RData")

