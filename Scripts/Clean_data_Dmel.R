library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(DoubletFinder)
library(SoupX)
library(glmGamPoi)
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

### Decontam with SoupX ###

sc_Dmel1 <- load10X('Dmel/S1_full_intron/outs/')
sc_Dmel1 = autoEstCont(sc_Dmel1)
out_Dmel1 = adjustCounts(sc_Dmel1)
Dmel_rep1_decon <- CreateSeuratObject(counts = out_Dmel1, project = "Dmel_rep1")

sc_Dmel2 <- load10X('Dmel/S2_full_intron/outs/')
sc_Dmel2 = autoEstCont(sc_Dmel2)
out_Dmel2 = adjustCounts(sc_Dmel2)
Dmel_rep2_decon <- CreateSeuratObject(counts = out_Dmel2, project = "Dmel_rep2")

sc_Dmel3 <- load10X('Dmel/S3_full_intron/outs/')
sc_Dmel3 = autoEstCont(sc_Dmel3)
out_Dmel3 = adjustCounts(sc_Dmel3)
Dmel_rep3_decon <- CreateSeuratObject(counts = out_Dmel3, project = "Dmel_rep3")

sc_Dmel4 <- load10X('Dmel/S4_full_intron/outs/')
sc_Dmel4 = autoEstCont(sc_Dmel4)
out_Dmel4 = adjustCounts(sc_Dmel4)
Dmel_rep4_decon <- CreateSeuratObject(counts = out_Dmel4, project = "Dmel_rep4")

sc_Dmel5 <- load10X('Dmel/S5_full_intron/outs/')
sc_Dmel5 = autoEstCont(sc_Dmel5)
out_Dmel5 = adjustCounts(sc_Dmel5)
Dmel_rep5_decon <- CreateSeuratObject(counts = out_Dmel5, project = "Dmel_rep5")

sc_Dmel6 <- load10X('Dmel/S6_full_intron/outs/')
sc_Dmel6 = autoEstCont(sc_Dmel6)
out_Dmel6 = adjustCounts(sc_Dmel6)
Dmel_rep6_decon <- CreateSeuratObject(counts = out_Dmel6, project = "Dmel_rep6")

# Merge dataset
Dmel_decon_combined12<- merge(Dmel_rep1_decon, y = Dmel_rep2_decon, add.cell.ids = c("rep1", "rep2"), project = "Dmel_full_intron_decon")
Dmel_decon_combined34<- merge(Dmel_rep3_decon, y = Dmel_rep4_decon, add.cell.ids = c("rep3", "rep4"), project = "Dmel_full_intron_decon")
Dmel_decon_combined1234<- merge(Dmel_decon_combined12, y = Dmel_decon_combined34, project = "Dmel_full_intron_decon")
Dmel_decon_combined56<- merge(Dmel_rep5_decon, y = Dmel_rep6_decon, add.cell.ids = c("rep5", "rep6"), project = "Dmel_full_intron_decon")
Dmel_decon_combined<- merge(Dmel_decon_combined1234, y = Dmel_decon_combined56, project = "Dmel_full_intron_decon")

#saveRDS(Dmel_decon_combined, file = "Processed_Data/Dmel_combined_full_intron_decon_unfiltered.rds")
#Dmel_decon_combined <- readRDS(file = "Processed_Data/Dmel_combined_full_intron_decon_unfiltered.rds")

### QC check ###

Dmel_decon_combined[["percent.mt"]] <- PercentageFeatureSet(Dmel_decon_combined, pattern = "^mt:")

VlnPlot(Dmel_decon_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) ## mostly clean 
FeatureScatter(Dmel_decon_combined, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.2)
FeatureScatter(Dmel_decon_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1)

### SCT / clustering ###

## SCTransform

Dmel_decon_combined <- SCTransform(Dmel_decon_combined, method = "glmGamPoi", verbose = FALSE)

### Perform linear dimensional reduction

Dmel_decon_combined <- RunPCA(Dmel_decon_combined, features = VariableFeatures(object = Dmel_decon_combined))

### Clustering

Dmel_decon_combined <- RunUMAP(Dmel_decon_combined, dims = 1:50)
Dmel_decon_combined <- RunTSNE(Dmel_decon_combined, dims = 1:50)
Dmel_decon_combined <- FindNeighbors(Dmel_decon_combined, dims = 1:50)
Dmel_decon_combined <- FindClusters(Dmel_decon_combined, resolution = 0.8)

DimPlot(Dmel_decon_combined, reduction = "umap",  label = TRUE, pt.size = 0.01)

#saveRDS(Dmel_decon_combined, file = "Processed_Data/Dmel_combined_full_intron_decon.rds")
#Dmel_decon_combined <- readRDS(file = "Processed_Data/Dmel_combined_full_intron_decon.rds")


### Doublet identification ###

## personalize SCT x paramSweep_v3

parallel_paramSweep_v3_sct <- function (n, n.real.cells, real.cells, pK, pN, data, orig.commands, 
                                        PCs, sct) 
{
  sweep.res.list = list()
  list.ind = 0
  print(paste("Creating artificial doublets for pN = ", pN[n] * 
                100, "%", sep = ""))
  n_doublets <- round(n.real.cells/(1 - pN[n]) - n.real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[, real.cells1] + data[, real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  data_wdoublets <- cbind(data, doublets)
  if (sct == FALSE) {
    print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
    print("Normalizing Seurat object...")
    seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
                                   scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor, 
                                   margin = orig.commands$NormalizeData.RNA@params$margin)
    print("Finding variable genes...")
    seu_wdoublets <- FindVariableFeatures(seu_wdoublets, 
                                          selection.method = orig.commands$FindVariableFeatures.RNA$selection.method, 
                                          loess.span = orig.commands$FindVariableFeatures.RNA$loess.span, 
                                          clip.max = orig.commands$FindVariableFeatures.RNA$clip.max, 
                                          mean.function = orig.commands$FindVariableFeatures.RNA$mean.function, 
                                          dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function, 
                                          num.bin = orig.commands$FindVariableFeatures.RNA$num.bin, 
                                          binning.method = orig.commands$FindVariableFeatures.RNA$binning.method, 
                                          nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures, 
                                          mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff, 
                                          dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
    print("Scaling data...")
    seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                               model.use = orig.commands$ScaleData.RNA$model.use, 
                               do.scale = orig.commands$ScaleData.RNA$do.scale, 
                               do.center = orig.commands$ScaleData.RNA$do.center, 
                               scale.max = orig.commands$ScaleData.RNA$scale.max, 
                               block.size = orig.commands$ScaleData.RNA$block.size, 
                               min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
    print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                            npcs = length(PCs), rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                            weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, 
                            verbose = FALSE)
  }
  if (sct == TRUE) {
    require(sctransform)
    print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
    print("Running SCTransform...")
    seu_wdoublets <- SCTransform(seu_wdoublets, method = "glmGamPoi")
    print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
  }
  print("Calculating PC distance matrix...")
  nCells <- nrow(seu_wdoublets@meta.data)
  pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                            PCs]
  rm(seu_wdoublets)
  gc()
  dist.mat <- fields::rdist(pca.coord)[, 1:n.real.cells]
  print("Defining neighborhoods...")
  for (i in 1:n.real.cells) {
    dist.mat[, i] <- order(dist.mat[, i])
  }
  ind <- round(nCells * max(pK)) + 5
  dist.mat <- dist.mat[1:ind, ]
  print("Computing pANN across all pK...")
  for (k in 1:length(pK)) {
    print(paste("pK = ", pK[k], "...", sep = ""))
    pk.temp <- round(nCells * pK[k])
    pANN <- as.data.frame(matrix(0L, nrow = n.real.cells, 
                                 ncol = 1))
    colnames(pANN) <- "pANN"
    rownames(pANN) <- real.cells
    list.ind <- list.ind + 1
    for (i in 1:n.real.cells) {
      neighbors <- dist.mat[2:(pk.temp + 1), i]
      pANN$pANN[i] <- length(which(neighbors > n.real.cells))/pk.temp
    }
    sweep.res.list[[list.ind]] <- pANN
  }
  return(sweep.res.list)
} ### SCT with method = "glmGamPoi"

paramSweep_v3_SCT <- function (seu, PCs = 1:10, sct = T, num.cores = 1) 
{
  require(Seurat)
  require(fields)
  pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)
  min.cells <- round(nrow(seu@meta.data)/(1 - 0.05) - nrow(seu@meta.data))
  pK.test <- round(pK * min.cells)
  pK <- pK[which(pK.test >= 1)]
  orig.commands <- seu@commands
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 
                                                 10000, replace = FALSE)]
    data <- seu@assays$RNA@counts[, real.cells]
    n.real.cells <- ncol(data)
  }
  if (nrow(seu@meta.data) <= 10000) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }
  if (num.cores > 1) {
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)), FUN = parallel_paramSweep_v3_sct, 
                        n.real.cells, real.cells, pK, pN, data, orig.commands, 
                        PCs, sct, mc.cores = num.cores)
    stopCluster(cl)
  }
  else {
    output2 <- lapply(as.list(1:length(pN)), FUN = parallel_paramSweep_v3_sct, 
                      n.real.cells, real.cells, pK, pN, data, orig.commands, 
                      PCs, sct)
  }
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2)) {
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, 
                                  sep = "_"))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}

doubletFinder_v3_SCT <- function (seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, 
                                  sct = T) 
{
  require(Seurat)
  require(fields)
  require(KernSmooth)
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    return(seu)
  }
  if (reuse.pANN == FALSE) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...", 
                sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    orig.commands <- seu@commands
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor, 
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets, 
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method, 
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span, 
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max, 
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function, 
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function, 
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin, 
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method, 
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures, 
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff, 
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                                 model.use = orig.commands$ScaleData.RNA$model.use, 
                                 do.scale = orig.commands$ScaleData.RNA$do.scale, 
                                 do.center = orig.commands$ScaleData.RNA$do.center, 
                                 scale.max = orig.commands$ScaleData.RNA$scale.max, 
                                 block.size = orig.commands$ScaleData.RNA$block.size, 
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                              npcs = length(PCs), rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, 
                              verbose = FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets, method = "glmGamPoi")
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                 ncol = 1))
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      neighbor.names <- rownames(dist.mat)[neighbors]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 
                                                                    1]
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    return(seu)
  }
}

## pK Identification (no ground-truth)
sweep.res.list_Dmel_decon_combined <- paramSweep_v3_SCT(Dmel_decon_combined, PCs = 1:50, sct = T)
sweep.stats_Dmel_decon_combined <- summarizeSweep(sweep.res.list_Dmel_decon_combined, GT = FALSE)
bcmvn_Dmel_decon_combined <- find.pK(sweep.stats_Dmel_decon_combined)  ## pk=0.12

## Homotypic Doublet Proportion Estimate : 15%
annotations <- Dmel_decon_combined@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.15*nrow(Dmel_decon_combined@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies

Dmel_decon_combined <- doubletFinder_v3_SCT(Dmel_decon_combined, PCs = 1:50, pN = 0.25, pK = 0.12, nExp = nExp_poi.adj, reuse.pANN = F, sct = T)

DimPlot(Dmel_decon_combined,
        pt.size = 0.01, label=F, reduction = "umap",
        group.by = "DF.classifications_0.25_0.12_8043" )+theme(aspect.ratio = 1)
 
saveRDS(Dmel_decon_combined, file = "Processed_Data/Dmel_combined_full_intron_decon.rds")
Dmel_decon_combined <- readRDS(file = "Processed_Data/Dmel_combined_full_intron_decon.rds")

## reclustering after doublet removal ##

Dmel_decon_combined_DF <- subset(Dmel_decon_combined, subset = DF.classifications_0.25_0.12_8043 == "Singlet")

Dmel_decon_combined_DF <- SCTransform(Dmel_decon_combined_DF, method = "glmGamPoi", verbose = FALSE)
Dmel_decon_combined_DF <- RunPCA(Dmel_decon_combined_DF, npcs = 50, verbose = FALSE)

Dmel_decon_combined_DF <- FindNeighbors(Dmel_decon_combined_DF, dims = 1:50)
Dmel_decon_combined_DF <- FindClusters(Dmel_decon_combined_DF, resolution = 0.8)

Dmel_decon_combined_DF <- RunUMAP(Dmel_decon_combined_DF, dims = 1:50)

DimPlot(Dmel_decon_combined_DF, reduction = "umap",  label = TRUE, pt.size = 0.01)
FeaturePlot(Dmel_decon_combined_DF, feature='Pka-R2', order=T)

VlnPlot(Dmel_decon_combined_DF, features=c('nCount_RNA', 'nFeature_RNA'), pt.size=0, ncol=1) + NoLegend()

saveRDS(Dmel_decon_combined_DF, file = "Processed_Data/Dmel_combined_full_intron_decon_DF.rds")
Dmel_decon_combined_DF <- readRDS(file = "Processed_Data/Dmel_combined_full_intron_decon_DF.rds")
  
###########  

## Subset KC and check doublet ##

DimPlot(Dmel_decon_combined_doublet, label=T)

Dmel_celltype_KC_doublet <- subset(Dmel_decon_combined_doublet, idents = c(8,9,14))

Dmel_celltype_KC_doublet <- SCTransform(Dmel_celltype_KC_doublet, method = "glmGamPoi", verbose = FALSE)

Dmel_celltype_KC_doublet <- RunPCA(Dmel_celltype_KC_doublet, npcs = 50, verbose = FALSE)
# t-SNE and Clustering
Dmel_celltype_KC_doublet <- RunUMAP(Dmel_celltype_KC_doublet, reduction = "pca", dims = 1:50)
Dmel_celltype_KC_doublet <- FindNeighbors(Dmel_celltype_KC_doublet, reduction = "pca", dims = 1:50)
Dmel_celltype_KC_doublet <- FindClusters(Dmel_celltype_KC_doublet, resolution = 0.1)

DimPlot(object=Dmel_celltype_KC_doublet, reduction="umap")

DimPlot(Dmel_celltype_KC_doublet,
        pt.size = 0.3, label=F, reduction = "umap",
        group.by = "DF.classifications_0.25_0.12_8043" )+theme(aspect.ratio = 1)

FeaturePlot(object=Dmel_celltype_KC_doublet, feature="rn")

### re-clustering after doublet filtering ###

Dmel_celltype_KC_DF <- subset(Dmel_celltype_KC_doublet, subset = DF.classifications_0.25_0.12_8043 == "Singlet")

Dmel_celltype_KC_DF <- SCTransform(Dmel_celltype_KC_DF, method = "glmGamPoi", verbose = FALSE)
Dmel_celltype_KC_DF <- RunPCA(Dmel_celltype_KC_DF, npcs = 50, verbose = FALSE)

Dmel_celltype_KC_DF <- FindNeighbors(Dmel_celltype_KC_DF, dims = 1:50)
Dmel_celltype_KC_DF <- FindClusters(Dmel_celltype_KC_DF, resolution = 0.15)

Dmel_celltype_KC_DF <- RunUMAP(Dmel_celltype_KC_DF, dims = 1:50)
DimPlot(Dmel_celltype_KC_DF, reduction = "umap",  label = TRUE, pt.size = 0.3)

## frequencies ##

Dmel_celltype_KC_DF_metadata <- Dmel_celltype_KC_DF@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))

Dmel_celltype_KC_DF_metadata_summary_rep <- Dmel_celltype_KC_DF_metadata %>%
  dplyr::group_by(species, seurat_clusters, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

Dmel_celltype_KC_DF_metadata_summary <- Dmel_celltype_KC_DF_metadata_summary_rep %>%
  dplyr::group_by(species, seurat_clusters) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

Dmel_celltype_KC_DF_metadata_summary %>%
  ggplot(.) +
  geom_col(aes(x=seurat_clusters, y=percent_combined, fill=species), alpha = 1, position = position_dodge(width=0.9)) + 
  geom_errorbar(aes(x=seurat_clusters, y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,0,percent_combined-sem), 
                    ymax=percent_combined+sem, group=species), position = position_dodge(width=0.9), width=0.7) +
  geom_point(data=Dmel_celltype_KC_DF_metadata_summary_rep, size = 1, alpha =0.8,  
             aes(x=seurat_clusters, y=percent, group=species), fill='red', shape=21, position = position_dodge(width=0.9)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=seurat_clusters, y=n/total*100)) +
  theme_bw() +
  theme(axis.title.y = element_text(size=11, color='black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10, color='black'),
        axis.text.x=element_text(size=10, angle=90, hjust = 1, vjust = 0.5, color='black'),
        panel.grid = element_blank()) +
  labs(x="Cluster", y="Percent of cell type (%)", fill="") +
  scale_fill_brewer(palette = 'Set1')

## Subset glia and check doublet ##

DimPlot(Dmel_decon_combined_doublet, label=T)

FeaturePlot(Dmel_decon_combined_doublet, feature=c("axo","trol")) #ENS - 8,20
FeaturePlot(Dmel_decon_combined_doublet, feature=c("CG42342")) #AST - 10
FeaturePlot(Dmel_decon_combined_doublet, feature=c("Tret1-1")) #PRN - 14
FeaturePlot(Dmel_decon_combined_doublet, feature="CG40470") #CTX/CHI - 28
FeaturePlot(Dmel_decon_combined_doublet, feature=c("baz", "wat","fus"), order=T) #SUB - 36

Dmel_celltype_glia_doublet <- subset(Dmel_decon_combined_doublet, idents = c(8,10,14,20,28,36)) # 

Dmel_celltype_glia_doublet <- SCTransform(Dmel_celltype_glia_doublet, method = "glmGamPoi", verbose = FALSE)

Dmel_celltype_glia_doublet <- RunPCA(Dmel_celltype_glia_doublet, npcs = 50, verbose = FALSE)
# t-SNE and Clustering
Dmel_celltype_glia_doublet <- RunUMAP(Dmel_celltype_glia_doublet, reduction = "pca", dims = 1:50)
Dmel_celltype_glia_doublet <- FindNeighbors(Dmel_celltype_glia_doublet, reduction = "pca", dims = 1:50)
Dmel_celltype_glia_doublet <- FindClusters(Dmel_celltype_glia_doublet, resolution = 0.1)

### Cell type classification using module score ###

df_KD_cluster_markers <- read.csv(file="Processed_Data/KD(2018)_cluster_markers.csv")
df_KD_cluster_annotation <- read.csv(file="Processed_Data/KD(2018)_cluster_annotation.csv") %>%
  dplyr::rename(cluster=Cluster_ID)

df_KD_cluster <- df_KD_cluster_markers %>%
  dplyr::left_join(., df_KD_cluster_annotation, by='cluster') %>%
  na.omit()

df_KD_cluster_select <- df_KD_cluster %>%
  dplyr::select(gene, celltype=Annotation)

glia_markers_KD <- df_KD_cluster_markers %>%
  dplyr::left_join(., df_KD_cluster_annotation, by='cluster') %>%
  na.omit() %>%
  dplyr::filter(grepl("glia", Annotation) | grepl("Glia", Annotation) | grepl("Astrocyte", Annotation)) %>%
  #dplyr::filter(gene %in% Dmel.combined.markers_filter$gene)  %>%
  dplyr::group_by(Annotation) %>%
  top_n(n = 50, wt = -p.value) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene, celltype=Annotation)

ENS_markers <- dplyr::filter(glia_markers_KD, grepl("Ensheathing", celltype))$gene
AST_markers <- dplyr::filter(glia_markers_KD, grepl("Astrocyte", celltype))$gene
PRN_markers <- dplyr::filter(glia_markers_KD, grepl("Perineu", celltype))$gene
SUB_markers <- dplyr::filter(glia_markers_KD, grepl("Subperi", celltype))$gene
CTX_markers <- dplyr::filter(glia_markers_KD, grepl("Cortex", celltype))$gene
CHI_markers <- dplyr::filter(glia_markers_KD, grepl("Chiasm", celltype))$gene

Dmel_celltype_glia_doublet <- AddModuleScore(object = Dmel_celltype_glia_doublet, features = list(ENS_markers), name = 'ens')
Dmel_celltype_glia_doublet <- AddModuleScore(object = Dmel_celltype_glia_doublet, features = list(AST_markers), name = 'ast')
Dmel_celltype_glia_doublet <- AddModuleScore(object = Dmel_celltype_glia_doublet, features = list(PRN_markers), name = 'prn')
Dmel_celltype_glia_doublet <- AddModuleScore(object = Dmel_celltype_glia_doublet, features = list(SUB_markers), name = 'sub')
Dmel_celltype_glia_doublet <- AddModuleScore(object = Dmel_celltype_glia_doublet, features = list(CTX_markers), name = 'ctx')
Dmel_celltype_glia_doublet <- AddModuleScore(object = Dmel_celltype_glia_doublet, features = list(CHI_markers), name = 'chi')

VlnPlot(Dmel_celltype_glia_doublet, features=c("ens1","ast1", "prn1","sub1","ctx1","chi1"), pt.size=0) + NoLegend()

DimPlot(object=Dmel_celltype_glia_doublet, reduction="umap", label=T)

DimPlot(Dmel_celltype_glia_doublet,
        pt.size = 0.3, label=F, reduction = "umap",
        group.by = "DF.classifications_0.25_0.12_8043" )+theme(aspect.ratio = 1)

### re-clustering after doublet filtering ###

Dmel_celltype_glia_DF <- subset(Dmel_celltype_glia_doublet, subset = DF.classifications_0.25_0.12_8043 == "Singlet")

Dmel_celltype_glia_DF <- SCTransform(Dmel_celltype_glia_DF, method = "glmGamPoi", verbose = FALSE)
Dmel_celltype_glia_DF <- RunPCA(Dmel_celltype_glia_DF, npcs = 50, verbose = FALSE)

Dmel_celltype_glia_DF <- FindNeighbors(Dmel_celltype_glia_DF, dims = 1:50)
Dmel_celltype_glia_DF <- FindClusters(Dmel_celltype_glia_DF, resolution = 0.15)

Dmel_celltype_glia_DF <- RunUMAP(Dmel_celltype_glia_DF, dims = 1:50)
DimPlot(Dmel_celltype_glia_DF, reduction = "umap",  label = TRUE, pt.size = 0.3)

Dmel_celltype_glia_DF <- AddModuleScore(object = Dmel_celltype_glia_DF, features = list(ENS_markers), name = 'ens')
Dmel_celltype_glia_DF <- AddModuleScore(object = Dmel_celltype_glia_DF, features = list(AST_markers), name = 'ast')
Dmel_celltype_glia_DF <- AddModuleScore(object = Dmel_celltype_glia_DF, features = list(PRN_markers), name = 'prn')
Dmel_celltype_glia_DF <- AddModuleScore(object = Dmel_celltype_glia_DF, features = list(SUB_markers), name = 'sub')
Dmel_celltype_glia_DF <- AddModuleScore(object = Dmel_celltype_glia_DF, features = list(CTX_markers), name = 'ctx')
Dmel_celltype_glia_DF <- AddModuleScore(object = Dmel_celltype_glia_DF, features = list(CHI_markers), name = 'chi')

VlnPlot(Dmel_celltype_glia_DF, features=c("ens1","ast1", "prn1","sub1","ctx1","chi1"), pt.size=0) + NoLegend()

## frequencies ##

Dmel_celltype_glia_DF_metadata <- Dmel_celltype_glia_DF@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))

Dmel_celltype_glia_DF_metadata_summary_rep <- Dmel_celltype_glia_DF_metadata %>%
  dplyr::group_by(species, seurat_clusters, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

Dmel_celltype_glia_DF_metadata_summary <- Dmel_celltype_glia_DF_metadata_summary_rep %>%
  dplyr::group_by(species, seurat_clusters) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

Dmel_celltype_glia_DF_metadata_summary %>%
  ggplot(.) +
  geom_col(aes(x=seurat_clusters, y=percent_combined, fill=species), alpha = 1, position = position_dodge(width=0.9)) + 
  geom_errorbar(aes(x=seurat_clusters, y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,0,percent_combined-sem), 
                    ymax=percent_combined+sem, group=species), position = position_dodge(width=0.9), width=0.7) +
  geom_point(data=Dmel_celltype_glia_DF_metadata_summary_rep, size = 1, alpha =0.8,  
             aes(x=seurat_clusters, y=percent, group=species), fill='red', shape=21, position = position_dodge(width=0.9)) +
  #stat_summary(fun = "mean", colour = "red", size = 2, geom = "point", shape=17,  aes(x=seurat_clusters, y=n/total*100)) +
  theme_bw() +
  theme(axis.title.y = element_text(size=11, color='black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10, color='black'),
        axis.text.x=element_text(size=10, angle=90, hjust = 1, vjust = 0.5, color='black'),
        panel.grid = element_blank()) +
  labs(x="Cluster", y="Percent of cell type (%)", fill="") +
  scale_fill_brewer(palette = 'Set1')


## Subset Cluster 0 and check doublet ##

DimPlot(Dmel_decon_combined_doublet, label=T)

Dmel_celltype_Cluster0_doublet <- subset(Dmel_decon_combined_doublet, idents = 0)

Dmel_celltype_Cluster0_doublet <- SCTransform(Dmel_celltype_Cluster0_doublet, method = "glmGamPoi", verbose = FALSE)

Dmel_celltype_Cluster0_doublet <- RunPCA(Dmel_celltype_Cluster0_doublet, npcs = 50, verbose = FALSE)
# t-SNE and Clustering
Dmel_celltype_Cluster0_doublet <- RunUMAP(Dmel_celltype_Cluster0_doublet, reduction = "pca", dims = 1:50)
Dmel_celltype_Cluster0_doublet <- FindNeighbors(Dmel_celltype_Cluster0_doublet, reduction = "pca", dims = 1:50)
Dmel_celltype_Cluster0_doublet <- FindClusters(Dmel_celltype_Cluster0_doublet, resolution = 0.5)

DimPlot(object=Dmel_celltype_Cluster0_doublet, reduction="umap")

DimPlot(Dmel_celltype_Cluster0_doublet,
        pt.size = 0.3, label=F, reduction = "umap",
        group.by = "DF.classifications_0.25_0.12_8043" )+theme(aspect.ratio = 1)

VlnPlot(Dmel_celltype_Cluster0_doublet, features=c('nCount_RNA', 'nFeature_RNA'), pt.size=0, ncol=1) + NoLegend()

### re-clustering after doublet filtering ###

Dmel_celltype_Cluster0_DF <- subset(Dmel_celltype_Cluster0_doublet, subset = DF.classifications_0.25_0.12_8043 == "Singlet")

Dmel_celltype_Cluster0_DF <- SCTransform(Dmel_celltype_Cluster0_DF, method = "glmGamPoi", verbose = FALSE)
Dmel_celltype_Cluster0_DF <- RunPCA(Dmel_celltype_Cluster0_DF, npcs = 50, verbose = FALSE)

Dmel_celltype_Cluster0_DF <- FindNeighbors(Dmel_celltype_Cluster0_DF, dims = 1:50)
Dmel_celltype_Cluster0_DF <- FindClusters(Dmel_celltype_Cluster0_DF, resolution = 0.15)

Dmel_celltype_Cluster0_DF <- RunUMAP(Dmel_celltype_Cluster0_DF, dims = 1:50)
DimPlot(Dmel_celltype_Cluster0_DF, reduction = "umap",  label = TRUE, pt.size = 0.3)

FeaturePlot(Dmel_decon_combined, feature='Gad1', order=T)
VlnPlot(Dmel_celltype_Cluster0_DF, features=c('nCount_RNA', 'nFeature_RNA'), pt.size=0, ncol=1) + NoLegend()
VlnPlot(Dmel_celltype_KC_DF, features=c('nCount_RNA', 'nFeature_RNA'), pt.size=0, ncol=1) + NoLegend()





