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

sc_Dsec1 <- load10X('Dsec_to_Dsim/S1_full_intron/outs/')
sc_Dsec1 = autoEstCont(sc_Dsec1)
out_Dsec1 = adjustCounts(sc_Dsec1)
Dsec_rep1_decon <- CreateSeuratObject(counts = out_Dsec1, project = "Dsec_rep1")

sc_Dsec2 <- load10X('Dsec_to_Dsim/S2_full_intron/outs/')
sc_Dsec2 = autoEstCont(sc_Dsec2)
out_Dsec2 = adjustCounts(sc_Dsec2)
Dsec_rep2_decon <- CreateSeuratObject(counts = out_Dsec2, project = "Dsec_rep2")

sc_Dsec3 <- load10X('Dsec_to_Dsim/S3_full_intron/outs/')
sc_Dsec3 = autoEstCont(sc_Dsec3)
out_Dsec3 = adjustCounts(sc_Dsec3)
Dsec_rep3_decon <- CreateSeuratObject(counts = out_Dsec3, project = "Dsec_rep3")

sc_Dsec4 <- load10X('Dsec_to_Dsim/S4_full_intron/outs/')
sc_Dsec4 = autoEstCont(sc_Dsec4)
out_Dsec4 = adjustCounts(sc_Dsec4)
Dsec_rep4_decon <- CreateSeuratObject(counts = out_Dsec4, project = "Dsec_rep4")

sc_Dsec5 <- load10X('Dsec_to_Dsim/S5_full_intron/outs/')
sc_Dsec5 = autoEstCont(sc_Dsec5)
out_Dsec5 = adjustCounts(sc_Dsec5)
Dsec_rep5_decon <- CreateSeuratObject(counts = out_Dsec5, project = "Dsec_rep5")

sc_Dsec6 <- load10X('Dsec_to_Dsim/S6_full_intron/outs/')
sc_Dsec6 = autoEstCont(sc_Dsec6)
out_Dsec6 = adjustCounts(sc_Dsec6)
Dsec_rep6_decon <- CreateSeuratObject(counts = out_Dsec6, project = "Dsec_rep6")

# Merge dataset
Dsec_decon_combined12<- merge(Dsec_rep1_decon, y = Dsec_rep2_decon, add.cell.ids = c("rep1", "rep2"), project = "Dsec_full_intron_decon")
Dsec_decon_combined34<- merge(Dsec_rep3_decon, y = Dsec_rep4_decon, add.cell.ids = c("rep3", "rep4"), project = "Dsec_full_intron_decon")
Dsec_decon_combined1234<- merge(Dsec_decon_combined12, y = Dsec_decon_combined34, project = "Dsec_full_intron_decon")
Dsec_decon_combined56<- merge(Dsec_rep5_decon, y = Dsec_rep6_decon, add.cell.ids = c("rep5", "rep6"), project = "Dsec_full_intron_decon")
Dsec_decon_combined<- merge(Dsec_decon_combined1234, y = Dsec_decon_combined56, project = "Dsec_full_intron_decon")

saveRDS(Dsec_decon_combined, file = "Processed_Data/Dsec_combined_full_intron_decon_unfiltered.rds")
Dsec_decon_combined <- readRDS(file = "Processed_Data/Dsec_combined_full_intron_decon_unfiltered.rds")

### QC check ###

Dsec_decon_combined[["percent.mt"]] <- PercentageFeatureSet(Dsec_decon_combined, pattern = "^mt:")

VlnPlot(Dsec_decon_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) ## mostly clean 
FeatureScatter(Dsec_decon_combined, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.2)
FeatureScatter(Dsec_decon_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1)

### SCT / clustering ###

## SCTransform

Dsec_decon_combined <- SCTransform(Dsec_decon_combined, method = "glmGamPoi", verbose = FALSE)

### Perform linear dimensional reduction

Dsec_decon_combined <- RunPCA(Dsec_decon_combined, features = VariableFeatures(object = Dsec_decon_combined))

### Clustering

Dsec_decon_combined <- RunUMAP(Dsec_decon_combined, dims = 1:50)
Dsec_decon_combined <- RunTSNE(Dsec_decon_combined, dims = 1:50)
Dsec_decon_combined <- FindNeighbors(Dsec_decon_combined, dims = 1:50)
Dsec_decon_combined <- FindClusters(Dsec_decon_combined, resolution = 0.8)

DimPlot(Dsec_decon_combined, reduction = "umap",  label = TRUE, pt.size = 0.01)

saveRDS(Dsec_decon_combined, file = "Processed_Data/Dsec_combined_full_intron_decon.rds")
Dsec_decon_combined <- readRDS(file = "Processed_Data/Dsec_combined_full_intron_decon.rds")

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

sweep.res.list_Dsec_decon_combined <- paramSweep_v3_SCT(Dsec_decon_combined, PCs = 1:50, sct = T)
sweep.stats_Dsec_decon_combined <- summarizeSweep(sweep.res.list_Dsec_decon_combined, GT = FALSE)
bcmvn_Dsec_decon_combined <- find.pK(sweep.stats_Dsec_decon_combined)  ## pk=0.22

## Homotypic Doublet Proportion Estimate : 15%
annotations <- Dsec_decon_combined@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.15*nrow(Dsec_decon_combined@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies

Dsec_decon_combined <- doubletFinder_v3_SCT(Dsec_decon_combined, PCs = 1:50, pN = 0.25, pK = 0.22, nExp = nExp_poi.adj, reuse.pANN = F, sct = T)

DimPlot(Dsec_decon_combined,
        pt.size = 0.01, label=F, reduction = "umap",
        group.by = "DF.classifications_0.25_0.22_7508" )+theme(aspect.ratio = 1)

saveRDS(Dsec_decon_combined, file = "Processed_Data/Dsec_combined_full_intron_decon.rds")
Dsec_decon_combined <- readRDS(file = "Processed_Data/Dsec_combined_full_intron_decon.rds")

## reclustering after doublet removal ##

Dsec_decon_combined_DF <- subset(Dsec_decon_combined, subset = DF.classifications_0.25_0.22_7508 == "Singlet")

Dsec_decon_combined_DF <- SCTransform(Dsec_decon_combined_DF, method = "glmGamPoi", verbose = FALSE)
Dsec_decon_combined_DF <- RunPCA(Dsec_decon_combined_DF, npcs = 50, verbose = FALSE)

Dsec_decon_combined_DF <- FindNeighbors(Dsec_decon_combined_DF, dims = 1:50)
Dsec_decon_combined_DF <- FindClusters(Dsec_decon_combined_DF, resolution = 0.8)

Dsec_decon_combined_DF <- RunUMAP(Dsec_decon_combined_DF, dims = 1:50)

DimPlot(Dsec_decon_combined_DF, reduction = "umap",  label = TRUE, pt.size = 0.01)

VlnPlot(Dsec_decon_combined_DF, features=c('nCount_RNA', 'nFeature_RNA'), pt.size=0, ncol=1) + NoLegend()

saveRDS(Dsec_decon_combined_DF, file = "Processed_Data/Dsec_combined_full_intron_decon_DF.rds")


###########  

## Subset KC and check doublet ##

FeaturePlot(Dsec_decon_combined, feature='Pka-R2', order=T)

Dsec_celltype_KC_doublet <- subset(Dsec_decon_combined, idents = c(5,8,16))

Dsec_celltype_KC_doublet <- SCTransform(Dsec_celltype_KC_doublet, method = "glmGamPoi", verbose = FALSE)

Dsec_celltype_KC_doublet <- RunPCA(Dsec_celltype_KC_doublet, npcs = 50, verbose = FALSE)
# t-SNE and Clustering
Dsec_celltype_KC_doublet <- RunUMAP(Dsec_celltype_KC_doublet, reduction = "pca", dims = 1:50)
Dsec_celltype_KC_doublet <- FindNeighbors(Dsec_celltype_KC_doublet, reduction = "pca", dims = 1:50)
Dsec_celltype_KC_doublet <- FindClusters(Dsec_celltype_KC_doublet, resolution = 0.1)

DimPlot(object=Dsec_celltype_KC_doublet, reduction="umap")

DimPlot(Dsec_celltype_KC_doublet,
        pt.size = 0.3, label=F, reduction = "umap",
        group.by = "DF.classifications_0.25_0.22_7508" )+theme(aspect.ratio = 1)

FeaturePlot(object=Dsec_celltype_KC_doublet, feature="rn")

### re-clustering after doublet filtering ###

Dsec_celltype_KC_DF <- subset(Dsec_celltype_KC_doublet, subset = DF.classifications_0.25_0.22_7508 == "Singlet")

Dsec_celltype_KC_DF <- SCTransform(Dsec_celltype_KC_DF, method = "glmGamPoi", verbose = FALSE)
Dsec_celltype_KC_DF <- RunPCA(Dsec_celltype_KC_DF, npcs = 50, verbose = FALSE)

Dsec_celltype_KC_DF <- FindNeighbors(Dsec_celltype_KC_DF, dims = 1:50)
Dsec_celltype_KC_DF <- FindClusters(Dsec_celltype_KC_DF, resolution = 0.2)

Dsec_celltype_KC_DF <- RunUMAP(Dsec_celltype_KC_DF, dims = 1:50)
DimPlot(Dsec_celltype_KC_DF, reduction = "umap",  label = TRUE, pt.size = 0.3)

## frequencies ##

Dsec_celltype_KC_DF_metadata <- Dsec_celltype_KC_DF@meta.data %>%
  dplyr::mutate(species=ifelse(orig.ident %in% c("Dmel_rep1", "Dmel_rep2", "Dmel_rep3", "Dmel_rep4", "Dmel_rep5", "Dmel_rep6"), "Dmel",
                               ifelse(orig.ident %in% c("Dsim_rep1", "Dsim_rep2", "Dsim_rep3", "Dsim_rep4", "Dsim_rep5", "Dsim_rep6"), "Dsim",
                                      ifelse(orig.ident %in% c("Dsec_rep1", "Dsec_rep2", "Dsec_rep3", "Dsec_rep4", "Dsec_rep5", "Dsec_rep6"), "Dsec",
                                             "DsecNoni")))) %>%
  dplyr::mutate(species=factor(species, levels = c("Dmel", "Dsim", "Dsec","DsecNoni")))

Dsec_celltype_KC_DF_metadata_summary_rep <- Dsec_celltype_KC_DF_metadata %>%
  dplyr::group_by(species, seurat_clusters, orig.ident) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(species, orig.ident) %>%
  dplyr::mutate(total=sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(percent = n/total*100)

Dsec_celltype_KC_DF_metadata_summary <- Dsec_celltype_KC_DF_metadata_summary_rep %>%
  dplyr::group_by(species, seurat_clusters) %>%
  dplyr::summarise(percent_combined = sum(n)/sum(total)*100, sem=sd(percent)/sqrt(6)) %>%
  dplyr::ungroup()

Dsec_celltype_KC_DF_metadata_summary %>%
  ggplot(.) +
  geom_col(aes(x=seurat_clusters, y=percent_combined, fill=species), alpha = 1, position = position_dodge(width=0.9)) + 
  geom_errorbar(aes(x=seurat_clusters, y=percent_combined,
                    ymin=ifelse(percent_combined-sem<0,0,percent_combined-sem), 
                    ymax=percent_combined+sem, group=species), position = position_dodge(width=0.9), width=0.7) +
  geom_point(data=Dsec_celltype_KC_DF_metadata_summary_rep, size = 1, alpha =0.8,  
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
