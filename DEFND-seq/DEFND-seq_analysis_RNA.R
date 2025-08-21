##This is the file to analyse the RNA data of the DEFND-seq sequencing data. 
#The variable names will carry the numbers, the project name will give clues 
#about the sample itself
#all samples are from mice of the 36 wks cohort
#39 -> CasB6-DEN
#42 -> B6Cas-DEN
#45 -> B6Cas-neg
#48 -> CasB6-neg
#The goal of this script is to get a first look at the data, do some QC, rule 
#outliers and create a first UMAP plot to visualize the data
#for this we follow the Seurat tutorial at:
#https://satijalab.org/seurat/articles/pbmc3k_tutorial

##Libraries and Environment
save.image("/omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq/02_analysis/DEFND_eff_env.RData")
load("/omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq/02_analysis/DEFND_eff_env.RData")

library(QuasR)
library(Rsamtools)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


##import data from all samples----
rna.paths <- c('/omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq/AS-1617751-LR-80745/outs/filtered_feature_bc_matrix',
               '/omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq/AS-1617753-LR-80745/outs/filtered_feature_bc_matrix',
               '/omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq/AS-1617755-LR-80745/outs/filtered_feature_bc_matrix',
               '/omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq/AS-1617757-LR-80745/outs/filtered_feature_bc_matrix')

rna.cts <- lapply(rna.paths, Read10X)


#create seurat object
sample.names <- 
  c('#39 CasB6-DEN', '#42 B6Cas-DEN', '#45 B6Cas-neg', '#48 CasB6-neg')

rna.seur <- mapply(function(counts, samples){
  CreateSeuratObject(counts = counts, 
                     project = samples, 
                     min.cells = 3, 
                     min.features = 200)
}, rna.cts, sample.names, SIMPLIFY = FALSE)


##basic QC##----
#We look at the total amount of RNA (nCount_RNA), the number of detected 
#features and number of mt-RNAs per single cell
#Generally, nFeatures between 200 and 2500 per cell are regarded as reliable
#cells with over 5% of mitochondria RNA are regarded as unreliable

rna.seur <- lapply(rna.seur, function(x){
  x[['percent.mt']] <- PercentageFeatureSet(x, pattern = '^mt-')
  return(x)
})

#1: CasB6-DEN, #2: B6Cas-DEN, #3: B6Cas-neg, #4: CasB6-neg
VlnPlot(rna.seur[[4]], 
        features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), 
        col = 3)

FeatureScatter(rna.seur[[4]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(rna.seur[[4]], feature1 = "nCount_RNA", feature2 = "percent.mt")


#check the amount of mt RNA
hist(rna.seur[[4]]$percent.mt, breaks = 100)

#subset the data, so that only reliable cells are included
rna.seur <- lapply(rna.seur, function(x){
  subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
})

##Normalization##----
#Here we the feature expression measurements for each cell are normalized by the 
#total expression, multiplied by a scale factor (10,000 by default),
#and log-transformed

rna.seur <- lapply(rna.seur, function(x){
  NormalizeData(x, normalization.method = 'LogNormalize', scale.factor = 10000)
})


##Scaling##----
#standard pre-processing step prior to dimensional reduction techniques like PCA
#includes:- Shifting the expression of each gene, 
#           so that the mean expression across cells is 0
#         - Scales the expression of each gene -> variance across cells is 1 
all.genes <- lapply(rna.seur, rownames)
rna.seur <- mapply(function(data, features){
  ScaleData(data, features = features)
}, rna.seur, all.genes, SIMPLIFY = FALSE)

#This was computationally too intense - this is the same, but simpler to ease 
#the computational burden
top.genes.1 <- head(VariableFeatures(rna.seur[[1]]), 1000)
top.genes.2 <- head(VariableFeatures(rna.seur[[2]]), 1000)
top.genes.3 <- head(VariableFeatures(rna.seur[[3]]), 1000)
top.genes.4 <- head(VariableFeatures(rna.seur[[4]]), 1000)
rna.seur[[1]] <- ScaleData(rna.seur[[1]], top.genes.1)
rna.seur[[2]] <- ScaleData(rna.seur[[2]], top.genes.2)
rna.seur[[3]] <- ScaleData(rna.seur[[3]], top.genes.3)
rna.seur[[4]] <- ScaleData(rna.seur[[4]], top.genes.4)

##PCA##----
rna.seur <- lapply(rna.seur, function(x){
  RunPCA(x, features = VariableFeatures(object = x))
})

#1: CasB6-DEN, #2: B6Cas-DEN, #3: B6Cas-neg, #4: CasB6-neg
VizDimLoadings(rna.seur[[4]], 
               dims = 1:2, 
               reduction = "pca", 
               nfeatures = 20, 
               balanced = TRUE)
DimPlot(rna.seur[[4]], 
        reduction = 'pca') 
DimHeatmap(rna.seur[[1]], 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)

##dimensionality##----
ElbowPlot(rna.seur[[4]])

#adjust PCA cluster dimensions for each sample: 
rna.dims <- list(c(1:10), c(1:10), c(1:10), c(1:10))
rna.seur <- mapply(function(x, dims){
  FindNeighbors(x, dims = dims)
}, rna.seur, rna.dims, SIMPLIFY = FALSE)

rna.seur <- lapply(rna.seur, function(x){
  FindClusters(x, resolution = 0.5)
})

##UMAP##----
rna.seur <- mapply(function(x, dims){
  RunUMAP(x, dims = dims)
}, rna.seur, rna.dims, SIMPLIFY = FALSE)

#1: CasB6-DEN, #2: B6Cas-DEN, #3: B6Cas-neg, #4: CasB6-neg
DimPlot(rna.seur[[3]], reduction = 'umap')


##Distinguishing Features of clusters
#for 42 clusters 0, 5, 6, 7 and 8 look interesting 
#compared to 45 clusters 5, 8, 9, 6 and 2
markers.42 <- FindAllMarkers(rna.seur[[2]])
markers.42 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

#let's quickly visualize the interesting ones:
#cluster 0
FeaturePlot(rna.seur[[2]], 
            features = c('Ptprb', 'Stab2', "St6galnac3", 'Meis2', 'Fbxl7',
                         'Sema6a', 'F8', 'Cyyr1', 'Adgrl4', 'Plekhg1'))                           
#cluster 1
FeaturePlot(rna.seur[[2]], 
            features = c('Aox3', 'Cmss1', "Gm29966", 'Lars2', 'Dpyd',
                         'Gm19951', 'Egfr', 'Cdk8', 'Camk1d', 'Jarid2'))

#cluster 5
FeaturePlot(rna.seur[[3]], 
            features = c('Bmp5', 'Ank3', "Pde3a", 'Ntm', 'Reln',
                         'Nrxn1', 'Hand2os1', 'Nr1h5', 'Rbms3', 'Rbms3'))

#cluster 6
FeaturePlot(rna.seur[[3]], 
            features = c('Bmp5', 'Ank3', "Pde3a", 'Ntm', 'Reln',
                         'Nrxn1', 'Hand2os1', 'Nr1h5', 'Rbms3', 'Rbms3'))


#cluster 8
FeaturePlot(rna.seur[[3]], 
            features = c('Anxa13', 'Spp1', "Nipal2", '2610307P16Rik', 'Tmem45a',
                         'Ddit4l', 'A230001M10Rik', 'Glis3', 'Bicc1', 'Dcdc2a'))

##Fututre analysis could include:
#    - tweaking UMAP features
#    - plotting the samples in the same UMAP space to see differences
#    - use this to identify differences between neg and DEN samples
#    - identify features of clusters and assign cell types
#    - identify festures of potential tumor-cell-clusters


