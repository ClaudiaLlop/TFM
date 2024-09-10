set.seed(10)
library(biomaRt)
library(Seurat)
library(glmGamPoi)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(scales)

goodcounts = read.csv("/Users/claudia/Documents/TFM/data/counts/goodcounts.csv", row.names=1)
ensembl <- useMart(host="https://feb2021.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
gene_symbols <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                      filters = 'ensembl_gene_id',
                      values = rownames(goodcounts),
                      mart = ensembl)

not_in = rownames(goodcounts)[!rownames(goodcounts) %in% gene_symbols$ensembl_gene_id]
# "ENSMUSG00000000325" "ENSMUSG00000021745" "ENSMUSG00000029386" "ENSMUSG00000035310"
# "ENSMUSG00000039530" "ENSMUSG00000051412" "ENSMUSG00000063897" "ENSMUSG00000079065"
# "ENSMUSG00000079834" "ENSMUSG00000097971" "ENSMUSG00000098178"
data_feat = goodcounts[which(!rownames(goodcounts) %in% not_in),]
data_feat$ensembl_gene_id = rownames(data_feat)
data_feat = merge(data_feat, gene_symbols, by="ensembl_gene_id")
data_feat = data_feat[!duplicated(data_feat$external_gene_name), ] #hay 1 repetido
rownames(data_feat) = data_feat$external_gene_name
data_feat <- data_feat[, !colnames(data_feat) %in% c('external_gene_name', 'ensembl_gene_id')]
rm(ensembl, gene_symbols, not_in)

# Keeping samples from Zone PF
data_feat = data_feat[,which(str_sub(colnames(data_feat),start=-1) %in% c(3,4))]

# Create Seurat Object
seu.object <- CreateSeuratObject(data_feat, min.cells = 0, min.features = 0)

seu.object@meta.data$sample=str_sub(rownames(seu.object@meta.data),start=-1)
summary(factor(seu.object@meta.data$sample))
#   3    4
#5338 4618
seu.object[["RNA"]] <- split(seu.object[["RNA"]], f = seu.object$sample)

seu.object[["percent.mt"]] <- PercentageFeatureSet(seu.object, pattern = "^mt-")

nHVG<-2000
seu.object <- SCTransform(seu.object, method = "glmGamPoi", vst.flavor ="v2", vars.to.regress ="percent.mt", variable.features.n = nHVG)
seu.object <- RunPCA(seu.object, npcs = 50, verbose = F)
ElbowPlot(seu.object, ndims = 50, reduction = "pca")

seu.object<-IntegrateLayers(object = seu.object, method = CCAIntegration, normalization.method = "SCT",new.reduction = "integrated.cca")

dim=30
res=0.75
seu.object<- FindNeighbors(seu.object, dims= 1:dim, reduction = "integrated.cca")
seu.object<- FindClusters(seu.object, resolution = res, reduction = "integrated.cca")
seu.object<- RunUMAP(seu.object, dims = 1:dim, reduction = "integrated.cca")
set.seed(10)
seu.object<- RunTSNE(seu.object, dims = 1:dim, reduction = "integrated.cca",perplexity=30)

DimPlot(seu.object, reduction = 'tsne', label=TRUE,group.by="seurat_clusters")

# Quick check of markers expression. Not definitive labels.

## Cluster 18 is PF
FeaturePlot(object = seu.object,  reduction = 'tsne', features = "Cacna2d2")

save(file="zone3.Rdata",seu.object)

