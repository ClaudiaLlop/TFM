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

goodcounts = read.csv("goodcounts_update.csv", row.names=1)
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

chec=str_sub(colnames(data_feat),start=-1)
# Keeping samples from Zone SAN
data_feat = data_feat[,which(str_sub(colnames(data_feat),start=-1) == 1)]

# Create Seurat Object
seu.object <- CreateSeuratObject(data_feat, min.cells = 0, min.features = 0)

seu.object@meta.data$sample=1
seu.object[["percent.mt"]] <- PercentageFeatureSet(seu.object, pattern = "^mt-")

nHVG<-2000
seu.object <- SCTransform(seu.object, method = "glmGamPoi", vst.flavor ="v2", vars.to.regress ="percent.mt", variable.features.n = nHVG)
seu.object <- RunPCA(seu.object, npcs = 50, verbose = F)
ElbowPlot(seu.object, ndims = 50, reduction = "pca")

dim=25
res=0.5
seu.object<- FindNeighbors(seu.object, dims= 1:dim)
seu.object<- FindClusters(seu.object, resolution = res)
seu.object<- RunUMAP(seu.object, dims = 1:dim)
set.seed(10)

seu.object<- RunTSNE(seu.object, dims = 1:dim,perplexity=30)

DimPlot(seu.object, reduction = 'tsne', label=TRUE,group.by="seurat_clusters")

# Quick check of markers expression. Not definitive labels.

FeaturePlot(object = seu.object,  reduction = 'tsne', features = "Hcn4")
## Cluster 10 is SAN
FeaturePlot(object = seu.object,  reduction = 'tsne', features = "Actn2")
## Clusters 0, 2, 4, 5 7, 15 Cardiomyocytes (atrial, since this is sino-atrial node, in fact we check:
FeaturePlot(object = seu.object,  reduction = 'tsne', features = "Myh6")
##  Cluster 11 epicardial
FeaturePlot(object = seu.object,  reduction = 'tsne', features = "Upk3b")
##  Cluster 6 endocardial,
FeaturePlot(object = seu.object,  reduction = 'tsne', features = "Npr3")
##  Cluster 9 and maybe 8 endothelial
FeaturePlot(object = seu.object,  reduction = 'tsne', features = "Tgfb1")
## Cluster 8 WBCs
FeaturePlot(object = seu.object,  reduction = 'tsne', features = "Csf1r")
##  Clusters 1, 3 fibroblasts
FeaturePlot(object = seu.object,  reduction = 'tsne', features = "Tcf21")
FeaturePlot(object = seu.object,  reduction = 'tsne', features = "Pdgfra")
##  Cluster 15 RBCs
FeaturePlot(object = seu.object,  reduction = 'tsne', features = "Hba-a1")
##  Cluster 16 neurons.
FeaturePlot(object = seu.object,  reduction = 'tsne', features = "Ngfr")
## Less obvious labelling that need more analysis: 12, 13, 14 

save(file="zone1.Rdata",seu.object)

