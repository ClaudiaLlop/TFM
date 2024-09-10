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

# Create Seurat object
seu.object <- CreateSeuratObject(data_feat, min.cells = 0, min.features = 0)

seu.object@meta.data$file=str_sub(rownames(seu.object@meta.data),start=-1)
seu.object@meta.data$sample=1
seu.object@meta.data$sample[which(seu.object@meta.data$file %in% c(3,4))]=2
seu.object@meta.data$sample[which(seu.object@meta.data$file %in% c(5,6))]=3
seu.object@meta.data$sample[which(seu.object@meta.data$file %in% c(7,8))]=4

seu.object[["RNA"]] <- split(seu.object[["RNA"]], f = seu.object$sample)

seu.object[["percent.mt"]] <- PercentageFeatureSet(seu.object, pattern = "^mt-")

nHVG<-2000
seu.object <- SCTransform(seu.object, method = "glmGamPoi", vst.flavor ="v2", vars.to.regress ="percent.mt", variable.features.n = nHVG)
seu.object <- RunPCA(seu.object, npcs = 30, verbose = F)


seu.object <- IntegrateLayers(
  object = seu.object, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE
)


seu.object<-IntegrateLayers(object = seu.object, method = CCAIntegration, normalization.method = "SCT",new.reduction = "integrated.cca")

dim=20
res=0.4
seu.object<- FindNeighbors(seu.object, dims= 1:20, reduction = "integrated.dr")
seu.object<- FindClusters(seu.object, resolution = res)
seu.object<- RunUMAP(seu.object, dims = 1:dim, reduction = "integrated.dr")

DimPlot(seu.object, reduction = 'umap', label=TRUE,group.by="sample")
DimPlot(seu.object, reduction = 'umap', label=TRUE,group.by="seurat_clusters")+ ggtitle("PF Zone Clusters") + NoLegend()

save(file="integrated_object.Rdata",seu.object)

#################################################
## Plot markers expression and enumerate clusters by expression levels

mar11 = FindMarkers(seu.object, ident.1 = 11)
FeaturePlot(seu.object, features = c("Col1a1"), cols = c("#fcec03", "#fc030b"), label=TRUE)
cluster_avg_expression <- AggregateExpression(seu.object, features = mm, group.by = "seurat_clusters")
rev(order(cluster_avg_expression$RNA["Ngfr", ])) - 1

#################################################
# Plot comparing expression from two markers

cells13 = WhichCells(seu.object, idents = 13)
exp13 = FetchData(seu.object, vars= c("Pdgfra", "Kdr"), cells = cells13, layer = "scale.data")
ggplot(exp13, aes(x = Pdgfra, y = Kdr)) + geom_point(color = "blue", alpha=0.6) +
  labs(title="Expression of Pdgfra vs Kdr in Cluster 13", 
       x = "Pdgfra expression", y = "Kdr expression") +
  theme_minimal()

#################################################
# List of markers
cm_markers = c("Tnni3","Tnnt2","Actn2")
sa_nodal_markers = c("Hcn4","Hcn1","Gjc1","Shox2","Isl1","Tbx3","Tbx18","Bmp2","Cacna2d2")
endo_markers = c("Npr3", "Plvap", "Cdh5", "Pecam1", "Kdr", "Fabp4")
fibro_markers = c("Tcf21","Pdgfra","Col1a1","Col1a2","Postn","Dcn")
coronarySMC_markers = c("Pdgfrb", "Acta2", "Tagln", "Des", "Myh11")
epicardial_markers = c("Upk3b","Wt1","Tbx18","Tcf21","Wnt5a","Prrx1")
wbc_markers = c("Csf1r","Ncam1","Cd34","Cd38","Cd83")
rbc_markers = c("Hba-a1","Hba-a2","Phox2b","Sox10")

#################################################
# Table showing cluster with max expression of cell type markers 
avg_exp_smc = AverageExpression(seu.object, features=epicardial_markers,slot="counts",  return.seurat = FALSE)
avg_exp_smc_df = as.data.frame(avg_exp_smc$RNA)
max_exp_smc = apply(avg_exp_smc_df, 1, function(x){
  cluster = names(which.max(x))
  value = max(x)
  return(c(cluster=cluster, max_exp=value))
})
max_exp_smc_df = as.data.frame(t(max_exp_smc))

#################################################
# Plot violin, tsne and umap 
for (marker in epicardial_markers){
  path_mark = paste0("/Users/.../TFM/data/markers/zona3/", marker)
  violin_plot <- VlnPlot(seu.object, features = marker, group.by = "seurat_clusters", pt.size=0.1) +
    theme(legend.position = "none") +
    labs(title = marker, x= "Cluster", y="Expression level")
  tsne_plot <- FeaturePlot(seu.object, features = marker, reduction = "tsne") + 
    labs(title = marker, x= "tSNE_1", y="tSNE_2")
  umap_plot <- FeaturePlot(seu.object, features = c(marker), cols = c("#fcec03", "#fc030b"), label=TRUE)
  ggsave(filename = paste0(path_mark, "_violin.png"), plot= violin_plot, width = 8, height = 6)
  ggsave(filename = paste0(path_mark, "_tsne.png"), plot= tsne_plot, width = 6, height = 4)
  ggsave(filename = paste0(path_mark, "_umap.png"), plot= umap_plot, width = 6, height = 4)
}

#################################################
# Subclustering of cluster 12 zone 1

seu.object = seu.object1
cluster12 = subset(seu.object, idents = 12)
cluster12 = RunPCA(cluster12, npcs = 50, verbose = F)
cluster12 = FindNeighbors(cluster12, dims=1:10)
cluster12 = FindClusters(cluster12, resolution=0.5)
cluster12 = RunUMAP(cluster12, dims=1:10)
DimPlot(cluster12, reduction = "umap", label = TRUE,group.by="seurat_clusters")

seu.object$subclusters = Idents(cluster12)
DimPlot(seu.object, group.by = "subclusters", label=TRUE)

#################################################
# Attempt to match our clustering results to Goodwin's by contingency table
# At the end we conclude that this is innecesary since Goodwin did not labeled
# their clusters by cell type so we don't need this check

clusters_mine = data.frame(seu.object@meta.data$SCT_snn_res.0.4)
rownames(clusters_mine) = rownames(seu.object@meta.data)
colnames(clusters_mine) = "cluster"
clusters_mine$cell = sub("_.*", "",  rownames(clusters_mine))
clusters_mine$sam = sub(".*_", "", rownames(clusters_mine))

clusters_paper = t(clusters)
clusters_paper = data.frame(clusters_paper[3:nrow(clusters_paper),])
colnames(clusters_paper) = "cluster"
clusters_paper$cell = sub(".*-", "", rownames(clusters_paper))
clusters_paper$sam = sub("-.*", "",  rownames(clusters_paper))

nrow(clusters_mine)
# 36965
nrow(clusters_paper)
# 25957
length(clusters_mine$cell %in% clusters_paper$cell)
# 36965
length(clusters_paper$cell %in% clusters_mine$cell)
# 25957
length(unique(clusters_mine$cell))
# 18811
length(unique(clusters_paper$cell))
# 25606

pair12 = colSums(contingency_table_2[c(1,2),c(2,3,4,5)])
pair34 = colSums(contingency_table_2[c(3,4),c(2,3,4,5)])
pair56 = colSums(contingency_table_2[c(5,6),c(2,3,4,5)])
pair78 = colSums(contingency_table_2[c(7,8),c(2,3,4,5)])
ct2 <- data.frame("pair12" = pair12, "pair34" = pair34, "pair56" = pair56, "pair78" = pair78)

combined_clusters <- clusters_mine %>%
  inner_join(clusters_paper, by = "cell", suffix = c("_mine", "_paper"), relationship = "many-to-many")

contingency_table <- combined_clusters %>%
  group_by(cluster_mine, cluster_paper) %>%
  summarize(count = n(), .groups = 'drop') %>%
  tidyr::spread(key = cluster_paper, value = count, fill = 0)

heatmap_data <- contingency_table %>%
  tidyr::gather(key = "cluster_paper", value = "count", -cluster_mine)
ggplot(heatmap_data, aes(x = cluster_mine, y = cluster_paper, fill = count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Clusters Obtained", y = "Clusters from Paper", fill = "Cell Count") +
  theme_minimal()


vec_rowsums <- rowSums(contingency_table)

normalized_table <- contingency_table %>%
  mutate(across(-cluster_mine, ~ . / sum(.), .names = "normalized_{.col}"))
heatmap_data_nm <- normalized_table %>%
  tidyr::pivot_longer(-cluster_mine, names_to = "cluster_paper", values_to = "normalized_count")
heatmap_data_nm <- heatmap_data_nm %>%
  mutate(cluster_paper = as.integer(gsub("normalized_", "", cluster_paper)))


ggplot(heatmap_data_nm, aes(x = factor(cluster_mine), y = factor(cluster_paper), fill = normalized_count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", trans = "log10", na.value = "white") +
  labs(x = "Cluster Obtained", y = "Clusters from Paper", fill = "Normalized Cell Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#################################################
# Check how many cells coincide in our dataset and in Goodwin's after different QC filtering

clusters_mine <- clusters_mine %>%
  mutate(in_paper = ifelse(cell %in% clusters_paper$cell, "In Paper", "Not in Paper"))
clusters_paper <- clusters_paper %>%
  mutate(in_mine = ifelse(cell %in% clusters_mine$cell, "In Mine", "Not in Mine"))

counts_mine <- clusters_mine %>%
  group_by(cluster, in_paper) %>%
  summarise(count = n(), .groups = 'drop')

ggplot(counts_mine, aes(x = factor(cluster), y = count, fill = in_paper)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), 
            size = 3, angle = 90) +
  labs(title = "Cells from paper in my clusters", 
       x = "Cluster", 
       y = "# of cells", 
       fill = "Presence") +
  theme_minimal() + theme()

counts_paper <- clusters_paper %>%
  group_by(cluster, in_mine) %>%
  summarise(count = n(), .groups = 'drop')

ggplot(counts_paper, aes(x = factor(cluster), y = count, fill = in_mine)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), 
            size = 3, angle = 90) +
  labs(title = "Cells from my dataset in the paper's clusters", 
       x = "Cluster", 
       y = "# of cells", 
       fill = "Presence") +
  theme_minimal()
