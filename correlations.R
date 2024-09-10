library(ggplot2)
library(tidyr)
load("~/Documents/TFM/Rdatas/expression_df.RData")
mito_genes <-c("mt-Co2","mt-Cytb","mt-Co3","mt-Nd2","mt-Nd1","mt-Nd4","mt-Atp6" ,"mt-Co1","mt-Nd5",
               "mt-Nd3" ,"mt-Nd4l" ,"mt-Atp8", "mt-Nd6" )
keep_genes = read.csv("/Users/.../TFM/data/keep_genes.csv", row.names=1)

################################################
## correlation mtgenes vs all rest

mito_genes2 <-c("mt-Co2","mt-Cytb","mt-Co3","mt-Nd2","mt-Nd1","mt-Nd4","mt-Atp6" ,"mt-Co1","mt-Nd5",
               "mt-Nd3" )

mito_data <- expression_df[, mito_genes2]
rest_genes_data <- expression_df[, keep_genes$x]
rest_genes_data <- rest_genes_data[, !colnames(rest_genes_data) %in% mito_genes2]

r_matrix <- matrix(NA, nrow = 10, ncol = ncol(rest_genes_data))
pvalue_matrix <- matrix(NA, nrow = 10, ncol = ncol(rest_genes_data))
fdr_matrix <- matrix(NA, nrow = 10, ncol = ncol(rest_genes_data))

rownames(r_matrix) <- mito_genes2
colnames(r_matrix) <- colnames(rest_genes_data)
rownames(pvalue_matrix) <- mito_genes2
colnames(pvalue_matrix) <- colnames(rest_genes_data)
rownames(fdr_matrix) <- mito_genes2
colnames(fdr_matrix) <- colnames(rest_genes_data)

for (i in 1:10) {
  for (j in 1:ncol(rest_genes_data)) {
    test_result <- cor.test(mito_data[, i], rest_genes_data[, j], method = "pearson")
    r_matrix[i, j] <- test_result$estimate
    pvalue_matrix[i, j] <- test_result$p.value
  }
}

fdr_matrix <- apply(pvalue_matrix, 1, p.adjust, method = "BH")
fdr_matrix <- t(fdr_matrix)  

r_df <- as.data.frame(r_matrix)
pvalue_df <- as.data.frame(pvalue_matrix)
fdr_df <- as.data.frame(fdr_matrix)

path_save = "/Users/.../TFM/data/CoExp/Coexp_mt_vs_todos"
write.csv(r_df, paste0(path_save, "/pearson_correlation_r_values_filtered.csv"))
write.csv(pvalue_df, paste0(path_save, "/pearson_correlation_p_values_filtered.csv"))
write.csv(fdr_df, paste0(path_save, "/pearson_correlation_fdr_values_filtered.csv"))

################################################
## ordered pearson acumulated mtgenes vs all rest
legenddata = c("mt-Nd1" = "firebrick4",
               "mt-Nd2" = "red",
               "mt-Nd3" = "deeppink2",
               "mt-Nd4" = "palevioletred",
               #"mt-Nd4l" = "magenta",
               "mt-Nd5" = "plum1",
               #"mt-Nd6" = "lightpink",
               "mt-Cytb" = "gold",
               "mt-Co1" = "blue",
               "mt-Co2" = "deepskyblue",
               "mt-Co3" = "cadetblue1",
               "mt-Atp6" = "olivedrab1"
               #,"mt-Atp8" = "seagreen"
               )
r_df_sorted <- as.data.frame(t(apply(abs(r_df), 1, sort)))
r_df_sorted <- r_df_sorted[match(names(legenddata), rownames(r_df_sorted)), ]
mat <- as.matrix(r_df_sorted)

matplot(t(mat), type = "l", lty = 1, col = legenddata, xlab = "Number of genes", ylab = "Pearson coefficient")
abline(h = 0, lty = 2)

legend("topright", legend = names(legenddata), col = legenddata, lty = 1)


################################################
## correlation mtgenes vs mtgenes

mito_data <- expression_df[, mito_genes]

r_matrix <- matrix(NA, nrow = 13, ncol = 13)
pvalue_matrix <- matrix(NA, nrow = 13, ncol = 13)
fdr_matrix <- matrix(NA, nrow = 13, ncol = 13)

rownames(r_matrix) <- mito_genes
colnames(r_matrix) <- mito_genes
rownames(pvalue_matrix) <- mito_genes
colnames(pvalue_matrix) <- mito_genes
rownames(fdr_matrix) <- mito_genes
colnames(fdr_matrix) <- mito_genes

for (i in 1:13) {
  for (j in 1:13) {
    if (i != j) { 
      test_result <- cor.test(mito_data[, i], mito_data[, j], method = "pearson")
      r_matrix[i, j] <- test_result$estimate
      pvalue_matrix[i, j] <- test_result$p.value
    } else {
      r_matrix[i, j] <- 1 
      pvalue_matrix[i, j] <- NA  
    }
  }
}

fdr_matrix <- apply(pvalue_matrix, 1, p.adjust, method = "BH")
fdr_matrix <- t(fdr_matrix)  

r_df <- as.data.frame(r_matrix)
pvalue_df <- as.data.frame(pvalue_matrix)
fdr_df <- as.data.frame(fdr_matrix)

path_save = "/Users/.../TFM/data/CoExp/Coexp_mt_vs_mt"
write.csv(r_df, paste0(path_save, "/pearson_correlation_r_values.csv"))
write.csv(pvalue_df, paste0(path_save, "/pearson_correlation_p_values.csv"))
write.csv(fdr_df, paste0(path_save, "/pearson_correlation_fdr_values.csv"))

library(pheatmap)
library(corrplot)
library(RColorBrewer)


################################################
## heatmap mtgenes vs mtgenes

r_df_heatmap = r_df[!rownames(r_df) %in% c("mt-Nd6", "mt-Nd4l", "mt-Atp8"),]
r_df_heatmap = r_df_heatmap[,!colnames(r_df_heatmap) %in% c("mt-Nd6", "mt-Nd4l", "mt-Atp8")]

min_val <- min(r_df_heatmap, na.rm = TRUE)
max_val <- max(r_df_heatmap, na.rm = TRUE)

breaks <- unique(c(seq(min_val, 0, length.out = 50), 
                   seq(0, max_val, length.out = 51)))


clustering_result = pheatmap(
  r_df_heatmap,                       
  cluster_cols = FALSE,
  clustering_method = "average",      
  color = colorRampPalette(c("white", "lemonchiffon","yellow", "orange", "red", "maroon"))(100),  
  breaks = breaks,
  show_rownames = TRUE,               
  show_colnames = TRUE,               
  treeheight_row = 150,                 
  treeheight_col = 150                 
)

row_order <- clustering_result$tree_row$order
r_df_heatmap = r_df_heatmap[,row_order]

pheatmap(
  r_df_heatmap,                        
  cluster_cols = FALSE,
  clustering_method = "average",     
  color = colorRampPalette(c("white", "lemonchiffon","yellow", "orange", "red", "maroon"))(100), 
  breaks = breaks,
  show_rownames = TRUE,               
  show_colnames = TRUE,               
  treeheight_row = 150,                 
  treeheight_col = 150                 
)

################################################
## correlation average mtgene vs all rest

mito_data <- expression_df[, mito_genes]
mito_data_avg <- mito_data[,!colnames(mito_data) %in% c("mt-Atp8", "mt-Nd6", "mt-Nd4l")]
mito_data_avg <- rowMeans(mito_data_avg)
rest_genes_data <- expression_df[, keep_genes$x]
rest_genes_data <- rest_genes_data[, !colnames(rest_genes_data) %in% mito_genes]

r_matrix <- matrix(NA, nrow = 1, ncol = ncol(rest_genes_data))
pvalue_matrix <- matrix(NA, nrow = 1, ncol = ncol(rest_genes_data))
fdr_matrix <- matrix(NA, nrow = 1, ncol = ncol(rest_genes_data))

rownames(r_matrix) <- "AverageGene"
colnames(r_matrix) <- colnames(rest_genes_data)
rownames(pvalue_matrix) <-  "AverageGene"
colnames(pvalue_matrix) <- colnames(rest_genes_data)
rownames(fdr_matrix) <-  "AverageGene"
colnames(fdr_matrix) <- colnames(rest_genes_data)

for (j in 1:ncol(rest_genes_data)) {
  test_result <- cor.test(mito_data_avg, rest_genes_data[, j], method = "pearson")
  r_matrix[1, j] <- test_result$estimate
  pvalue_matrix[1, j] <- test_result$p.value
}

fdr_matrix <- apply(pvalue_matrix, 1, p.adjust, method = "BH")
fdr_matrix <- t(fdr_matrix)  # Volver a la forma original

r_df <- as.data.frame(r_matrix)
pvalue_df <- as.data.frame(pvalue_matrix)
fdr_df <- as.data.frame(fdr_matrix)

aux_r = r_matrix[1,]
names(aux_r) = colnames(r_matrix)

load("~/Documents/TFM/data/CoExp/sets_de_factores_ensamblaje.Rdata")
#> names(pathways_example)
#[1] "mito_genes"   "factores_I"   "factores_II"  "factores_III" "factores_IV"  "factores_V"  

plotEnrichment(pathways_example[["mito_genes"]], aux_r)

library(clusterProfiler)
library(org.Mm.eg.db)
library(fgsea)
gene_list <- sort(aux_r_vector, decreasing = TRUE)
names(gene_list) <- rownames(aux_r)  
pathways <- list(
  "Mitochondrial Genes" = rownames(aux_r)  
)
fgsea_results <- fgsea(pathways = pathways,
                       stats = gene_list)
fgsea_results <- fgsea_results[order(fgsea_results$NES, decreasing = TRUE), ]
plotEnrichment(fgsea_results, pathways[["Mitochondrial Genes"]])


