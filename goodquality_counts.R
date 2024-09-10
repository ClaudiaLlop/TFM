library(cellity)
library(ggplot2)

#########################################################

edited_assess_cell_quality_PCA <- function (features, file = "", quan = 1/2, alpha = 0.025) 
{
  pca <- prcomp(features, scale = TRUE, center = TRUE)
  pca_var_explained <- summary(pca)
  pcout_c <- mvoutlier::pcout(pca$x[, 1:2])
  dimens <- min(10, ncol(features))
  low_qual_i <- which(pcout_c$wfinal01 == 0)
  uni_2 <- (edited_uni_plot(pca$x[, 1:2], symb = FALSE, quan, alpha))
  low_qual_i <- which(uni_2$outliers == TRUE)
  mtdna <- NA
  mtdna_i <- grep("mtDNA", colnames(features))
  if (length(mtdna_i) > 0) {
    mtdna <- t.test(features[, mtdna_i][low_qual_i], features[, 
                                                              mtdna_i][-low_qual_i], alternative = "greater")$p.value
  }
  mapped_prop_i <- grep("Mapped", colnames(features))
  mapped_prop <- NA
  if (length(grep("Mapped", colnames(features))) > 0) {
    mapped_prop <- t.test(features[, mapped_prop_i][low_qual_i], 
                          features[, mapped_prop_i][-low_qual_i], alternative = "less")$p.value
  }
  types <- rep(1, nrow(features))
  if (!is.na(mapped_prop) && !is.na(mtdna)) {
    if (mapped_prop > 0.5 && mtdna > 0.5) {
      types[-low_qual_i] <- 0
    }else {
      types[low_qual_i] <- 0
    }
  }else {
    popul_1 <- length(low_qual_i)
    popul_2 <- nrow(features) - length(low_qual_i)
    if (popul_1 < popul_2) {
      types[low_qual_i] <- 0
    }else {
      types[-low_qual_i] <- 0
    }
  }
  annot <- data.frame(cell = rownames(features), quality = types)
  if (file != "") {
    col <- c(`0` = "red", `1` = "darkgreen")
    cellity:::plot_pca(features, as.character(types), pca, col, output_file = file)
  }
  return(annot)
}


#########################################################

edited_uni_plot <- function (x, symb = FALSE, quan = 1/2, alpha = 0.025) 
{
  if (!is.matrix(x) && !is.data.frame(x)) 
    stop("x must be matrix or data.frame")
  if (ncol(x) < 2) 
    stop("x must be at least two-dimensional")
  if (ncol(x) > 10) 
    stop("x should not be more than 10-dimensional")
  rob <- robustbase:::covMcd(x, alpha = quan)
  xarw <- mvoutlier::arw(x, rob$center, rob$cov, alpha = alpha)
  dist <- mahalanobis(x, center = rob$center, cov = rob$cov)
  sx <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) sx[, i] <- (x[, i] - xarw$m[i])/sqrt(xarw$c[i, 
                                                                   i])
  r <- range(sx)
  if (symb == FALSE) {
    for (i in 1:ncol(x)) {
      o <- (sqrt(dist) > min(sqrt(xarw$cn), sqrt(qchisq(0.975, 
                                                        dim(x)[2]))))
      l <- list(outliers = o, md = sqrt(dist))
    }
  }
  if (symb == TRUE) {
    rd <- sqrt(dist)
    lpch <- c(3, 3, 16, 1, 1)
    lcex <- c(1.5, 1, 0.5, 1, 1.5)
    xs <- scale(x) - min(scale(x))
    eucl <- sqrt(apply(xs^2, 1, sum))
    rbcol <- rev(rainbow(nrow(x), start = 0, end = 0.7))[as.integer(cut(eucl, 
                                                                        nrow(x), labels = 1:nrow(x)))]
    o <- (sqrt(dist) > min(sqrt(xarw$cn), sqrt(qchisq(0.975, 
                                                      dim(x)[2]))))
    l <- list(outliers = o, md = sqrt(dist), euclidean = eucl)
  }
  par(yaxt = "s")
  l
}

#########################################################

# As input we need the features matrix from the processed BAM files
features=read.csv("/Users/.../TFM/data/features/total_features.csv", row.names = 1)
colnames(features) = c("Mapped %", "Intergenic %", "Intronic %", "Exonic %", "Ambigious %", "#Detected genes",
                       "Cell-to-mean", "Transcriptome variance",
                       "High expr.var genes intv.1", "High expr.var genes intv.2", "High expr.var genes intv.3",
                       "High expr.var genes intv.4", "High expr.var genes intv.5", "#High exp + var genes",
                       "Intragenic %", "Apoptotic Process %", "Metabolic Process %",
                       "Ribosome %", "Membrane %", "Cytoplasm %", "Extracellular Region %", "mtDNA %",                    
                       "mitochondria_upreg %", "mitochondria_downreg %", "Actb %", "Gadph %",
                       "Volume-surface ratio" )

# We remove the cells labeled as RBC by Goodwin.
not5=read.csv("/Users/claudia/Documents/TFM/data/cells_not5.csv", row.names = 1)
colnames(not5) = "cell"
features_not5 = features[which(rownames(features) %in% not5$cell),]
features_scaled <- scale(features_not5)
pca_result <- prcomp(features_scaled, scale. = TRUE, center = TRUE)
pca_data <- data.frame(pca_result$x)
residuals_list <- list()

for (i in 1:ncol(features_scaled)) {
  feature_name <- colnames(features_scaled)[i]
  feature_data <- features_scaled[, i]
  regression_model <- lm(feature_data ~ PC1 + PC2 +PC3, data = pca_data)
  residuals_list[[feature_name]] <- residuals(regression_model)
}

residuals_df <- as.data.frame(residuals_list)
colnames(residuals_df) = colnames(features)

#write.csv(residuals_df, "/Users/.../TFM/data/residuals/regpc3_not5.csv")
annot = edited_assess_cell_quality_PCA(residuals_df, quan=1, alpha=0.025)
goodcells = annot$cell[which(annot$quality == 1)]
#write.csv(goodcells, "/Users/.../TFM/data/goodqualitycells.csv")

#load("~/.../TFM/countsnm_stats.RData")
counts12 <- counts12[,colnames(counts12) %in% goodcells]
counts34 <- counts34[,colnames(counts34) %in% goodcells]
counts56 <- counts56[,colnames(counts56) %in% goodcells]
counts78 <- counts78[,colnames(counts78) %in% goodcells]

goodcounts <- cbind(counts12, counts34)
goodcounts <- cbind(counts, counts56)
goodcounts <- cbind(counts, counts78)


#write.csv(goodcounts, "/Users/.../TFM/data/counts/goodcounts_update.csv")
