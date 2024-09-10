library(cellity)

load("~/Documents/TFM/countsnm_stats.RData")

#########################################################
## Some edited Cellity functions

edited_simple_cap <- function (x) 
{
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", 
        collapse = " ")
}

edited_sum_prop <- function (counts, genes_interest) 
{
  genes_interest_i <- which(rownames(counts) %in% unlist(genes_interest))
  genes_interest_counts_prop <- colSums(counts[genes_interest_i, 
  ])
  return(genes_interest_counts_prop)
}

edited_feature_generation <- function (counts_nm, read_metrics, GO_terms, extra_genes, organism) 
{
  features <- list()
  read_metrics <- data.frame(read_metrics)
  counts_nm <- data.frame(counts_nm)
  genes_mean <- rowMeans(counts_nm)
  genes_zero <- which(genes_mean == 0)
  if (length(genes_zero) > 0) {
    genes_mean <- genes_mean[-genes_zero]
    counts_nm_mean <- counts_nm[-genes_zero, ]/genes_mean
  }else {
    counts_nm_mean <- counts_nm/genes_mean
  }
  ercc_counts <- read_metrics$ercc
  if (is.null(ercc_counts)) {
    ercc_counts <- 0
  }
  number_mapped_reads_prop <- ((read_metrics$mapped - ercc_counts)/read_metrics$total)
  detected_genes <- apply(counts_nm, 2, function(x) {
    return(length(which(x > 0)))
  })
  counts_nm_mean_log <- log(counts_nm_mean + 0.001)
  genes_var <- apply(counts_nm_mean_log, 1, var)
  genes_means_log <- log(genes_mean)
  cell_to_mean_corr_spearman <- cor(counts_nm, rowMeans(counts_nm), 
                                    method = "spearman")
  i <- which(genes_var > quantile(genes_var)[4] & genes_means_log > 
               quantile(genes_means_log)[4])
  counts_nm_mean_log_high_var_mean <- counts_nm_mean_log[i, 
  ]
  transcriptome_variance <- matrix(0, ncol(counts_nm_mean_log_high_var_mean))
  number_of_highly_expressed_variable_genes <- matrix(0, ncol(counts_nm_mean_log_high_var_mean))
  num_of_high_var_exp_genes_interval <- matrix(0, ncol(counts_nm_mean_log_high_var_mean))
  if (length(i) > 100) {
    transcriptome_variance <- apply(counts_nm_mean_log_high_var_mean, 
                                    2, var)
    m <- colMeans(counts_nm_mean[i, ])
    number_of_highly_expressed_variable_genes <- colSums(counts_nm_mean[i, ] > m)
    
    num_of_high_var_exp_genes_interval <- apply(counts_nm_mean_log_high_var_mean, 
                                                2, function(x) {
                                                  hst <- hist(x, breaks = c(-100, -4, -2, 0, 2, 
                                                                            100), plot = FALSE)
                                                  hst$counts
                                                })
    num_of_high_var_exp_genes_interval <- t(num_of_high_var_exp_genes_interval)
    colnames(num_of_high_var_exp_genes_interval) <- paste0("num_of_high_var_exp_genes_interval_", 
                                                           1:ncol(num_of_high_var_exp_genes_interval))
  }
  mean_ex <- apply(counts_nm, 1, mean)
  i <- order(mean_ex, decreasing = FALSE)
  mean_ex <- mean_ex[i]
  lowl_expr <- mean_ex[1:(length(mean_ex) * 0.01)]
  l_i <- which(rownames(counts_nm) %in% names(lowl_expr))
  cell_to_mean_corr_spearman_low_ex = matrix(0, ncol(counts_nm))
  if (length(l_i) > 100) {
    cell_to_mean_corr_spearman_low_ex <- cor(counts_nm[l_i,], rowMeans(counts_nm[l_i, ]), method = "spearman")
  }
  
  read_metrics$intronic <- read_metrics$intragenic
  read_metrics$intragenic <- read_metrics$exonic + read_metrics$intronic
  
  techincal_features <- cbind(number_mapped_reads_prop,
                              read_metrics[, c("intergenic", "intragenic", "exonic", "intronic", "ambigious")], 
                              detected_genes, cell_to_mean_corr_spearman, cell_to_mean_corr_spearman_low_ex, 
                              transcriptome_variance, num_of_high_var_exp_genes_interval, 
                              number_of_highly_expressed_variable_genes)
  tech_names <- c("Mapped %", "Intergenic %", 
                  "Intragenic %", "Exonic %", "Intronic %", "Ambigious %","#Detected genes", "Cell-to-mean", "Cell-to-mean lowE", 
                  "Transcriptome variance", paste0("High expr.var genes intv.", 
                                                   1:ncol(num_of_high_var_exp_genes_interval)), "#High exp + var genes")
  colnames(techincal_features) <- tech_names
  GO_BP <- topGO::annFUN.org("BP", mapping = organism, ID = "ensembl")
  GO_CC <- topGO::annFUN.org("CC", mapping = organism, ID = "ensembl")
  GO <- c(GO_BP, GO_CC)
  go_prop <- sapply(unlist(GO_terms), function(go_id) {
    prop <- edited_sum_prop(counts_nm, unlist(GO[go_id]))
    return(prop)
  }, simplify = FALSE)
  go_prop <- do.call(cbind, go_prop)
  go_names <- Term(unlist(GO_terms))
  go_names <- sapply(go_names, edited_simple_cap)
  colnames(go_prop) <- go_names
  m_i <- which(GO_terms[, 1] == "GO:0016020")
  c_i <- which(GO_terms[, 1] == "GO:0005737")
  volume_surface_ratio <- matrix(0, ncol(counts_nm))
  if (length(m_i) > 0 && length(c_i) > 0 && sum(go_prop[, c_i]) > 
      0) {
    volume_surface_ratio <- go_prop[, m_i]/go_prop[, c_i]
  }
  extra_genes_prop <- sapply(extra_genes, function(extra_g) {
    prop <- edited_sum_prop(counts_nm, extra_g)
    return(prop)
  }, simplify = FALSE)
  extra_genes_prop <- do.call(cbind, extra_genes_prop)
  colnames(extra_genes_prop) <- unlist(names(extra_genes))
  biological_features <- cbind(go_prop, extra_genes_prop)
  colnames(biological_features) <- paste0(colnames(biological_features), 
                                          " %")
  features <- data.frame(techincal_features, biological_features, 
                         volume_surface_ratio)
  colnames(features) <- c(colnames(techincal_features), colnames(biological_features), 
                          "Volume-surface ratio")
  rownames(features) <- colnames(counts_nm)
  return(features)
}

edited_extract_features <- function (counts_nm, read_metrics, prefix = "", output_dir = "", common_features = NULL, GO_terms = NULL, extra_genes = NULL, 
                                     organism = "mouse") 
{
  feature_info <- get("feature_info")
  if (is.null(common_features)) {
    common_features <- feature_info[[2]]
  }
  if (is.null(GO_terms)) {
    GO_terms <- feature_info[[1]]
  }
  if (organism == "human" || organism == "org.Hs.eg.db") {
    organism <- "org.Hs.eg.db"
  }else{
    if (organism == "mouse" || organism == "org.Mm.eg.db") {
      organism <- "org.Mm.eg.db"
    }else {
      print.warnings("You have specified a different organism than mouse or human.\n\n           This might work, but you have to make sure you have specified the appropiate database as organism (e.g. org.Hs.eg.db), and you have it also installed.\n \n           Also, pleae note that extra_genes need to match the organism of interest.")
    }
  }
  if (is.null(extra_genes)) {
    if (organism == "org.Hs.eg.db") {
      extra_genes <- get("extra_human_genes")
    }
    else if (organism == "org.Mm.eg.db") {
      extra_genes <- get("extra_mouse_genes")
    }
  }
  print("Extracting features")
  genes <- rownames(counts_nm)
  if (is.null(genes) | length(genes) == 0) {
    print("Please annotate your expression matrix with genes identifiers as rownames")
    return(NULL)
  }
  features_all <- edited_feature_generation(counts_nm, read_metrics, 
                                            GO_terms, extra_genes, organism)
  print(paste0("Features extracted."))
  sds <- apply(features_all, 2, sd)
  features_all <- features_all[, !is.na(sds)]
  sds <- sds[!is.na(sds)]
  features_all <- features_all[, sds != 0]
  types <- c("all", "common")
  features_common <- features_all[, which(colnames(features_all) %in% 
                                            common_features)]
  if (prefix != "" && output_dir != "") {
    o <- paste(output_dir, prefix, sep = "/")
    print(output_dir)
    dir.create(o, showWarnings = TRUE, recursive = TRUE)
    f_all <- file.path(o, paste0(prefix, ".", types[1], ".features"))
    f_common <- file.path(o, paste0(prefix, ".", types[2], 
                                    ".features"))
    write.table(features_common, f_common)
    write.table(features_all, f_all)
    .info(paste0("Features saved: ", f_all))
    .info(paste0("Features saved: ", f_common))
  }
  return(list(features_all, features_common))
}

#########################################################
# Apply edited Cellity functions to our BAM processed data and obtain necessary 
# biological and technical features to asses cell quality

sample_features12 <- edited_extract_features(counts_nm12, stats12)
sample_features34 <- edited_extract_features(counts_nm34, stats34)
sample_features56 <- edited_extract_features(counts_nm56, stats56)
sample_features78 <- edited_extract_features(counts_nm78, stats78)

features = rbind(sample_features12[[1]], sample_features34[[1]])
features = rbind(features, sample_features56[[1]])
features = rbind(features, sample_features78[[1]])

#########################################################

write.csv(features, "/Users/.../TFM/data/features/total_features_update.csv")

