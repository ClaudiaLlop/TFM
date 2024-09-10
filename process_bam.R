library(Rsamtools)
library(stringr)
library(reshape2)
library(dplyr)
library(tidyr)
library(biomaRt)

# Function to process a BAM file and generate the statistics necessary to run Cellity
process_BAM <- function(file_path) {

  p4 <- ScanBamParam(tag=c("CB", "RE", "GN"), what="flag")
  bam <- scanBam(file_path, param=p4)
  
  df <- data.frame(cell_barcode = bam[[1]]$tag$CB,
                   region_type = bam[[1]]$tag$RE,
                   gene_names = bam[[1]]$tag$GN)
  
  df <- df[-is.na(df$cell_barcode),]
  df$region_type[is.na(df$region_type)] <- "unmapped"
  
  df$read_type <- "unmapped"
  df$read_type[which(str_detect(df$gene_names, ";"))] <- "ambigious"
  df$read_type[which((df$read_type != "ambigious") &
                       (!is.na(df$gene_names)) &
                       (df$region_type == "E"))] <- "exonic"
  df$read_type[which((df$read_type != "ambigious") &
                       (!is.na(df$gene_names)) &
                       (df$region_type == "N"))] <- "intronic"
  df$read_type[which(df$region_type == "I")] <- "intergenic"
  
  ensembl <- useMart(host="https://feb2021.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
  gene_names_unique <- unique(df$gene_names)
  ensembl_names <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
                         filters = "external_gene_name", 
                         values = gene_names_unique, 
                         mart = ensembl)
  mt_genes <- unique(ensembl_names$external_gene_name[grep("^mt-", ensembl_names$external_gene_name)])
  
  df$mt_mapped <- 0
  df$mt_mapped[which(df$gene_names %in% mt_genes)] <- 1
  
  metrics <- df %>%
    group_by(cell_barcode, read_type) %>%
    summarise(count = n()) %>%
    pivot_wider(names_from = read_type, values_from = count, values_fill = 0) %>%
    ungroup()
  
  metrics <- metrics %>%
    mutate(mapped = exonic + intronic + intergenic + ambigious)
  
  mapped_reads <- metrics$exonic + metrics$intronic + metrics$intergenic + metrics$ambigious
  metrics$exonic <- metrics$exonic/mapped_reads
  metrics$intronic <- metrics$intronic/mapped_reads
  metrics$intergenic <- metrics$intergenic/mapped_reads
  metrics$ambigious <- metrics$ambigious/mapped_reads
  
  mt_df <- df %>%
    group_by(cell_barcode) %>%
    summarize(mt_mapped_sum = sum(mt_mapped))
  
  metrics <- merge(metrics, mt_df, by = "cell_barcode", all.x = TRUE)
  metrics$pc_mt_mapped <- metrics$mt_mapped_sum/mapped_reads

  return(metrics)
}

# Function to process all BAM files in a folder
process_BAM_folder <- function(folder_path) {

  bam_files <- list.files(path = folder_path, pattern = "\\.bam$", full.names = TRUE)
  stats_list <- list()
  i = 1
  for (file in bam_files) {
    print(i)
    i = i+1
    stats_df <- process_BAM(file)
    stats_list[[file]] <- stats_df
  }


  combined_stats_df <- do.call(rbind, stats_list)
  
  return(combined_stats_df)
}

carpeta_bam <- "/Users/.../TFM/data/BAMcellranger"

resultados <- process_BAM_folder(carpeta_bam)

