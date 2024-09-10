library(fgsea)
library(biomaRt)
load("datafit.RData")
keep_genes = read.csv("keep_genes.csv", row.names=1) # genes that have passed QC 

cell_types <- c("Cardiomyocytes", "ConductionSystem", "CoronarySMC", 
                "Endocardial", "Epicardial", "Fibroblasts")
mito_genes <-c("mt-Co2","mt-Cytb","mt-Co3","mt-Nd2","mt-Nd1","mt-Nd4","mt-Atp6" ,"mt-Co1","mt-Nd5",
               "mt-Nd3" ,"mt-Nd4l" ,"mt-Atp8", "mt-Nd6" )
zones_list <- c("ZoneSAN", "ZoneAVN", "ZonePF")

genes = colnames(expression_df)
rm(expression_df)

mouse <- useMart(host="https://feb2021.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl")
human <- useMart(host="https://feb2021.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

load("~/Documents/TFM/data/CoExp/sets_de_factores_ensamblaje.Rdata")
DE_folder="DE_fv/"

contrasts_list=list.files(paste0(DE_folder),recursive=TRUE)

results_gsea = as.data.frame(matrix(NA, ncol = 3, nrow=63))
colnames(results_gsea) = c("NES", "pvalue", "Bonferroni")
rownames(results_gsea)=contrasts_list
results_gsea_mito=results_gsea
results_gsea_I=results_gsea
results_gsea_II=results_gsea
results_gsea_III=results_gsea
results_gsea_IV=results_gsea
results_gsea_V=results_gsea

for(i in 1:length(contrasts_list))
{
    print(i)
    contrast = read.csv(paste0(DE_folder,contrasts_list[i]), row.names=1)
    contrast = as.data.frame(contrast[rownames(contrast) %in% c(keep_genes$x, factores),])
    stats_contrast = abs(contrast$t)
    names(stats_contrast)=rownames(contrast)
    fgseaRes <- fgsea(pathways = pathways_example, stats = stats_contrast, eps=0,scoreType = "pos",  nproc = 1 )
    
    results_gsea_mito[i,] <- c(fgseaRes$NES[1], fgseaRes$pval[1], fgseaRes$pval[1]*63)
    results_gsea_I[i,] <- c(fgseaRes$NES[2], fgseaRes$pval[2], fgseaRes$pval[2]*63)
    results_gsea_II[i,] <- c(fgseaRes$NES[3], fgseaRes$pval[3], fgseaRes$pval[3]*63)
    results_gsea_III[i,] <- c(fgseaRes$NES[4], fgseaRes$pval[4], fgseaRes$pval[4]*63)
    results_gsea_IV[i,] <- c(fgseaRes$NES[5], fgseaRes$pval[5], fgseaRes$pval[5]*63)
    results_gsea_V[i,] <- c(fgseaRes$NES[6], fgseaRes$pval[6], fgseaRes$pval[6]*63)

}

results_gsea_mito$BH=p.adjust(results_gsea_mito$p,method="BH")
results_gsea_I$BH=p.adjust(results_gsea_I$p,method="BH")
results_gsea_II$BH=p.adjust(results_gsea_II$p,method="BH")
results_gsea_III$BH=p.adjust(results_gsea_III$p,method="BH")
results_gsea_IV$BH=p.adjust(results_gsea_IV$p,method="BH")
results_gsea_V$BH=p.adjust(results_gsea_V$p,method="BH")

save(results_gsea_mito,file="results_gsea_mito.Rdata")
save(results_gsea_I,file="results_gsea_I.Rdata")
save(results_gsea_II,file="results_gsea_II.Rdata")
save(results_gsea_III,file="results_gsea_III.Rdata")
save(results_gsea_IV,file="results_gsea_IV.Rdata")
save(results_gsea_V,file="results_gsea_V.Rdata")

results_gsea_mito[order(results_gsea_mito$BH),]
results_gsea_I[order(results_gsea_I$BH),]
results_gsea_II[order(results_gsea_II$BH),]
results_gsea_III[order(results_gsea_III$BH),]
results_gsea_IV[order(results_gsea_IV$BH),]
results_gsea_V[order(results_gsea_V$BH),]


#############################################
### correlogram contrasts NES y BH from results_gsea
#############################################
plot_data_color = as.data.frame(matrix(NA, ncol=6, nrow=63))
plot_data_size = as.data.frame(matrix(NA, ncol=6, nrow=63))
colnames(plot_data_color) = c("MitochondrialGenes", "ComplexI_nucleus",
                              "ComplexII_nucleus", "ComplexIII_nucleus",
                              "ComplexIV_nucleus", "ComplexV_nucleus")
colnames(plot_data_size) = c("MitochondrialGenes", "ComplexI_nucleus",
                             "ComplexII_nucleus", "ComplexIII_nucleus",
                             "ComplexIV_nucleus", "ComplexV_nucleus")

load("~/Documents/TFM/data/CoExp/Coexp_mt_vs_todos/gsea/results_gsea_mito.Rdata")
gsea_plot =results_gsea_mito
gsea_plot$neg_log10_BH = -log10(gsea_plot$BH)
gsea_plot = gsea_plot[,c("NES", "neg_log10_BH")]

rownames(plot_data_color) = rownames(gsea_plot)
rownames(plot_data_size) = rownames(gsea_plot)

plot_data_color[,"MitochondrialGenes"] = gsea_plot$NES
plot_data_size[,"MitochondrialGenes"] = gsea_plot$neg_log10_BH


{load("~/Documents/TFM/data/CoExp/Coexp_mt_vs_todos/gsea/results_gsea_I.Rdata")
  gsea_plot =results_gsea_I
  gsea_plot$neg_log10_BH = -log10(gsea_plot$BH)
  gsea_plot = gsea_plot[,c("NES", "neg_log10_BH")]
  plot_data_color[,"ComplexI_nucleus"] = gsea_plot$NES
  plot_data_size[,"ComplexI_nucleus"] = gsea_plot$neg_log10_BH
  
  load("~/Documents/TFM/data/CoExp/Coexp_mt_vs_todos/gsea/results_gsea_II.Rdata")
  gsea_plot =results_gsea_II
  gsea_plot$neg_log10_BH = -log10(gsea_plot$BH)
  gsea_plot = gsea_plot[,c("NES", "neg_log10_BH")]
  plot_data_color[,"ComplexII_nucleus"] = gsea_plot$NES
  plot_data_size[,"ComplexII_nucleus"] = gsea_plot$neg_log10_BH
  
  load("~/Documents/TFM/data/CoExp/Coexp_mt_vs_todos/gsea/results_gsea_III.Rdata")
  gsea_plot =results_gsea_III
  gsea_plot$neg_log10_BH = -log10(gsea_plot$BH)
  gsea_plot = gsea_plot[,c("NES", "neg_log10_BH")]
  plot_data_color[,"ComplexIII_nucleus"] = gsea_plot$NES
  plot_data_size[,"ComplexIII_nucleus"] = gsea_plot$neg_log10_BH
  
  load("~/Documents/TFM/data/CoExp/Coexp_mt_vs_todos/gsea/results_gsea_IV.Rdata")
  gsea_plot =results_gsea_IV
  gsea_plot$neg_log10_BH = -log10(gsea_plot$BH)
  gsea_plot = gsea_plot[,c("NES", "neg_log10_BH")]
  plot_data_color[,"ComplexIV_nucleus"] = gsea_plot$NES
  plot_data_size[,"ComplexIV_nucleus"] = gsea_plot$neg_log10_BH
  
  load("~/Documents/TFM/data/CoExp/Coexp_mt_vs_todos/gsea/results_gsea_V.Rdata")
  gsea_plot =results_gsea_V
  gsea_plot$neg_log10_BH = -log10(gsea_plot$BH)
  gsea_plot = gsea_plot[,c("NES", "neg_log10_BH")]
  plot_data_color[,"ComplexV_nucleus"] = gsea_plot$NES
  plot_data_size[,"ComplexV_nucleus"] = gsea_plot$neg_log10_BH}

library(reshape2)
library(ggplot2)

color_type = as.matrix(plot_data_color[1:45,])
size_type = as.matrix(plot_data_color[1:45,])
color_zone = as.matrix(plot_data_color[46:63,])
size_zone = as.matrix(plot_data_color[46:63,])

color_type_long <- melt(color_type)
colnames(color_type_long) <- c("Contrast", "Set", "NES")
size_type_long <- melt(size_type)
colnames(size_type_long) <- c("Contrast", "Set", "neg_log10_BH")
plot_data_type <- merge(color_type_long, size_type_long, by = c("Contrast", "Set"))
colnames(plot_data_type) <- c("Contrast", "Set", "NES",  "neg_log10_BH")

plot_data_type$Contrast <- as.character(plot_data_type$Contrast)
plot_data_type$Contrast <- gsub("^DE1/|\\.csv$", "", plot_data_type$Contrast)
plot_data_type$Contrast <- gsub("zona1/", "ZoneSAN_", plot_data_type$Contrast)
plot_data_type$Contrast <- gsub("zona2/", "ZoneAVN_", plot_data_type$Contrast)
plot_data_type$Contrast <- gsub("zona3/", "ZonePF_", plot_data_type$Contrast)
plot_data_type$Contrast <- gsub("contrast_", "", plot_data_type$Contrast)
plot_data_type$Contrast <- factor(plot_data_type$Contrast, levels = c( 
  "ZonePF_Cardiomyocytes_vs_ConductionSystem", 
  "ZonePF_Cardiomyocytes_vs_CoronarySMC",       "ZonePF_Cardiomyocytes_vs_Endocardial" ,     
  "ZonePF_Cardiomyocytes_vs_Epicardial"  ,      "ZonePF_Cardiomyocytes_vs_Fibroblasts"  ,    
  "ZonePF_ConductionSystem_vs_CoronarySMC",     "ZonePF_ConductionSystem_vs_Endocardial" ,   
  "ZonePF_ConductionSystem_vs_Epicardial"  ,    "ZonePF_ConductionSystem_vs_Fibroblasts"  ,  
  "ZonePF_CoronarySMC_vs_Endocardial"       ,   "ZonePF_CoronarySMC_vs_Epicardial"         , 
  "ZonePF_CoronarySMC_vs_Fibroblasts"        ,  "ZonePF_Endocardial_vs_Epicardial"          ,
  "ZonePF_Endocardial_vs_Fibroblasts"         , "ZonePF_Epicardial_vs_Fibroblasts" ,
  "ZoneAVN_Cardiomyocytes_vs_ConductionSystem" ,"ZoneAVN_Cardiomyocytes_vs_CoronarySMC",     
  "ZoneAVN_Cardiomyocytes_vs_Endocardial" ,     "ZoneAVN_Cardiomyocytes_vs_Epicardial"  ,    
  "ZoneAVN_Cardiomyocytes_vs_Fibroblasts"  ,    "ZoneAVN_ConductionSystem_vs_CoronarySMC",   
  "ZoneAVN_ConductionSystem_vs_Endocardial" ,   "ZoneAVN_ConductionSystem_vs_Epicardial"  ,  
  "ZoneAVN_ConductionSystem_vs_Fibroblasts"  ,  "ZoneAVN_CoronarySMC_vs_Endocardial"       , 
  "ZoneAVN_CoronarySMC_vs_Epicardial"      ,    "ZoneAVN_CoronarySMC_vs_Fibroblasts"        ,
  "ZoneAVN_Endocardial_vs_Epicardial"       ,   "ZoneAVN_Endocardial_vs_Fibroblasts"        ,
  "ZoneAVN_Epicardial_vs_Fibroblasts"        ,           
  "ZoneSAN_ConductionSystem_vs_Cardiomyocytes", "ZoneSAN_ConductionSystem_vs_CoronarySMC"   ,
  "ZoneSAN_ConductionSystem_vs_Endocardial"    ,"ZoneSAN_ConductionSystem_vs_Epicardial"    ,
  "ZoneSAN_ConductionSystem_vs_Fibroblasts",    "ZoneSAN_CoronarySMC_vs_Cardiomyocytes"     ,
  "ZoneSAN_CoronarySMC_vs_Endocardial"      ,   "ZoneSAN_CoronarySMC_vs_Epicardial"         ,
  "ZoneSAN_CoronarySMC_vs_Fibroblasts"       ,  "ZoneSAN_Endocardial_vs_Cardiomyocytes"     ,
  "ZoneSAN_Endocardial_vs_Epicardial"         , "ZoneSAN_Endocardial_vs_Fibroblasts"        ,
  "ZoneSAN_Epicardial_vs_Cardiomyocytes"       ,"ZoneSAN_Epicardial_vs_Fibroblasts"         ,
  "ZoneSAN_Fibroblasts_vs_Cardiomyocytes" ))

plot_data_type$NES <- ifelse(plot_data_type$neg_log10_BH < -log10(0.05), NA, plot_data_type$NES)
plot_data_type = plot_data_type[order(plot_data_type$Contrast),]

ggplot(plot_data_type, aes(x = Set, y = Contrast)) +
  geom_point(aes(size = neg_log10_BH, color = NES, fill = NES), shape = 21,  color = "black", stroke = 0.5) +
  scale_size(name = "-log10(FDR)") + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, name = "NES", na.value = "grey") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, name = "NES", na.value = "gray") +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


color_zone_long <- melt(color_zone)
colnames(color_zone_long) <- c("Contrast", "Set", "NES")
size_zone_long <- melt(size_zone)
colnames(size_zone_long) <- c("Contrast", "Set", "neg_log10_BH")
plot_data_zone <- merge(color_zone_long, size_zone_long, by = c("Contrast", "Set"))
colnames(plot_data_zone) <- c("Contrast", "Set", "NES",  "neg_log10_BH")

plot_data_zone$Contrast <- as.character(plot_data_zone$Contrast)
plot_data_zone$Contrast <- gsub("^contrast_", "", gsub("\\.csv$", "", basename(plot_data_zone$Contrast)))

invert_format <- function(s) {
  parts <- unlist(strsplit(s, "_"))
  # Obtener el último elemento
  last_part <- parts[length(parts)]
  # Obtener el resto
  rest_parts <- parts[-length(parts)]
  # Reorganizar el string
  inverted <- paste(c(last_part, paste(rest_parts, collapse = "_")), collapse = "_")
  return(inverted)
}
plot_data_zone$Contrast <- sapply(plot_data_zone$Contrast, invert_format)
plot_data_zone$NES <- ifelse(plot_data_zone$neg_log10_BH < -log10(0.05), NA, plot_data_zone$NES)
plot_data_zone$Contrast = factor(plot_data_zone$Contrast, levels = rev(c(
  "Cardiomyocytes_ZoneSAN_vs_ZoneAVN", "Cardiomyocytes_ZoneSAN_vs_ZonePF", "Cardiomyocytes_ZoneAVN_vs_ZonePF",
  "ConductionSystem_ZoneSAN_vs_ZoneAVN", "ConductionSystem_ZoneSAN_vs_ZonePF", "ConductionSystem_ZoneAVN_vs_ZonePF",
  "CoronarySMC_ZoneSAN_vs_ZoneAVN", "CoronarySMC_ZoneSAN_vs_ZonePF", "CoronarySMC_ZoneAVN_vs_ZonePF",
  "Endocardial_ZoneSAN_vs_ZoneAVN", "Endocardial_ZoneSAN_vs_ZonePF", "Endocardial_ZoneAVN_vs_ZonePF",
  "Epicardial_ZoneSAN_vs_ZoneAVN", "Epicardial_ZoneSAN_vs_ZonePF", "Epicardial_ZoneAVN_vs_ZonePF",
  "Fibroblasts_ZoneSAN_vs_ZoneAVN", "Fibroblasts_ZoneSAN_vs_ZonePF", "Fibroblasts_ZoneAVN_vs_ZonePF"
)))

ggplot(plot_data_zone, aes(x = Set, y = Contrast)) +
  geom_point(aes(size = neg_log10_BH, color = NES, fill = NES), shape = 21, color = "black", stroke = 0.5) +
  scale_size(name = "-log10(FDR)") + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, name = "NES", na.value = "gray") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, name = "NES", na.value = "gray") +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


