library(stringr)
library(limma)
library(ggplot2)
library(viridis)
library(ggrepel)
load("~/Documents/TFM/Rdatas/expression_df.RData")
mito_genes <-c("mt-Co2","mt-Cytb","mt-Co3","mt-Nd2","mt-Nd1","mt-Nd4","mt-Atp6" ,"mt-Co1","mt-Nd5",
               "mt-Nd3" ,"mt-Nd4l" ,"mt-Atp8", "mt-Nd6" )

### filtering of genes present in at least 10% cells pf at least 1 cell_group ###
##################################################################

{X <- as.data.frame(t(expression_df[, !colnames(expression_df) %in% c("Cell_type", "Zone","cell_group")]))

metadata=data.frame(Cell_group=expression_df$cell_group)
metadata$Cell_type=str_sub(metadata$Cell_group, end = -2)
metadata$Zone=str_sub(metadata$Cell_group, start = -1)

metadata$Zone[which(metadata$Zone==1)]="SAN"
metadata$Zone[which(metadata$Zone==2)]="AVN"
metadata$Zone[which(metadata$Zone==3)]="PF"
metadata$Tomm20=expression_df[,"Tomm20"]

auxiliar_design=model.matrix(~0+Cell_group,data=metadata)
gene_detections=X[,1:ncol(auxiliar_design)]
colnames(gene_detections)=colnames(auxiliar_design)
fraction_detections=gene_detections

count=function(x){length(which(x>0))}

for(cell_group in  colnames(gene_detections)){
  print(cell_group)
  chunk=X[,which(auxiliar_design[,cell_group]==1)]
  gene_detections[,cell_group]=apply(chunk,1,count)
  fraction_detections[,cell_group]=gene_detections[,cell_group]/ncol(chunk)
}

gene_detections$max=apply(gene_detections,1,max)
fraction_detections$max=apply(fraction_detections,1,max)

keep_genes = rownames(fraction_detections)[fraction_detections$max>0.1]
#write.csv(keep_genes, "/Users/.../TFM/data/keep_genes.csv")
genes = colnames(expression_df)}

##################################################################
###### assembly factors and subunits

{factores_I = c("Ndufv1", "Ndufv2", "Ndufs1", "Ndufa2", "Ndufa6", "Ndufa7", "Ndufa12", "Ndufs4", "Ndufs6", "Ndufv3", "Ndufa5", "Ndufs2", "Ndufs3", "Ndufs7", "Ndufs8",
                "Ndufa3", "C3orf60", "NdufAf4", "C6orf66", "Ndufaf6", "C8orf38", "Ndufaf3", "Ndufaf4", "Ndufaf5", "Ndufaf6", "Ndufaf7", "Nubpl", "Ind1", "Timmdc1", "C3orf1", "Ndufaf1", "Cia30", "Ecsit", "Acad9", "Tmem126b", "Tmem186", "Coa1", "Mitrac15", "Foxred1", "Atp5sl", "Tmem70", "Dmac1", "Tmem261", "Ndufaf2", "Ndufa12l")
  factores_I = factores_I[factores_I %in% genes]
  factores_Nd1 = c("Ndufa3", "Ndufa8", "Ndufa13")
  factores_Nd1 = factores_Nd1[factores_Nd1 %in% genes]
  factores_Nd2364L = c("Ndufc1", "Ndufc2", "Ndufa1", "Ndufa10", "Ndufs5")
  factores_Nd2364L = factores_Nd2364L[factores_Nd2364L %in% genes]
  factores_Nd4L = c("Ndufb1", "Ndufb4", "Ndufb5", "Ndufb6", "Ndufb10", "Ndufb11")
  factores_Nd4L = factores_Nd4L[factores_Nd4L %in% genes]
  factores_Nd5 = c("Ndufb2", "Ndufb3", "Ndufb7", "Ndufb8", "Ndufb9", "Ndufab1")
  factores_Nd5 = factores_Nd5[factores_Nd5 %in% genes]
  
  #complex 2
  factores_II = c("Sdha", "Sdh1", "Sdhb", "Sdh2", "Sdhc", "Sdh3", "Sdhd", "Sdh4")
  factores_II = factores_II[factores_II %in% genes]
  
  #complex 3
  factores_III = c("Uqcrc1", "Core1", "Bcs1l", "Uqcrc2", "Core2",	"Lyrm7", "Mzm1l","Cyc1",	"Uqcc1", "Uqcrfs1", "Risp", "Uqcc2", "Mt-Cyb",	"Uqcc3","Uqcrb", "Ttc19", "Uqcrq",	"C12orf73", "Brawnin", "Uqcrh", "C16orf91", "Uqcc4","Uqcr10",	"Ociad1", "Uqcr11",	"Ociad2", "Stmp1", "Smim4", "Sfxn1")
  factores_III = factores_III[factores_III %in% genes]
  
  #complex 4
  factores_Co1 = c("Taco1", "Mitrac12", "Coa3", "Ccdc56", "C12orf62", "Cox14", "Mitrac7", "Smim20", "Cmc1", "Coa1", "Mitrac15", "Surf1", "Cox10", "Cox11", "Cox15", "Cox17", "Cox19")
  factores_Co1 = factores_Co1[factores_Co1 %in% genes]
  factores_Co2 = c("Cox20", "Cox18", "Coa", "Sco1", "Sco2", "Cox17", "Cox16", "Pet100", "Pet177", "Mr1s")
  factores_Co2 = factores_Co2[factores_Co2 %in% genes]
  factores_Co = c("Cox4i1", "Cox5a", "Cox5b", "Cox6c", "Cox7c", "Cox8a", "Cox7b", "Cox6a1", "Cox6b1", "Cox7a2", "Ndufa4")
  factores_Co = factores_Co[factores_Co %in% genes]
  
  #complex 5
  factores_Atp = c("Atp5a1", "Atp5b", "Atp5c1", "Atp5d", "Atp5e", "Atpaf1", "Atp11", "Atpg1", "Atpaf2", "Atp12", "Atpg2", "Atpg3", "Atp5f1", "Atph", "Atp5j", "Atp5o", "Atp5i", "Atp5l", "Atpj2", "Tmem70", "C14orf2", "Usmg5")
  factores_Atp = factores_Atp[factores_Atp %in% genes]}

factores =c(factores_I, factores_Nd1, factores_Nd2364L, factores_Nd4L, factores_Nd5, factores_II, factores_III, factores_Co1, factores_Co2, factores_Co)
factores = unique(factores)
rm(factores_I, factores_Nd1, factores_Nd2364L, factores_Nd4L, factores_Nd5, factores_II, factores_III, factores_Co1, factores_Co2, factores_Co)

##################################################################
########## LMFit

genes_good =c(keep_genes, factores)
write.csv(genes_good, "/Users/.../TFM/data/keep_genes.csv")

expression_df_filtered = expression_df[,colnames(expression_df) %in% c(genes_good, "cell_group")]
expression_df_filtered$cell_group <- factor(expression_df_filtered$cell_group, 
                                   levels = c("Cardiomyocytes1", 
                                              "Fibroblasts1", 
                                              "CoronarySMC1", 
                                              "Epicardial1", 
                                              "ConductionSystem1", 
                                              "Endocardial1",
                                              "Cardiomyocytes2", 
                                              "Fibroblasts2", 
                                              "CoronarySMC2", 
                                              "Epicardial2", 
                                              "ConductionSystem2", 
                                              "Endocardial2",
                                              "Cardiomyocytes3", 
                                              "Fibroblasts3", 
                                              "CoronarySMC3", 
                                              "Epicardial3", 
                                              "ConductionSystem3", 
                                              "Endocardial3"))

y <- model.matrix(~ cell_group, data =expression_df_filtered)
X <- as.data.frame(t(expression_df_filtered[, !colnames(expression_df_filtered) %in% c("cell_group")]))

fit <- lmFit(X, design = y)
fit <- eBayes(fit)

##################################################################
########## contrasts


get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


########## contrasts of zones in same cell type
cell_types = c("Cardiomyocytes", "ConductionSystem", "CoronarySMC", "Endocardial","Epicardial","Fibroblasts") 
vec_indices <- data.frame("ident" = colnames(y), "index" = c(1:18))
output_dir <- "/Users/.../TFM/data/DE/DE2"

list_data_z = list()
list_density_z = list()
for (celltype_n in 1:6){
  celltype = cell_types[celltype_n]
  barplot_data_z = as.data.frame(matrix(0, nrow = 3, ncol = 3))
  rownames(barplot_data_z) = c("ZoneSAN", "ZoneAVN", "ZonePF")
  colnames(barplot_data_z) = c("ZoneSAN", "ZoneAVN", "ZonePF")
  density_data_z = as.data.frame(matrix(0, nrow = 1000, ncol = 6))
  colnames(density_data_z) = c("ZoneSAN_ZoneAVN", "ZoneSAN_ZonePF", 
                               "ZoneAVN_ZoneSAN", "ZoneAVN_ZonePF",
                               "ZonePF_ZoneSAN", "ZonePF_ZoneAVN")
  
  for (i in 1:2) {
    if (i==1){
      ident1 = "ZoneSAN"
    }else if(i == 2){
      ident1 = "ZoneAVN"
    }else if(i == 3){
      ident1 = "ZonePF"
    }else{
      return("Error")
    }
    for (j in (i + 1):3) {
      if (j==1){
        ident2 = "ZoneSAN"
      }else if(j == 2){
        ident2 = "ZoneAVN"
      }else if(j == 3){
        ident2 = "ZonePF"
      }else{
        return("Error")
      }
      
      plot_path <- paste0(output_dir, "/", celltype)
      compare1 = paste0("cell_group",celltype, i)
      compare2 = paste0("cell_group",celltype, j)
      if ((celltype == 1) && (i == 1)){
        vec = rep(0, 18)
        vec[vec_indices$index[vec_indices$ident == compare2]] = 1 
        
        c = contrasts.fit(fit, vec)
        c =eBayes(c)
        result = topTable(c, number = Inf, sort.by = "none", adjust.method = "fdr")
        
        #write.csv(result,paste0(plot_path, "/contrast_", ident2,"_vs_", ident1, ".csv"))
        
      }else{
        vec = rep(0, 18)
        vec[vec_indices$index[vec_indices$ident == compare1]] = 1 
        vec[vec_indices$index[vec_indices$ident == compare2]] = -1 
        
        c = contrasts.fit(fit, vec)
        c =eBayes(c)
        result = topTable(c, number = Inf, sort.by = "none", adjust.method = "fdr")
        
        #write.csv(result, paste0(plot_path, "/contrast_", ident1,"_vs_", ident2,  "_", celltype, ".csv"))
      }
      
      result_clean <- result[result$adj.P.Val > 0, ]
      result_clean$density = get_density(x=result_clean$logFC, y = -log10(result_clean$adj.P.Val),h=1, n=100)
      
      # Volcano plot
      volcano_plot <- ggplot(result_clean, aes(x = logFC, y = -log10(adj.P.Val))) +
        geom_point(aes(color = density), size = 0.7) + scale_color_viridis() +
        geom_text_repel(data = head(result_clean[order(abs(result_clean$logFC)),], 10),
                        aes(label = rownames(head(result_clean[order(abs(result_clean$logFC)),], 10))),
                        size = 3.3, max.overlaps = 20, colour = "red", segment.colour="green") +
        theme_minimal() + labs(title = paste0(ident1, " vs ", ident2, " in ", celltype))
      
      #ggsave(filename = paste0(plot_path, "/volcano_", ident1, "_vs_", ident2, "_", celltype, ".png"), plot = volcano_plot)
      
      barplot_data_z[ident1, ident2] = nrow(result[which((result$logFC > 0.2) & (result$adj.P.Val < 0.05)),])
      barplot_data_z[ident2, ident1] = nrow(result[which((result$logFC > 0.2) & (result$adj.P.Val < 0.05)),])
      
      aux = result[order(result$adj.P.Val), "logFC"]
      colname1 = paste0(ident1, "_", ident2)
      colname2 = paste0(ident2, "_", ident1)
      density_data_z[,colname1] = head(aux, n = 1000)
      density_data_z[,colname2] = -head(aux, n = 1000)
      
      #print(paste0(celltype, ": max logFC ", colname1, ": ", max(result$logFC)))
    }
  }
  list_data_z[[celltype_n]] = barplot_data_z
  list_density_z[[celltype_n]] = density_data_z
}

########## contrasts of cell types in same zone

output_dir <- "/Users/.../TFM/data/DE/DE1"
zona1_dir <- file.path(output_dir, "zona1")
zona2_dir <- file.path(output_dir, "zona2")
zona3_dir <- file.path(output_dir, "zona3")

list_data = list()
list_density = list()
for (index in 1:3){
  if (index == 1){
    zone="Zone SAN"
    plot_path = zona1_dir
  }else if(index == 2){
    zone="Zone AVN"
    plot_path = zona2_dir
  }else if(index == 3){
    zone="Zone PF"
    plot_path = zona3_dir
  }else{
    print("Error")
  }
  barplot_data = as.data.frame(matrix(0, nrow = 6, ncol = 6))
  rownames(barplot_data) = cell_types
  colnames(barplot_data) = cell_types
  
  density_data = as.data.frame(matrix(0, nrow = 1000, ncol = 50))
  colnames(density_data) = c("Cardiomyocytes_ConductionSystem", "Cardiomyocytes_CoronarySMC", "Cardiomyocytes_Endocardial", "Cardiomyocytes_Epicardial", "Cardiomyocytes_Fibroblasts",
                             "ConductionSystem_Cardiomyocytes", "ConductionSystem_CoronarySMC", "ConductionSystem_Endocardial", "ConductionSystem_Epicardial", "ConductionSystem_Fibroblasts",
                             "CoronarySMC_Cardiomyocytes", "CoronarySMC_ConductionSystem", "CoronarySMC_Endocardial", "CoronarySMC_Epicardial", "CoronarySMC_Fibroblasts",
                             "Endocardial_Cardiomyocytes", "Endocardial_ConductionSystem", "Endocardial_CoronarySMC", "Endocardial_Epicardial", "Endocardial_Fibroblasts",
                             "Epicardial_Cardiomyocytes", "Epicardial_ConductionSystem", "Epicardial_CoronarySMC", "Epicardial_Endocardial", "Epicardial_Fibroblasts",
                             "Fibroblasts_Cardiomyocytes", "Fibroblasts_ConductionSystem", "Fibroblasts_CoronarySMC", "Fibroblasts_Epicardial", "Fibroblasts_Endocardial")
  
  for (i in 1:(length(cell_types) - 1)) {
    for (j in (i + 1):length(cell_types)) {
      celltype1 <- cell_types[i]
      celltype2 <- cell_types[j]
      
      compare1 = paste0("cell_group",celltype1, index)
      compare2 = paste0("cell_group",celltype2, index)
      if ((i == 1) && (index == 1)){
        vec = rep(0, 18)
        vec[vec_indices$index[vec_indices$ident == compare2]] = 1 
        
        c = contrasts.fit(fit, vec)
        c =eBayes(c)
        result = topTable(c, number = Inf, sort.by = "none", adjust.method = "fdr")
        
        write.csv(result, paste0(plot_path, "/contrast_", celltype2,"_vs_", celltype1, ".csv"))
        
      }else{
        vec = rep(0, 18)
        vec[vec_indices$index[vec_indices$ident == compare1]] = 1 
        vec[vec_indices$index[vec_indices$ident == compare2]] = -1 
        
        c = contrasts.fit(fit, vec)
        c =eBayes(c)
        result = topTable(c, number = Inf, sort.by = "none", adjust.method = "fdr")
        
        write.csv(result, paste0(plot_path, "/contrast_", celltype1,"_vs_", celltype2, ".csv"))
      }
      
      result_clean <- result[result$adj.P.Val > 0, ]
      result_clean$density = get_density(x=result_clean$logFC, y = -log10(result_clean$adj.P.Val),h=1, n=100)
      
      # Volcano plot
      volcano_plot <- ggplot(result_clean, aes(x = logFC, y = -log10(adj.P.Val))) +
        geom_point(aes(color = density), size = 0.7) + scale_color_viridis() +
        geom_text_repel(data = head(result_clean[order(abs(result_clean$logFC)),], 10),
                        aes(label = rownames(head(result_clean[order(abs(result_clean$logFC)),], 10))),
                        size = 3.3, max.overlaps = 20, colour = "red", segment.colour="green") +
        theme_minimal() + labs(title = paste0(ident1, " vs ", ident2, " in ", zone))
      
      #ggsave(filename = file.path(plot_path, paste0("volcano_", ident1, "_vs_", ident2, ".png")),plot = volcano_plot)
      
      
      barplot_data[celltype1, celltype2] = nrow(result[which((result$logFC > 0.2) & (result$adj.P.Val < 0.05)),])
      barplot_data[celltype2, celltype1] = nrow(result[which((result$logFC > 0.2) & (result$adj.P.Val < 0.05)),])
    
      aux = result[order(result$adj.P.Val), "logFC"]
      colname1 = paste0(celltype1, "_", celltype2)
      colname2 = paste0(celltype2, "_", celltype1)
      density_data[,colname1] = head(aux, n = 1000)
      density_data[,colname2] = -head(aux, n = 1000)
      
      #print(paste0(zone, ": max logFC ", colname1, ": ", max(result$logFC)))
    }
  }
  
  list_data[[index]] = barplot_data
  list_density[[index]] = density_data
}


##################################################################
######## barplots and density plots

library(ggplot2)
library(tidyr)
library(tibble)
library(gridExtra)

zones_list = c("ZoneSAN", "ZoneAVN", "ZonePF")
preparar_datos_z <- function(tabla, i) {
  tabla <- tabla[i, ]
  tabla_larga <- as.data.frame(tabla) %>%
    rownames_to_column("Zone1") %>%
    pivot_longer(cols = zones_list, names_to = "Zone2", values_to = "valor")
  return(tabla_larga)
}

for (celltype_n in 1:6){
  lista_graficos <- list()
  for(index in 1:3){
    barplot_data_long <- preparar_datos_z(list_data_z[[celltype_n]], index)
    barplot_data_long$Zone2 = factor(barplot_data_long$Zone2, levels=zones_list)
    
    grafico = ggplot(barplot_data_long, aes(x = Zone2, y = valor, fill = Zone2)) +
      geom_col(position = "dodge", color = c("steelblue", "chocolate1", "darkolivegreen"), fill = c("steelblue", "chocolate1", "darkolivegreen")) +
      geom_text(aes(label = valor), 
                position = position_dodge(width = 0.9),
                size = 10) +
      labs(x = "", y = "Number of DE Genes") +
      theme_minimal() + 
      ggtitle(zones_list[index]) + 
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 15),          
        axis.title.y = element_text(size = 15),       
        axis.text.x = element_text(size = 13),           
        axis.text.y = element_text(size = 13),    
        legend.position = 'none'
      )
    
    lista_graficos[[index]] <- grafico
  }
  
  png(paste0("/Users/.../TFM/imgs/", cell_types[celltype_n],"_bars_zones.png"), width = 800, height = 800, units = "px")
  grid.arrange(grobs = lista_graficos, nrow = 3)
  dev.off()
  
  
  lista_graficos <- list()
  for(index in 1:3){
    if(index ==1){
      data_plot_den = list_density_z[[celltype_n]][, 1:2]
      colors_pl = c("chocolate1", "darkolivegreen")
    }else if(index ==2){
      data_plot_den = list_density_z[[celltype_n]][, 3:4]
      colors_pl = c("steelblue", "darkolivegreen")
    }else if(index ==3){
      data_plot_den = list_density_z[[celltype_n]][, 5:6]
      colors_pl = c("steelblue", "chocolate1")
    }else{
      print("error")
    }
    data_plot_den_long <- data.frame(
      value = c(abs(data_plot_den[,1]), abs(data_plot_den[,2])),
      group = rep(c("Group 1", "Group 2"), each = nrow(data_plot_den))
    )
    grafico = ggplot(data_plot_den_long, aes(x = value, color = group)) +
      geom_density(size = 1) +
      scale_color_manual(values = colors_pl) +  
      labs(x = "abs(logFC)", y = "Density") +
      theme_minimal() + ylim(c(0, 10)) + xlim(c(0, 1.5)) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 15),          
        axis.title.y = element_text(size = 15),       
        axis.text.x = element_text(size = 13),           
        axis.text.y = element_text(size = 13),    
        legend.position = 'none'
      ) 
    
    lista_graficos[[index]] <- grafico
    
  }
  png(paste0("/Users/.../TFM/imgs/", cell_types[celltype_n],"_bars_zones_den.png"), width = 600, height = 800, units = "px")
  grid.arrange(grobs = lista_graficos, nrow = 3)
  dev.off()
  
}

#######################################

zones_list = c("ZoneSAN", "ZoneAVN", "ZonePF")
preparar_datos <- function(tabla, i) {
  tabla <- tabla[i, ]
  tabla_larga <- as.data.frame(tabla) %>%
    rownames_to_column("CellType1") %>%
    pivot_longer(cols = cell_types, names_to = "CellType2", values_to = "valor")
  return(tabla_larga)
}

for (zone_n in 1:3){
  lista_graficos <- list()
  for(index in 1:6){
    barplot_data_long <- preparar_datos(list_data[[zone_n]], index)
    
    grafico = ggplot(barplot_data_long, aes(x = CellType2, y = valor, fill = CellType2)) +
      geom_col(position = "dodge") +
      geom_text(aes(label = valor), 
                position = position_dodge(width = 0.9),
                size = 5) +
      labs(x = "", y = "Nº of DE genes") +
      theme_minimal() + ggtitle(cell_types[index]) + 
      theme(plot.title = element_text(hjust = 0.5), legend.position='none')
    
    lista_graficos[[index]] <- grafico
  }
  
  png(paste0("/Users/.../TFM/imgs/", zones_list[zone_n], "bars_celltype.png"), width = 500, height = 1000, units = "px")
  grid.arrange(grobs = lista_graficos, nrow = 6)
  dev.off()
  
  lista_graficos <- list()
  for(index in 1:6){
    if(index ==1){
      data_plot_den = list_density[[zone_n]][, 1:5]
      colors_pl = c("goldenrod", "green3", "cyan3", "dodgerblue", "hotpink")
    }else if(index ==2){
      data_plot_den = list_density[[zone_n]][, 6:10]
      colors_pl = c("salmon", "green3", "cyan3", "dodgerblue", "hotpink")
    }else if(index ==3){
      data_plot_den = list_density[[zone_n]][, 11:15]
      colors_pl = c("salmon","goldenrod", "cyan3", "dodgerblue", "hotpink")
    }else if(index ==4){
      data_plot_den = list_density[[zone_n]][, 16:20]
      colors_pl = c("salmon", "goldenrod", "green3", "dodgerblue", "hotpink")
    }else if(index ==5){
      data_plot_den = list_density[[zone_n]][, 21:25]
      colors_pl = c("salmon", "goldenrod", "green3","cyan3", "hotpink")
    }else if(index ==6){
      data_plot_den = list_density[[zone_n]][, 26:30]
      colors_pl = c("salmon","goldenrod", "green3","cyan3", "dodgerblue")
    }else{
      print("error")
    }
    data_plot_den_long <- data.frame(
      value = c(abs(data_plot_den[,1]), abs(data_plot_den[,2]), abs(data_plot_den[,3]), abs(data_plot_den[,4]), abs(data_plot_den[,5])),
      group = rep(c(colnames(data_plot_den)[1], colnames(data_plot_den)[2], colnames(data_plot_den)[3], colnames(data_plot_den)[4], colnames(data_plot_den)[5]), each = nrow(data_plot_den))
    )
    grafico <- ggplot(data_plot_den_long, aes(x = value, color = group)) +
      geom_density(size = 1) +
      scale_color_manual(values = colors_pl) +  
      labs(x = "abs(logFC)", y = "Density") +
      theme_minimal() +
      ylim(c(0, 10)) +
      xlim(c(0, 1.5)) +
      theme(legend.position = 'none')
    
    lista_graficos[[index]] <- grafico
  }
  png(paste0("/Users/.../TFM/imgs/", zones_list[zone_n], "bars_celltype_den.png"), width = 300, height = 1000, units = "px")
  grid.arrange(grobs = lista_graficos, nrow = 6)
  dev.off()
  
}






