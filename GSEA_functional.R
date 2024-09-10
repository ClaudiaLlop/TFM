library(fgsea)
library(msigdbr)
library(biomaRt)

########################################
##### Correlation average mito gene expression vs all rest
load("~/Documents/TFM/Rdatas/expression_df.RData")
mito_genes <-c("mt-Co2","mt-Cytb","mt-Co3","mt-Nd2","mt-Nd1","mt-Nd4","mt-Atp6" ,"mt-Co1","mt-Nd5",
               "mt-Nd3" ,"mt-Nd4l" ,"mt-Atp8", "mt-Nd6" )
keep_genes = read.csv("/Users/claudia/Documents/TFM/data/keep_genes.csv", row.names=1)

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

vec_corr_avg_mt = r_matrix[1,]
names(vec_corr_avg_mt) = colnames(r_matrix)


tab = data.frame(msigdbr(species = "Mus musculus", category = "H"))
pathways_Hallmark_msigDB = split(x = tab$ensembl_gene, f = tab$gs_name)
mouse_mart <- useMart(host="https://feb2021.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl")

genes_ensembl <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                       filters = "mgi_symbol",
                       values = names(vec_corr_avg_mt),
                       mart = mouse_mart)
genes_ensembl_unique <- genes_ensembl[!duplicated(genes_ensembl$mgi_symbol), ]
genes_ensembl_unique=genes_ensembl_unique[order(genes_ensembl_unique$mgi_symbol),]
vec_corr_avg_mt=vec_corr_avg_mt[order(names(vec_corr_avg_mt))]
vec_corr_avg_mt=vec_corr_avg_mt[which(names(vec_corr_avg_mt) %in% genes_ensembl_unique$mgi_symbol)]
length(which(genes_ensembl_unique[,2]==names(vec_corr_avg_mt)))

names(vec_corr_avg_mt)=genes_ensembl_unique[,1]
vec_corr_avg_mt=vec_corr_avg_mt[order(-vec_corr_avg_mt)]

fgsea_results <- fgsea(pathways = pathways_Hallmark_msigDB, stats= vec_corr_avg_mt, eps=0)
res_hallmark=data.frame(fgsea_results)[which(fgsea_results$padj<0.05),1:7]
res_hallmark=res_hallmark[order(res_hallmark$NES),]
#colecciones=data.frame(msigdbr_collections())
#gs_cat       gs_subcat num_genesets
#1      C1                          299
#2      C2             CGP         3384
#3      C2              CP           29
#4      C2     CP:BIOCARTA          292
#5      C2         CP:KEGG          186
#6      C2          CP:PID          196
#7      C2     CP:REACTOME         1615
#8      C2 CP:WIKIPATHWAYS          664
#9      C3       MIR:MIRDB         2377
#10     C3  MIR:MIR_Legacy          221
#11     C3        TFT:GTRD          518
#12     C3  TFT:TFT_Legacy          610
#13     C4             CGN          427
#14     C4              CM          431
#15     C5           GO:BP         7658
#16     C5           GO:CC         1006
#17     C5           GO:MF         1738
#18     C5             HPO         5071
#19     C6                          189
#20     C7     IMMUNESIGDB         4872
#21     C7             VAX          347
#22     C8                          700
#23      H                           50



tab = data.frame(msigdbr(species = "Mus musculus", category = "C2",subcategory = "CP:KEGG"))
pathways_KEGG_msigDB = split(x = tab$ensembl_gene, f = tab$gs_name)
fgsea_results_KEGG <- fgsea(pathways = pathways_KEGG_msigDB, stats= vec_corr_avg_mt, eps=0)
res_KEGG=data.frame(fgsea_results_KEGG)[which(fgsea_results_KEGG$padj<0.05),1:7]
res_KEGG=res_KEGG[order(res_KEGG$NES),]



tab = data.frame(msigdbr(species = "Mus musculus", category = "C5",subcategory = "GO:BP"))
pathways_GO_msigDB = split(x = tab$ensembl_gene, f = tab$gs_name)
fgsea_results_GO <- fgsea(pathways = pathways_GO_msigDB, stats= vec_corr_avg_mt, eps=0)
res_GO=data.frame(fgsea_results_GO)[which(fgsea_results_GO$padj<0.05),1:7]
res_GO=res_GO[order(res_GO$NES),]


save(res_GO,file="results_enrichments_GO.Rdata")
save(res_KEGG,file="results_enrichments_KEGG.Rdata")
save(res_hallmark,file="results_enrichments_HALLMARK.Rdata")

##########################################
library(ggplot2)

lowest_NES <- res_GO[order(res_GO$NES), ][1:10, ]
highest_NES <- res_GO[order(res_GO$NES, decreasing = TRUE), ][1:10, ]
res_GO_filtered <- rbind(lowest_NES, highest_NES)
res_GO_filtered$NES_sign <- ifelse(res_GO_filtered$NES > 0, "Positive", "Negative")

res_GO_filtered$pathway = gsub("GOBP_", "", res_GO_filtered$pathway)

ggplot(res_GO_filtered, aes(x = NES, y = reorder(pathway, NES), color = NES_sign)) +
  geom_segment(aes(xend = 0, yend = reorder(pathway, NES)), size = 1, show.legend = FALSE) + 
  geom_point(aes(size = -log10(padj)), shape = 21, fill = "black", color = "white") + 
  scale_color_manual(values = c("Positive" = "red", "Negative" = "blue")) + 
  scale_size_continuous(name = "-log10(padj)", range = c(3, 10)) +
  labs(x = "NES", y = "Pathway", title = "GO Top Pathways") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 14, face = "bold"))

#############################################









library(biomaRt)
load("datafit.RData")
keep_genes = read.csv("keep_genes.csv", row.names=1)

cell_types <- c("Cardiomyocytes", "ConductionSystem", "CoronarySMC", 
                "Endocardial", "Epicardial", "Fibroblasts")
mito_genes <-c("mt-Co2","mt-Cytb","mt-Co3","mt-Nd2","mt-Nd1","mt-Nd4","mt-Atp6" ,"mt-Co1","mt-Nd5",
               "mt-Nd3" ,"mt-Nd4l" ,"mt-Atp8", "mt-Nd6" )
zones_list <- c("ZoneSAN", "ZoneAVN", "ZonePF")

genes = colnames(expression_df)
rm(expression_df)

mouse <- useMart(host="https://feb2021.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl")
human <- useMart(host="https://feb2021.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

{
  factores_I_human = toupper(c("Ndufv1", "Ndufv2", "Ndufs1", "Ndufa2", "Ndufa6", "Ndufa7", "Ndufa12", "Ndufs4", "Ndufs6", "Ndufv3", "Ndufa5", "Ndufs2", "Ndufs3", "Ndufs7", "Ndufs8",
                "Ndufa3", "C3orf60", "NdufAf4", "C6orf66", "Ndufaf6", "C8orf38", "Ndufaf3", "Ndufaf4", "Ndufaf5", "Ndufaf6", "Ndufaf7", "Nubpl", "Ind1", "Timmdc1", "C3orf1", "Ndufaf1", "Cia30", "Ecsit", "Acad9", "Tmem126b", "Tmem186", "Coa1", "Mitrac15", "Foxred1", "Atp5sl", "Tmem70", "Dmac1", "Tmem261", "Ndufaf2", "Ndufa12l"))
                factores_I_human[which(factores_I_human=="C8ORF38")]="C8orf38"
                factores_I_human[which(factores_I_human=="C3ORF60")]="C3orf60"

                factores_I = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values=factores_I_human , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)

  factores_I = factores_I$MGI.symbol[factores_I$MGI.symbol %in% genes]
  
  #complejo 2
  factores_II_human = toupper(c("Sdha", "Sdh1", "Sdhb", "Sdh2", "Sdhc", "Sdh3", "Sdhd", "Sdh4"))
  
  factores_II = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values=factores_II_human , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  factores_II = factores_II$MGI.symbol[factores_II$MGI.symbol %in% genes]
  
  #complejo 3
  factores_III_human = toupper(c("Uqcrc1", "Core1", "Bcs1l", "Uqcrc2", "Core2",	"Lyrm7", "Mzm1l","Cyc1",	"Uqcc1", "Uqcrfs1", "Risp", "Uqcc2", "Uqcc3","Uqcrb", "Ttc19", "Uqcrq",	"C12orf73", "Brawnin", "Uqcrh", "C16orf91", "Uqcc4","Uqcr10",	"Ociad1", "Uqcr11",	"Ociad2", "Stmp1", "Smim4", "Sfxn1"))
  factores_I_human[which(factores_III_human=="C12ORF73")]="C12orf73"
  factores_I_human[which(factores_III_human=="C16ORF91")]="C16orf91"

  
  ## Quito Mt-Cyb, es mitocondrial
  factores_III = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values=factores_III_human , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  factores_III = factores_III$MGI.symbol[factores_III$MGI.symbol %in% genes]

  #complejo 4
  factores_IV_human = toupper(c("Cox4i1", "Cox5a", "Cox5b", "Cox6c", "Cox7c", "Cox8a", "Cox7b", "Cox6a1", "Cox6b1", "Cox7a2", "Ndufa4"))
  factores_IV = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values=factores_IV_human , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  factores_IV = factores_IV$MGI.symbol[factores_IV$MGI.symbol %in% genes]

  #complejo 5
  factores_V_human = toupper(c("Atp5a1", "Atp5b", "Atp5c1", "Atp5d", "Atp5e", "Atpaf1", "Atp11", "Atpg1", "Atpaf2", "Atp12", "Atpg2", "Atpg3", "Atp5f1", "Atph", "Atp5j", "Atp5o", "Atp5i", "Atp5l", "Atpj2", "Tmem70", "C14orf2", "Usmg5"))
  factores_V_human[which(factores_V_human=="C14ORF2")]="C14orf2"

  
  factores_V = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values=factores_V_human , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  factores_V = factores_V$MGI.symbol[factores_V$MGI.symbol %in% genes]

  factores =c(factores_I, factores_II, factores_III, factores_IV, factores_V)
  
  
  factores_I_alt = (c("Ndufv1", "Ndufv2", "Ndufs1", "Ndufa2", "Ndufa6", "Ndufa7", "Ndufa12", "Ndufs4", "Ndufs6", "Ndufv3", "Ndufa5", "Ndufs2", "Ndufs3", "Ndufs7", "Ndufs8",
                "Ndufa3", "C3orf60", "NdufAf4", "C6orf66", "Ndufaf6", "C8orf38", "Ndufaf3", "Ndufaf4", "Ndufaf5", "Ndufaf6", "Ndufaf7", "Nubpl", "Ind1", "Timmdc1", "C3orf1", "Ndufaf1", "Cia30", "Ecsit", "Acad9", "Tmem126b", "Tmem186", "Coa1", "Mitrac15", "Foxred1", "Atp5sl", "Tmem70", "Dmac1", "Tmem261", "Ndufaf2", "Ndufa12l"))
                
                factores_I_alt=factores_I_alt[factores_I_alt %in% genes]
                
  factores_II_alt = c("Sdha", "Sdh1", "Sdhb", "Sdh2", "Sdhc", "Sdh3", "Sdhd", "Sdh4")
  
  factores_II_alt=factores_II_alt[factores_II_alt %in% genes]
  
  factores_III_alt = c("Uqcrc1", "Core1", "Bcs1l", "Uqcrc2", "Core2",    "Lyrm7", "Mzm1l","Cyc1",    "Uqcc1", "Uqcrfs1", "Risp", "Uqcc2", "Uqcc3","Uqcrb", "Ttc19", "Uqcrq",    "C12orf73", "Brawnin", "Uqcrh", "C16orf91", "Uqcc4","Uqcr10",    "Ociad1", "Uqcr11",    "Ociad2", "Stmp1", "Smim4", "Sfxn1")
  factores_III_alt=factores_III_alt[factores_III_alt %in% genes]
  
  factores_IV_alt = c("Cox4i1", "Cox5a", "Cox5b", "Cox6c", "Cox7c", "Cox8a", "Cox7b", "Cox6a1", "Cox6b1", "Cox7a2", "Ndufa4")
  factores_IV_alt=factores_IV_alt[factores_IV_alt %in% genes]
  
  factores_V_alt = c("Atp5a1", "Atp5b", "Atp5c1", "Atp5d", "Atp5e", "Atpaf1", "Atp11", "Atpg1", "Atpaf2", "Atp12", "Atpg2", "Atpg3", "Atp5f1", "Atph", "Atp5j", "Atp5o", "Atp5i", "Atp5l", "Atpj2", "Tmem70", "C14orf2", "Usmg5")
  
  factores_V_alt=factores_V_alt[factores_V_alt %in% genes]

  factores_alt =c(factores_I_alt, factores_II_alt, factores_III_alt, factores_IV_alt, factores_V_alt)
  
  factores[!factores %in% factores_alt]

  factores_alt[!factores_alt %in% factores]
  length(factores_alt[factores_alt %in% factores])
  
  factores_I=unique(c(factores_I,factores_I_alt))
  factores_II=unique(c(factores_II,factores_II_alt))
  factores_III=unique(c(factores_III,factores_III_alt))
  factores_IV=unique(c(factores_IV,factores_IV_alt))
  factores_V=unique(c(factores_V,factores_V_alt))

  factores =c(factores_I, factores_II, factores_III, factores_IV, factores_V)
  }

pathways_example=list(mito_genes, factores_I, factores_II, factores_III, factores_IV, factores_V)
names(pathways_example)=c("mito_genes", "factores_I", "factores_II", "factores_III", "factores_IV", "factores_V")
save(pathways_example,file="sets_de_factores_ensamblaje.Rdata")

DE_folder="DE_fv/"
## Tu codigo tenia problemas, lo simplifico
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



