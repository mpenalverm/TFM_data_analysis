## SCRIPT 3_Heatmaps: Script for plotting the heatmap of interesting regulons for L.monocytogenes cold growth.
## script developed as a part of final master's dissertation of Marcos Peñalver (2021) for analysing L. monocytogenes regulation under cold growth condition.
## Author: MARCOS PEÑALVER MEDINA
## Contact and code availability: github: mpenalverm

# ********************************************************************************************************************************
# GENERATING HEATMAPS OF INTERESTING REGULONS
# ********************************************************************************************************************************
# Load libraries: 
library(gplots)

#Load the list of manually supervised regulons
heat.table <- read.csv("data/significant_regulons/regulon_manually_curated.csv", sep="\t", header = FALSE)

#Load the differential expression analysis data: 
gene <- read.csv("data/DEGs/DEGs.tsv", sep="\t")
gene <- gene[,c(1,2,4,6,8)]
colnames(gene) <- c("locus_tag", "T_4:0","T_24:0","T_48:0","T_120:0")


#Load Listeriomics annotations:
anno <- read.csv("data/Annotation/Listeria_Annotation_Table_Listeria monocytogenes EGD-e.csv", sep="\t")
anno.short <- anno[,c(1,8,14)]

#Plot heatmap of each regulon:
dir.create("results/regulatory_networks/heatmaps")
for(i in levels(heat.table[,1])){
  regulator.name <- as.vector(heat.table[heat.table$V1 == as.character(i) & heat.table$V3 == "regulator", 2])
  regulated.name <- as.vector(heat.table[heat.table$V1 == as.character(i) & heat.table$V3 == "regulated_gene", 2])
  regulator <- gene[gene$locus_tag == regulator.name,]
  row.names(regulator) <- paste (regulator[,1], as.character(i), sep= " ")
  na.omit(regulator)
  regulator <- as.matrix(regulator[,2:5])
  regulated <- gene[gene$locus_tag %in% regulated.name,]
  label.names <- merge(regulated,anno.short, by.x=1, by.y=1, all.x=T)
  row.names(regulated) <- paste(label.names[,1],label.names[,6], sep= " ")
  na.omit(regulated)
  regulated <- as.matrix(regulated[,2:5])
  total <- rbind(regulator,regulated)
  height <- as.numeric(nrow(total)*300)
  jpeg(paste("results/regulatory_networks/heatmaps/",i, sep=""), height = height, width = 7000, quality = 100, res=300)
  heatmap.2(total,col=redblue, dendrogram='row', Rowv=TRUE, Colv=FALSE, breaks=seq(-4, 4,length.out=100), na.color = "grey", cexRow = 1.5, cexCol = 2, margins = c(12, 15))
  dev.off()
}

## R VERSION INFO:
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS

# gplots_3.1.1
