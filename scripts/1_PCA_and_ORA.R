## SCRIPT1_PCA_and_ORA: Script for functional analysis of L.monocytogenes cold Growth data: Principal Component Analysis and OverRepresentation Analysis. 
## script developed as a part of final master's dissertation of Marcos Peñalver (2021) for analysing L. monocytogenes regulation under cold growth condition.
## Author: MARCOS PEÑALVER MEDINA
## Contact and code availability: github: mpenalverm

# ********************************************************************************************************************************
# Setup
# ********************************************************************************************************************************

dir.create("results") #create a results directory.

#Load necessary libraries: 
library("DESeq2")
library(pheatmap)
library("FactoMineR")
library("factoextra")
library(corrplot)
library(clusterProfiler)
library(org.Lmonocytogenes.eg.db)
library(dplyr)
library(ggplot2)

#Load counts data:
deseq.data <- readRDS("rna-seq-star-deseq/deseq2/all.rds")
counts <- vst(deseq.data, blind=FALSE) # Normalization using VST (varianceStabilizingTransformation)
vst.values <- counts@assays@data@listData[[1]] # Access to vst.values

#Samples correlation
# load libraries pheatmap to create the heatmap plot
# calculate between-sample distance matrix
sampleDistMatrix <- as.matrix(dist(t(assay(counts))))

# create figure in PNG format
png("results/pca/sample_distance_matrix.png")
pheatmap(sampleDistMatrix)
dev.off() 

# ********************************************************************************************************************************
# PCA using FactoMineR package:
# ********************************************************************************************************************************

#Data preparation
vst.values.t <- t(as.matrix(vst.values)) # Obtain the transposed matrix

# Principal Component Analysis (PCA)
vst.values.pca <- PCA(vst.values.t, scale.unit = FALSE, graph = FALSE) # Calculate PCA

dir.create("results/pca")
png("results/pca/PCA_star_FactoMineR.png") #Get the PCA plot
fviz_pca_ind(vst.values.pca, repel=TRUE)
dev.off()

#Calculate those transcripts correlated with each principal component:
vst.values.desc <- dimdesc(vst.values.pca, axes = c(1,2), proba = 2) # Calculate descriptive paramethers of each principal component:

#For PC1
PC1 <-  as.data.frame(vst.values.desc$Dim.1$quanti)
pc1.pos = row.names(PC1[PC1$correlation >= 0 & PC1$p.value <= 0.05,])
pc1.neg = row.names(PC1[PC1$correlation <= 0 & PC1$p.value <= 0.05,])
write.table(pc1.pos,"results/pca/PC1_positive_correlated_transcripts.csv", sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(pc1.neg,"results/pca/PC1_negative_correlated_transcripts.csv", sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#For PC2
PC2 <- as.data.frame(vst.values.desc$Dim.2$quanti)
pc2.pos = row.names(PC2[PC2$correlation >= 0 & PC2$p.value <= 0.05,])
pc2.neg = row.names(PC2[PC2$correlation <= 0 & PC2$p.value <= 0.05,])
write.table(pc2.pos,"results/pca/PC2_positive_correlated_transcripts.csv", sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(pc2.neg,"results/pca/PC2_negative_correlated_transcripts.csv", sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


#Get a list of significantly correlated transcripts:
gene.pca <- list("pc1.t.pos" = pc1.pos,"pc1.t.neg" = pc1.neg,
                 "pc2.t.pos" = pc2.pos,"pc2.t.neg" = pc2.neg)



# ********************************************************************************************************************************
# Functional Over-representation analysis (ORA) of PCA correlated clusters using ClusterProfiler package and KEGG and GO annotations:
# ********************************************************************************************************************************

# KEGG
ck.kegg <- compareCluster(geneCluster = gene.pca, fun = "enrichKEGG",  organism="lmo") 
#Get the plot image:
dir.create("results/functional_analysis")
png("results/functional_analysis/KEGG_ORA.png") 
dotplot(ck.kegg,showCategory=NULL)
dev.off()

#GO BP
ck.go <- compareCluster(geneCluster = gene.pca, fun = "enrichGO", OrgDb="org.Lmonocytogenes.eg.db", keyType = 'GID',  ont = "BP")
ck.go.short <- clusterProfiler::simplify(ck.go, cutoff=0.5, by="p.adjust", select_fun=min)
#Get the plot image:
png("results/functional_analysis/GO_BP_ORA.png") 
dotplot(ck.go.short,showCategory=NULL) 
dev.off()

# ********************************************************************************************************************************
# REGULON Over-representation analysis (ORA) of PCA correlated clusters using ClusterProfiler package and KEGG and GO annotations:
# ********************************************************************************************************************************

#Load regulons database and transform into a compatible ClusterProfiler Format:
gmt <- read.gmt("data/regulons_db/regulon.gmt") #Load the database
wpid2gene <- gmt %>% dplyr::select(ont, gene) # Transform into a compatible format

#Perform ORA wih regulons gene-sets for each PC and save summary info: 
ck.regulon.pc1pos <- enricher(gene.pca[[1]], TERM2GENE = wpid2gene)
ck.regulon.pc1pos.summ <- data.frame("NAME" = ck.regulon.pc1pos@result$ID, "SIZE" = ck.regulon.pc1pos@result$Count, "p.value" = ck.regulon.pc1pos@result$p.adjust, "genes" = ck.regulon.pc1pos@result$geneID, "CORRELATION" ="PC1 pos")
ck.regulon.pc1pos.summ <- subset(ck.regulon.pc1pos.summ, ck.regulon.pc1pos.summ$p.value <= 0.05)
ck.regulon.pc1neg <- enricher(gene.pca[[2]], TERM2GENE = wpid2gene)
ck.regulon.pc1neg.summ <- data.frame("NAME" = ck.regulon.pc1neg@result$ID, "SIZE" = ck.regulon.pc1neg@result$Count, "p.value" = ck.regulon.pc1neg@result$p.adjust, "genes" = ck.regulon.pc1neg@result$geneID, "CORRELATION" ="PC1 neg")
ck.regulon.pc1neg.summ <- subset(ck.regulon.pc1neg.summ, ck.regulon.pc1neg.summ$p.value <= 0.05)
ck.regulon.pc2pos <- enricher(gene.pca[[3]], TERM2GENE = wpid2gene)
ck.regulon.pc2pos.summ <- data.frame("NAME" = ck.regulon.pc2pos@result$ID, "SIZE" = ck.regulon.pc2pos@result$Count, "p.value" = ck.regulon.pc2pos@result$p.adjust, "genes" = ck.regulon.pc2pos@result$geneID, "CORRELATION" ="PC2 pos")
ck.regulon.pc2pos.summ <- subset(ck.regulon.pc2pos.summ, ck.regulon.pc2pos.summ$p.value <= 0.05)
ck.regulon.pc2neg <- enricher(gene.pca[[4]], TERM2GENE = wpid2gene)
ck.regulon.pc2neg.summ <- data.frame("NAME" = ck.regulon.pc2neg@result$ID, "SIZE" = ck.regulon.pc2neg@result$Count, "p.value" = ck.regulon.pc2neg@result$p.adjust, "genes" = ck.regulon.pc2neg@result$geneID, "CORRELATION" ="PC2 neg")
ck.regulon.pc2neg.summ <- subset(ck.regulon.pc2neg.summ, ck.regulon.pc2neg.summ$p.value <= 0.05)
ck.regulon <- rbind(ck.regulon.pc1pos.summ,ck.regulon.pc1neg.summ,ck.regulon.pc2pos.summ,ck.regulon.pc2neg.summ)

#Result manually curation (ONLY COMPATIBLE WITH THIS SPECIFIC DATASET)
ck.regulon[,1] <- as.vector(ck.regulon[,1])
ck.regulon[13,1] <- "ccpA" #change ccpA name
ck.regulon <- ck.regulon[-c(1,2,5,7,8,14,15),] #Delete regulons repeated for more than one organism.

#Display graphic summary
S1<- ggplot(ck.regulon[,-4], aes(x= CORRELATION, y=NAME, size=SIZE, color=p.value, group=CORRELATION)) + geom_point(alpha = 0.8) + 
  theme_classic()
S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(0.000000000000000000, 0.05))
S1+scale_size(range = c(2, 8))

dir.create("results/regulatory_networks")
png("results/regulatory_networks/regulons.png")
S1
dev.off()

# ********************************************************************************************************************************
# CYTOSCAPE METADATA PREPARATION
# ********************************************************************************************************************************
# REGULON METADATA
#Get a dataframe of genes belonging to a each significantly enriched regulon in each PC: 
cytoscape.metadata.reg <- data.frame("regulator" = NA,"regulated_gene" = NA)
cytoscape.metadata.reg <- cytoscape.metadata.reg[-1,]
i <-1
while(i <= nrow(ck.regulon)){
  a <- strsplit(as.character(ck.regulon[i,4]),"/")
  a <- a[[1]]
  y <- 1
  aa <- data.frame("regulator" = NA, "regulated_gene" = NA)
  aa <- aa[-1,]
  while(y <=length(a)){
    bb <- data.frame("regulator" = ck.regulon[i,1], "regulated_gene" = a[y])
    aa <- rbind(aa,bb)
    y <- y + 1
  }
  cytoscape.metadata.reg <- rbind(cytoscape.metadata.reg, aa)
  i <- i + 1
}

#Collapse information of genes belonging to more than one regulon:
cytoscape.metadata.reg.final <- cytoscape.metadata.reg %>%
  group_by(regulated_gene) %>%
  summarise_each(funs(paste(., collapse = "; ")))

dir.create("results/cytoscape_metadata")
write.table(cytoscape.metadata.reg.final, file ="results/cytoscape_metadata/lmo2regulon.csv", row.names = FALSE, sep ="\t", quote=FALSE)

# KEGG METADATA
#Get a dataframe of genes belonging to a each significantly enriched KEGG pathway in each PC:
ck.kegg.summ <- data.frame("NAME" <- ck.kegg@compareClusterResult$Description, "SIZE" <- ck.kegg@compareClusterResult$Count, "p.value" <- ck.kegg@compareClusterResult$qvalue, "genes" <-ck.kegg@compareClusterResult$geneID)
ck.kegg.summ <- subset(ck.kegg.summ, ck.kegg.summ$X.p.value.....ck.kegg.compareClusterResult.qvalue <= 0.05)
cytoscape.metadata.kegg <- data.frame("regulator" = NA,"regulated_gene" = NA)
cytoscape.metadata.kegg <- cytoscape.metadata.kegg[-1,]
i <-1
while(i <= nrow(ck.kegg.summ)){
  a <- strsplit(as.character(ck.kegg.summ[i,4]),"/")
  a <- a[[1]]
  y <- 1
  aa <- data.frame("regulator" = NA, "regulated_gene" = NA)
  aa <- aa[-1,]
  while(y <=length(a)){
    bb <- data.frame("regulator" = ck.kegg.summ[i,1], "regulated_gene" = a[y])
    aa <- rbind(aa,bb)
    y <- y + 1
  }
  cytoscape.metadata.kegg <- rbind(cytoscape.metadata.kegg, aa)
  i <- i + 1
}

#Collapse information of genes belonging to more than one kegg:
cytoscape.metadata.kegg.final <- cytoscape.metadata.kegg %>%
  group_by(regulated_gene) %>%
  summarise_each(funs(paste(., collapse = "; ")))

write.table(cytoscape.metadata.kegg.final, file ="results/cytoscape_metadata/lmo2kegg.csv", row.names = FALSE, sep ="\t", quote=FALSE)

## R VERSION INFO:
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS

# dplyr_1.0.3                  org.Lmonocytogenes.eg.db_0.1 AnnotationDbi_1.48.0         clusterProfiler_3.14.3       corrplot_0.84               
# factoextra_1.0.7             ggplot2_3.3.3                FactoMineR_2.4               pheatmap_1.0.12              DESeq2_1.26.0               
# SummarizedExperiment_1.16.1  DelayedArray_0.12.3          BiocParallel_1.20.1          matrixStats_0.55.0           Biobase_2.46.0              
# GenomicRanges_1.38.0         GenomeInfoDb_1.22.0          IRanges_2.20.2               S4Vectors_0.24.4             BiocGenerics_0.32.0 
