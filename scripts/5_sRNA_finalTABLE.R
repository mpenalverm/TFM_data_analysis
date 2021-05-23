## SCRIPT 5_sRNA_finalTABLE: Script for evaluating the target candidates of sRNAs regulated under cold growth in L.monocytogenes
## script developed as a part of final master's dissertation of Marcos Peñalver (2021) for analysing L. monocytogenes regulation under cold growth condition.
## Author: MARCOS PEÑALVER MEDINA
## Contact and code availability: github: mpenalverm

#Load the libraries:

# ********************************************************************************************************************************
# Merge all sRNA info into a final table:
# ********************************************************************************************************************************
#Load data about significant targets prediction (p-value <= 0.005) according to TargetRNA2 and IntaRNA prediction tools
data.srna <- read.csv("data/sRNAs/Inta_targetprediction_pValue_005.csv", sep="\t")

#Load information log2FC differential expression analysis in proteomics and transcriptomics: 
#Proteomics:
protein <- read.csv("data/proteomics/DEPs.csv", sep="\t")
protein <- protein[,c(1,2,4,6,8)] #Keep only the columns with log2INFO
colnames(protein) <- c("locus_tag", "P_4:0","P_24:0","P_48:0","P_120:0")
#Transcriptomis:
gene <- read.csv("data/DEGs/DEGs.tsv", sep="\t")
gene <- gene[,c(1,2,4,6,8)]  #Keep only the columns with log2INFO
colnames(gene) <- c("locus_tag", "T_4:0","T_24:0","T_48:0","T_120:0")

#Merge the tables:
de.profile <- merge(data.srna, gene, by.x=2, by.y=1, all.x = TRUE)
de.profile <- merge(de.profile, protein, by.x=1, by.y=1, all.x = TRUE)
de.profile$corr <- "NA"
de.profile$corr_protein <- "NA"
de.profile <- de.profile[complete.cases(de.profile[ , 4:7]),]

#Calculate correlation between the sRNA and its putative target: 
i <- 1
while(i <= nrow(de.profile)){
  aa <- as.numeric(de.profile[i,4:7])
  bb <- as.numeric(de.profile[as.character(de.profile$target) == as.character(de.profile[i,"sRNA"]),4:7])
  de.profile[i,"corr"] <- cor(aa, bb, method = "pearson")
  i <- i+1
}

#Calculate correlation between the transcript and protein levels if possible 
i <- 1
while(i <= nrow(de.profile)){
  aa <- as.numeric(de.profile[i,4:7])
  bb <- as.numeric(de.profile[i,8:11])
  if(!is.na(aa) & !is.na(bb))
  de.profile[i,"corr_protein"] <- cor(aa, bb, method = "pearson")
  i <- i+1
}

write.table(de.profile, "results/sRNAs/sRNA_target_finaltable.csv", row.names = FALSE, quote = FALSE, sep="\t")


## R VERSION INFO:
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS
