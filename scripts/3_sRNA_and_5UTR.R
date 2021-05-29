## SCRIPT 3_sRNA_and_5UTR: Script for obtaning signficantly expressed sRNAs of L.monocytogenes and a list of known 5'UTRs for posterior analyses
## script developed as a part of final master's dissertation of Marcos Peñalver (2021) for analysing L. monocytogenes regulation under cold growth condition.
## Author: MARCOS PEÑALVER MEDINA
## Contact and code availability: github: mpenalverm

#Load the libraries:
library(readxl)
library(seqinr)

# ********************************************************************************************************************************
# Obtening the known 5'UTRs of Listeria monocytogenes EGD-e genome
# ********************************************************************************************************************************

#Load the information about known 5'UTRs in Listeria from Weizmann Browser: Wurtzel et al. Molecular Systems Biology, 8:583 (2012).
tss_file <- read_xlsx("data/five_utr/weizmann_broswer_TSS.xlsx")
tss_list <- tss_file[-c(1,2),c(1,4,7,8)]
colnames(tss_list) <- c("gene","strand","tss","start")
tss_list <- subset(tss_list, tss_list$tss != "NA")

#We also load the genome fasta
genome <- read.fasta("data/five_utr/EGD-e_genome.fasta")

#we generate a fasta for each described 5'UTR. If the gene is in the positive strand we take the length of the 5'UTR plus 50nts from the transcription start site. If it is in the negative strand we take the same lenght, but we substract from the TSS coordinate, and consider the complement strand.
fiveutr <- vector()
i <- 1
while (i <= nrow(tss_list)){
  if(tss_list[i,2] == "+"){
    coo.start <- as.numeric(tss_list[i,3])
    coo.end <- as.numeric(tss_list[i,3]) + as.numeric(tss_list[i,4])+50
    fiveutr <- c(fiveutr, paste(genome$NC_003210.1[coo.start:coo.end], collapse = ""))
  }
  else{
    coo.start <- as.numeric(tss_list[i,3]) -as.numeric(tss_list[i,4])-50
    coo.end <- as.numeric(tss_list[i,3])
    utr <- rev(comp(genome$NC_003210.1[coo.start:coo.end]))
    fiveutr <- c(fiveutr, paste(utr, collapse = ""))
  }
  i<- i+1
}

#We transformed it into FASTA FORMAT
five_utr_fasta <- data.frame(id = paste(">", tss_list$gene, sep=""), fasta = fiveutr)
dir.create("results/sRNAs")
write.table(five_utr_fasta, "results/sRNAs/fiveutrs.fna", row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\n")

# ********************************************************************************************************************************
# Obtening a list of sRNAs significantly changed in at least one DE analysis
# ********************************************************************************************************************************

#Load info about sRNA log2fc expression patters:
sRNA <- read.table("path/to/supplementary_material/sRNAS_log2FC.csv", sep="\t", header=TRUE)

#For each Differential Expression analysis (4h vs 0h, 24h vs 0h, 48h vs 0h, 120h vs 0h), we extract the names of those sRNAs that are known to be 
sRNA.4 <- as.vector(sRNA[(sRNA$l2fc_4v0h > 2 | sRNA$l2fc_4v0h < -2) & (sRNA$padj_4v0h < 0.05),1])
sRNA.24 <- as.vector(sRNA[(sRNA$l2fc_24v0h > 2 | sRNA$l2fc_24v0h < -2) & (sRNA$padj_24v0h < 0.05),1])
sRNA.48 <- as.vector(sRNA[(sRNA$l2fc_48v0h > 2 | sRNA$l2fc_48v0h < -2) & (sRNA$padj_48v0h < 0.05),1])
sRNA.120 <- as.vector(sRNA[(sRNA$l2fc_120v0h > 2 | sRNA$l2fc_120v0h < -2) & (sRNA$padj_120v0h < 0.05),1])

#We combine them: 
sRNA.sig <- c(sRNA.4,sRNA.24,sRNA.48,sRNA.120)
sRNA.sig <- unique(sRNA.sig)

#Load coordinates info
sRNA.listeriomics <- read.csv("data/Annotation/Listeria_Srna_Table.txt", sep="\t")
sRNA.list <- sRNA.listeriomics[,c(2,5,3,4)]
colnames(sRNA.list) <- c("gene","strand","start","end")

#As previously, extract info from coordinates
srna.fasta <- vector()
i <- 1
for(i in sRNA.sig){
  if(sRNA.list[sRNA.list$gene == i ,2] == "+"){
    coo.start <- as.numeric(sRNA.list[sRNA.list$gene == i,3])
    coo.end <- as.numeric(sRNA.list[sRNA.list$gene == i,4])
    srna.fasta <- c(srna.fasta, paste(genome$NC_003210.1[coo.start:coo.end], collapse = ""))
  }
  else{
    coo.start <- as.numeric(sRNA.list[sRNA.list$gene == i,3])
    coo.end <- as.numeric(sRNA.list[sRNA.list$gene == i,4])
    seq <- rev(comp(genome$NC_003210.1[coo.start:coo.end]))
    srna.fasta <- c(srna.fasta, paste(seq, collapse = ""))
  }
}

#Give a fasta format
srna.fasta.real <- data.frame(id = paste(">", sRNA.sig, sep=""), fasta = srna.fasta)
write.table(srna.fasta.real, "results/sRNAs/significant_sRNA.fna", row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\n")

## R VERSION INFO:
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS

# seqinr_4.2-5 readxl_1.3.1
