## SCRIPT 0_GO_ANNOTATION: Script for generating a Annotation package for Listeria monocytogenes EGD-e containing Biological Process GO terms
## script developed as a part of final master's dissertation of Marcos Peñalver (2021) for analysing L. monocytogenes regulation under cold growth condition.
## Author: MARCOS PEÑALVER MEDINA
## Contact and code availability: github: mpenalverm

# ********************************************************************************************************************************
# GENERATING A LIBRARY WITH GO BP ANNOTATIONS FOR LISTERIA
# ********************************************************************************************************************************
# Load libraries:
library(AnnotationForge)

dir.create("packages/lismoANNOFORGE") #Create needed directories

#Load GO annotations and generate the package:
fGO <- read.table("data/Annotation/GO_terms/lmo2go.gomap", header=FALSE)
fGO_2 <- unique(fGO)
fGO_2 <- cbind(fGO_2, "unknown")
colnames(fGO_2) <- c("GID","GO","EVIDENCE")
fGO_2 <- subset(fGO_2, fGO_2$GO != "")
makeOrgPackage( go=fGO_2,
                version="0.1",
                maintainer="Some One <so@someplace.org>",
                author="Some One <so@someplace.org>",
                outputDir = "packages/lismoANNOFORGE/",
                tax_id="169963",
                genus="Listeria",
                species="monocytogenes",
                goTable="go")


#Install the package in R
install.packages("packages/lismoANNOFORGE/org.Lmonocytogenes.eg.db", repos=NULL)

## R VERSION INFO:
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS

# AnnotationForge_1.28.0