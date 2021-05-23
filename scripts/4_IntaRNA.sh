## SCRIPT 4_intaRNA: Script for predicting sRNA targets in the 5'UTRs of L. monocytogenes genome using IntaRNA (DOI: DOI:10.1093/nar/gkx279 and DOI:10.1093/bioinformatics/btn544)
## script developed as a part of final master's dissertation of Marcos Peñalver (2021) for analysing L. monocytogenes regulation under cold growth condition.
## Author: MARCOS PEÑALVER MEDINA
## Contact and code availability: github: mpenalverm

IntaRNA -t results/sRNAs/fiveutrs.fna -q results/sRNAs/significant_sRNA.fna --out results/sRNAs/IntaRNA_targets_results --outMode=C --outCsvCols=id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E,E_norm

#Calculation of significance: 
Rscript --vanilla scripts/4_1_IntaRNA_CSV_p-value.R results/sRNAs/IntaRNA_targets_results results/sRNAs/IntaRNA_targets_results_withpvalue.csv E_norm

