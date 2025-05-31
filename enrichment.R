#======================================================
#enrichment analysis of the three clusters of interest
#======================================================

# This script takes a vector of Xenopus tropicalis marker gene \,
# maps them to both Xenopus Entrez IDs and human ortholog Entrez IDs via biomaRt,
# and then performs:
#   1) GO Biological Process enrichment using human orthologs (clusterProfiler + org.Hs.eg.db)
#   2) KEGG pathway enrichment directly for Xenopus (organism = "xtr")


library(biomaRt)          
library(clusterProfiler)  
library(org.Hs.eg.db)     
library(enrichplot)       
library(ggplot2)

#marker genes:
#1. From paired analysis gene significantly higher in cluster 20/12 than in cluster 11

markers_xt = row.names(pairedmarkers20vs11)
#set up biomart connection
