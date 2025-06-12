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
#enrichment of cluster 7 markers
makers <- read_excel("filtmarkers_resolution_0.7.xlsx")
makers = makers[makers$cluster == 7, 7]
markers_xt = makers
#set up biomart connection



mart_xt = useEnsembl(biomart = "genes", dataset = "xtropicalis_gene_ensembl")


## Map Xenopus tropicalis SYMBOLS → Xenopus Entrez IDs 



bm_xt2hs = getBM(attributes = c("external_gene_name", "entrezgene_id"),
                 filters = "external_gene_name",
                 values = markers_xt,
                 mart = mart_xt
)
                   
# Extract unique, non-NA Xenopus Entrez IDs
xtr_entrez <- unique(na.omit(bm_xt2hs$entrezgene_id))
if (length(xtr_entrez) == 0) {
  warning("None of the Xenopus SYMBOLS mapped to an Xt Entrez ID.")
}   


#============================================
#KEGG pathway enrichment
#============================================

if (length(xtr_entrez) > 0) {
  ekegg_xtr <- enrichKEGG(
    gene         = xtr_entrez,
    organism     = "xtr",           # KEGG code for Xenopus tropicalis
    keyType      = "ncbi-geneid",   # our vector is Entrez IDs
    pvalueCutoff = 0.05,
    pAdjustMethod= "BH",
    qvalueCutoff = 0.10
  )}


if (nrow(as.data.frame(ekegg_xtr)) == 0) {
  message("No KEGG pathways passed the p-value/q-value thresholds for Xenopus tropicalis.")
} else {
  message("Top 10 KEGG pathways for X. tropicalis cluster-7 markers:")
  print(head(ekegg_xtr, n = 10))
  
  # visualize top 10 KEGG pathways
  barplot(
    ekegg_xtr,
    showCategory = 10,
    title        = "KEGG Enrichment (X. tropicalis –cluster7 markers )"
  )
  
  dotplot(
    ekegg_xtr,
    showCategory = 10,
    title        = "KEGG Dotplot (X. tropicalis –paired markers 12vs11)"
  )
}


#save output files

if (exists("ekegg_xtr") && nrow(as.data.frame(ekegg_xtr)) > 0) {
  kegg_xtr_df <- as.data.frame(ekegg_xtr)
  write.csv(
    kegg_xtr_df,
    file = "markerssubcluster7_KEGG_enrichment.csv",
    row.names = FALSE
  )
}
