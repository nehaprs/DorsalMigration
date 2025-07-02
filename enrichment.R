#======================================================
#enrichment analysis of the three clusters of interest
#part 1: KEGG
#part2: GSEA
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
library(dplyr)
#marker genes:
#1. From paired analysis gene significantly higher in cluster 20/12 than in cluster 11
#enrichment of cluster 7 markers
makers <- read_excel("filtmarkers_resolution_2.4.xlsx")
#makers = filtfiltmarkers_resolution_0_5
makers_all = makers[makers$cluster == 35,]
makers = makers[makers$cluster == 35,7]
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
    qvalueCutoff = 0.05
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
    title        = "KEGG Dotplot (X. tropicalis –res 2.4 cluster 35)"
  )
}

#add gene id column
kegg_df <- as.data.frame(ekegg_xtr)

kegg_genes_long <- kegg_df %>%
  separate_rows(geneID, sep = "/") 

#  Map Entrez IDs to gene symbols
xtrop <- useMart("ensembl", dataset = "xtropicalis_gene_ensembl")

entrez2symbol <- getBM(
  attributes = c("entrezgene_id", "external_gene_name"),
  filters    = "entrezgene_id",
  values     = unique(kegg_genes_long$geneID),
  mart       = xtrop
)

entrez2symbol$entrezgene_id = as.character(entrez2symbol$entrezgene_id)
#Merge gene symbols into KEGG result
kegg_genes_with_symbols <- kegg_genes_long %>% 
  left_join(entrez2symbol, by = c("geneID" = "entrezgene_id")) %>% 
  as.data.frame() %>%
  dplyr::select(ID, Description, geneID, external_gene_name, pvalue, p.adjust, qvalue)


kegg_collapsed <- kegg_genes_with_symbols %>%
  group_by(ID, Description, pvalue, p.adjust, qvalue) %>%
  summarise(
    EntrezIDs = paste(unique(geneID), collapse = ","),
    GeneSymbols = paste(unique(external_gene_name), collapse = ","),
    .groups = "drop"
  )

write_xlsx(kegg_collapsed,"kegg_res2.4cl35.xlsx")



#save output files
'
if (exists("ekegg_xtr") && nrow(as.data.frame(ekegg_xtr)) > 0) {
  kegg_xtr_df <- as.data.frame(ekegg_xtr)
  write.csv(
    kegg_xtr_df,
    file = "~/BINF/scrnaseq general/dorsal migration/full head/version2/kegg_res2.4cl35.csv",
    row.names = FALSE
  )
}
'


############################################
#GSEA
#############################################
library(clusterProfiler)
library(org.Hs.eg.db) 
markers = makers_all

gene_list <- markers %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::pull(avg_log2FC)
rownames(markers) = markers$gene
names(gene_list) <- rownames(markers)[order(markers$avg_log2FC, decreasing = TRUE)]
head(gene_list)
xtrop <- useMart("ensembl", dataset = "xtropicalis_gene_ensembl")
gene_names = names(gene_list)

# Separate into Ensembl-style and non-Ensembl (likely symbols)
ensembl_like <- grep("^ENSXETG", gene_names, value = TRUE)
symbol_like  <- setdiff(gene_names, ensembl_like)

# Get mapping for Ensembl IDs → themselves (plus symbols if needed)
ens_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = ensembl_like,
  mart = xtrop
)
'
# Get mapping for gene symbols → Ensembl IDs
sym_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = symbol_like,
  mart = xtrop
)

# Combine both mappings
conversion_table <- bind_rows(
  ens_map %>% rename(input = ensembl_gene_id),
  sym_map %>% rename(input = external_gene_name)
)
'

tail(gene_names)

conversion = ens_map
# Filter out NAs and duplicates
#conversion <- conversion[!is.na(conversion$entrezgene_id), ] %>%
 # distinct(external_gene_name, .keep_all = TRUE)

# Rebuild gene_list with Entrez IDs as names
gene_list_named_entrez <- gene_list[names(gene_list) %in% conversion$ensembl_gene_id]


kegg_gsea <- gseKEGG(
  geneList     = sort(gene_list_named_entrez, decreasing = TRUE),
  organism     = "xtr",
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = TRUE
)





###############

xtrop <- useMart("ensembl", dataset = "xtropicalis_gene_ensembl")

conversion <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters    = "ensembl_gene_id",
  values     = names(gene_list_named_entrez),
  mart       = xtrop
)

# Filter out NAs
conversion <- conversion %>%
  filter(!is.na(entrezgene_id)) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)
# Merge to get Ensembl IDs in order
gene_list <- gene_list[names(gene_list) %in% conversion$external_gene_name]
names(gene_list) <- converted$ensembl_gene_id[match(names(gene_list), converted$external_gene_name)]
head(gene_list)
'


# Check if Xenopus is supported
organism <- "xtr"  # KEGG code for Xenopus tropicalis

# Run GSEA using KEGG
kegg_gsea <- gseKEGG(
  geneList     = gene_list,
  organism     = "xtr",
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = TRUE
)