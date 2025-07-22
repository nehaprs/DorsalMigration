library(clusterProfiler)
library(biomaRt)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(enrichplot)

# ────────────────────────────────────────────────────────────────
# STEP 1: Prepare biomart mapping
xtrop <- useMart("ensembl", dataset = "xtropicalis_gene_ensembl")

makers <- read_excel("~/BINF/scrnaseq general/dorsal migration/full head/version4/CCScoring/filt_markers/filtmarkers_resolution_2.1.xlsx")
#makers = filtfiltmarkers_resolution_0_5
makers_all = makers[makers$cluster == 22,]
makers = makers[makers$cluster == 22,7]
markers_xt = makers





markers = makers_all
#which(makers_all == "wnt1")
gene_list <- markers %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::pull(avg_log2FC)
rownames(markers) = markers$gene
names(gene_list) <- rownames(markers)[order(markers$avg_log2FC, decreasing = TRUE)]
head(gene_list)
# Convert Ensembl IDs → Entrez IDs + symbols
conversion <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = names(gene_list_named_entrez),
  mart = xtrop
) %>%
  filter(!is.na(entrezgene_id)) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Rebuild gene list with Entrez IDs as names
gene_list_mapped <- gene_list_named_entrez[names(gene_list_named_entrez) %in% conversion$ensembl_gene_id]
names(gene_list_mapped) <- conversion$entrezgene_id[match(names(gene_list_mapped), conversion$ensembl_gene_id)]
head(gene_list_mapped)
gene_list_entrez <- sort(gene_list_mapped, decreasing = TRUE)

# ────────────────────────────────────────────────────────────────
#  GSEA using KEGG
kegg_gsea <- gseKEGG(
  geneList     = gene_list_entrez,
  organism     = "xtr",
  nPerm        = 1000,
  minGSSize    = 5,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = TRUE
)

#  Extract genes per pathway (with gene symbols)
kegg_df <- as.data.frame(kegg_gsea) %>%
  separate_rows(core_enrichment, sep = "/")
# Add gene symbols
conversion$entrezgene_id = as.character(conversion$entrezgene_id)

kegg_df <- kegg_df %>%
  left_join(conversion, by = c("core_enrichment" = "entrezgene_id")) %>%
  group_by(ID, Description, NES, pvalue, p.adjust, qvalue)


kegg_df = kegg_df[,-c(12,13,14,15)]

kegg_df = kegg_df %>%
  summarise(
    EntrezIDs = paste(unique(core_enrichment), collapse = ","),
    GeneSymbols = paste(unique(external_gene_name), collapse = ","),
    .groups = "drop"
  )
#NES = Normalized Enrichment Score

write_xlsx(kegg_df,"gsea_kegg_res2.4cl35.xlsx")
# View top results
head(kegg_df)

# ────────────────────────────────────────────────────────────────
#GSEA like GO with gprofiler, since org.Xt... not available in bioconducter
#___________________________________________________________________
install.packages("gprofiler2")

library(gprofiler2)

# Gene names (symbols or Ensembl IDs)
genes_ranked <- names(gene_list_named_entrez)  

gost_res <- gost(
  query = genes_ranked,
  organism = "xtropicalis",
  ordered_query = TRUE,
  significant = TRUE,
  sources = c("GO:BP", "GO:MF", "GO:CC"),
  evcodes = TRUE 
)

# Extract result
go_results <- gost_res$result
colnames(go_results)


# Extract and flatten Ensembl IDs from all intersections
ensembl_ids <- go_results$intersection %>%
  str_split(",") %>%
  unlist() %>%
  unique()

xtrop <- useMart("ensembl", dataset = "xtropicalis_gene_ensembl")

ensembl_to_symbol <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters    = "ensembl_gene_id",
  values     = ensembl_ids,
  mart       = xtrop
)

# Create named vector for fast replacement
id2symbol <- setNames(ensembl_to_symbol$external_gene_name, ensembl_to_symbol$ensembl_gene_id)

# Replace Ensembl IDs with symbols in each row
go_results <- go_results %>%
  mutate(
    GeneSymbols = str_split(intersection, ",") %>%
      lapply(function(x) id2symbol[x]) %>%  # map each Ensembl ID to symbol
      lapply(function(x) paste(na.omit(x), collapse = ", ")) %>%
      unlist()
  )

write_xlsx(go_results,"gsea_gores2cl25.xlsx")
