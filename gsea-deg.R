library(dplyr)
library(tibble)
library(clusterProfiler)
library(msigdbr)
library(gprofiler2)
library(enrichplot)
library(ggplot2)

setwd("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/integrated/all4/high_res/diffexp/filt-r1")
library(readxl)
df <- read_excel("filt_r1_s19d_v_non-nc.xlsx")
df2 = df[,c(2,5,6)]

## Map Xenopus genes to human orthologs (gprofiler)

orth <- gorth(
  query = df2$gene,
  source_organism = "xtropicalis",
  target_organism = "hsapiens"
)

orth2 <- orth %>%
  select(input, ortholog_name) %>%
  
  rename(
    gene = input,
    human_gene = ortholog_name
  )

df_human <- df2 %>%
  inner_join(orth2, by = "gene") %>%
  filter(!is.na(human_gene))


df_ranked = df_human %>% 
  mutate(rank_score = sign(avg_log2FC)* -log10(p_val_adj + 1e-300)) %>% # a common ranking formula
  arrange(desc(rank_score))



gene_list = df_ranked$rank_score
names(gene_list) <- df_ranked$human_gene
names(gene_list) = df_ranked$human_gene
gene_list = sort(gene_list, decreasing = TRUE)


msig_reactome <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CP:REACTOME"
) %>%
  select(gs_name, gene_symbol)


gsea_reactome<- GSEA(
  geneList = gene_list,
  TERM2GENE = msig_reactome,
  pvalueCutoff = 0.05,
  verbose = FALSE
)



gsea_reactome@result$ID = sub("^REACTOME_", "", gsea_reactome@result$ID)
gsea_reactome@result$Description = sub("^REACTOME_", "", gsea_reactome@result$Description)
rownames(gsea_reactome@result) = sub("^REACTOME_", "", rownames(gsea_reactome@result))
writexl::write_xlsx( as.data.frame(gsea_reactome),"filt_r1_s19d_v_non-nc.xlsx-gsea_reactome.xlsx" )


p = dotplot(gsea_reactome, showCategory = 20) +
  ggtitle("Reactome GSEA")

ggsave("filt_r1_s19d_v_non-nc.xlsx-gsea_reactome.png", plot = p, width = 10, height = 12)

