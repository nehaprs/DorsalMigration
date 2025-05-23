library(biomaRt)
library(dplyr)
#id_vector <- read_excel("enrichment/id_vector20n11.xlsx", col_names = FALSE)

versions = listEnsemblArchives()
versions = versions$version
versions = versions[-1]

id_vector = in11not20$gene

# e.g. use Ensembl v114
ensembl114 <- useEnsembl(
  biomart    = "genes", 
  dataset    = "xtropicalis_gene_ensembl", 
  version    = 114
)

# query
mapping <- getBM(
  attributes    = c("ensembl_gene_id", "external_gene_name"),
  filters       = "ensembl_gene_id",
  values        = id_vector,
  mart          = ensembl114
)













#find version with the least empty vectors
sum(mapping$external_gene_name == "")




mapped = full_join(in11not20, mapping, by = c("gene"="ensembl_gene_id" ))
write_xlsx(mapped,"enrichment/mapped11not20.xlsx")
