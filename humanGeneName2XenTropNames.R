#converts human gene names to gene names of X. tropicalis

humanGeneName2XenTropNames <- function(x) {
  # Load human Ensembl dataset
  human <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
  
  # Load Xenopus tropicalis Ensembl dataset
  xenopus <- biomaRt::useEnsembl("ensembl", dataset = "xtropicalis_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
  
  # Link both datasets and retrieve Xenopus genes from human genes
  genes.list <- biomaRt::getLDS(
    attributes = c("hgnc_symbol"),
    filters = "hgnc_symbol",
    values = x,
    mart = human,
    attributesL = c("external_gene_name"),
    martL = xenopus,
    uniqueRows = FALSE
  )
  
  # Extract unique Xenopus gene names
  xt.gene.list <- unique(genes.list[, 2])
  return(xt.gene.list)
}
