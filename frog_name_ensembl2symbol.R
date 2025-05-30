frog_name_ensembl2symbol = function(s.object){
  library(biomaRt)
  #convert ensembl ID to gene symbol
  # Extract and clean gene IDs
  original_ids = rownames(s.object)
  
  
  message("Connecting to Ensembl...")
  mart <- tryCatch({
    useEnsembl("ensembl", dataset = "xtropicalis_gene_ensembl", mirror = "useast")
  }, error = function(e) {
    stop("Failed to connect to Ensembl mirror. Try again later or with a different mirror.")
  })
  
  #  Query Ensembl for Ensembl gene IDs
  message("Querying Ensembl gene IDs...")
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = original_ids,
    mart = mart
  )
  
  # 4. Build a named vector: names = Ensembl ID → values = symbol
  id2symbol <- setNames(
    gene_map$external_gene_name,
    gene_map$ensembl_gene_id
  )
  
  # 5. Map and fall back:
  new_names <- ifelse(
    !is.na(id2symbol[original_ids]) & id2symbol[original_ids] != "",
    id2symbol[original_ids],  # use mapped symbol
    original_ids              # fallback
  )
  
  # 6. Ensure uniqueness
  new_names <- make.unique(new_names)
  
  # 7. Assign back and return
  rownames(s.object) <- new_names
  return(s.object)
  
  
}


