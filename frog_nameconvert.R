frog_nameconvert = function(s.object){
  library(biomaRt)
  
  # Extract and clean gene IDs
  original_ids = rownames(s.object)
  cleaned_ids = gsub("^Xetrov", "", original_ids)
  cleaned_ids = gsub("m$", "", cleaned_ids)
  
  message("Connecting to Ensembl...")
  mart <- tryCatch({
    useEnsembl("ensembl", dataset = "xtropicalis_gene_ensembl", mirror = "useast")
  }, error = function(e) {
    stop("Failed to connect to Ensembl mirror. Try again later or with a different mirror.")
  })
  
  # Step 3: Query Ensembl for Ensembl gene IDs
  message("Querying Ensembl gene IDs...")
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id"),
    filters = "entrezgene_id",
    values = cleaned_ids,
    mart = mart
  )
  
  # Step 4: Build mapping vector
  id_map <- setNames(gene_map$ensembl_gene_id, gene_map$entrezgene_id)
  
  # Step 5: Map and clean names
  new_rownames <- id_map[cleaned_ids]
  new_rownames_clean <- ifelse(is.na(new_rownames) | new_rownames == "",
                               original_ids,  # fallback to original
                               new_rownames)
  new_rownames_clean <- make.unique(new_rownames_clean)
  print(head(new_rownames_clean))
    
    # Assign names
    rownames(s.object) = new_rownames_clean
    
    return(s.object)
    
    
}


