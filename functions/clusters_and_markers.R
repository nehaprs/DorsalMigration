#================================
'
function does the following:
1. finds the clusters enriched in listed nc markers and zic genes
2. for clusters enriching both foxd3 and zic3, find the wnt and bmp lignds they enrich
'
#=================================



clusters_and_markers = function(folder = "filt_markers"){
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(stringr)
  
  #define nc and migrtory genes
  genes_nc <- c("foxd3", "sox10", "zic3", "zic2", "snai1", "twist1", "sox9", "tfap2a")
  
  #list all excel files in the folder
  excel_files = list.files(folder, pattern = "\\.xlsx$", full.names = TRUE)
  
  #create outbut workbooks
  nc_Candidates_wb = createWorkbook()
  wnt_ligands_wb = createWorkbook()
  bmp_ligands_wb = createWorkbook()
  
  for (file_path in excel_files) {
    sheet_name = sub(".*filt_markers_resolution_", "", tools::file_path_sans_ext(basename(file_path)))
    
    #read file
    df = read_excel(file_path)
    # Ensure relevant columns exist
    if (!all(c("cluster", "gene") %in% colnames(df))) {
      warning(paste("Skipping file:", file_path, "due to missing columns."))
      next
    }
    
    # 1: Find clusters expressing genes_nc
    df_nc = df %>% filter(tolower(gene) %in% genes_nc)
    nc_candidates = unique(df_nc$cluster)
    nc_expr = df_nc %>% select(cluster, gene)
    
    # Add to nc_candidates workbook
    addWorksheet(nc_Candidates_wb, sheet_name)
    writeData(nc_Candidates_wb,sheet_name, nc_expr)
    
    
    # 2: Find clusters expressing BOTH foxd3 and zic3
    foxd3_clusters <- df %>% filter(tolower(gene) == "foxd3") %>% pull(cluster)
    zic3_clusters <- df %>% filter(tolower(gene) == "zic3") %>% pull(cluster)
    common_clusters <- intersect(foxd3_clusters, zic3_clusters)
    
    # Function to extract ligand genes (starts with pattern) in those clusters
          extract_ligands = function(pattern){
        df %>% 
          filter(cluster %in% common_clusters) %>%
          filter(str_starts(tolower(gene), pattern)) %>%
          select(cluster, gene)
      }
    
          wnt_genes <- extract_ligands("wnt")
          bmp_genes <- extract_ligands("bmp")
          
          # Add to ligand workbooks
          addWorksheet(wnt_ligands_wb, sheet_name)
          writeData(wnt_ligands_wb, sheet_name, wnt_genes)
          
          addWorksheet(bmp_ligands_wb, sheet_name)
          writeData(bmp_ligands_wb, sheet_name, bmp_genes)
  }
  
  # Save output Excel files
  saveWorkbook(nc_Candidates_wb, file.path(folder, "nc_candidates.xlsx"), overwrite = TRUE)
  saveWorkbook(wnt_ligands_wb, file.path(folder, "wnt_ligands.xlsx"), overwrite = TRUE)
  saveWorkbook(bmp_ligands_wb, file.path(folder, "bmp_ligands.xlsx"), overwrite = TRUE)
  
  message("All output files saved to ", folder)
}