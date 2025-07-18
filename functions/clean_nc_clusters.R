'

cleans the excel sheet that lists nc and migratory genes in the cluster

'



clean_nc_clusters <- function(orig_file) {
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(tools)
  
  # Extract base name from original file and create output file name
  base_name <- file_path_sans_ext(basename(orig_file))
  new_file <- paste0(base_name, "_reshaped.xlsx")
  
  # Get all sheet names from the original file
  sheet_names <- excel_sheets(orig_file)
  
  # Create a new workbook to store reshaped data
  out_wb <- createWorkbook()
  
  for (sheet in sheet_names) {
    # Read sheet
    df <- read_excel(orig_file, sheet = sheet)
    
    # Check for required columns
    if (!all(c("cluster", "gene") %in% colnames(df))) {
      warning(paste("Skipping sheet", sheet, "— missing 'cluster' or 'gene' column"))
      next
    }
    
    # Reshape to wide format: one column per cluster
    df_wide <- df %>%
      mutate(cluster = as.character(cluster)) %>%
      group_by(cluster) %>%
      mutate(row = row_number()) %>%
      ungroup() %>%
      pivot_wider(names_from = cluster, values_from = gene)
    
    # Write reshaped sheet
    addWorksheet(out_wb, sheet)
    writeData(out_wb, sheet, df_wide)
  }
  
  # Save workbook
  saveWorkbook(out_wb, new_file, overwrite = TRUE)
  
  message("Reformatted data saved to: ", new_file)
}

