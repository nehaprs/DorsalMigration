resolution.range <- seq(from = 0, to = 2, by = 0.1)

# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  allCombined<- FindClusters(allCombined, resolution = res)
  
  # Find all markers for the clusters at this resolution
  allCombined.markers <- FindAllMarkers(allCombined, only.pos = TRUE)
  
  # Define the file name for saving the markers
  file_name <- paste0("markers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  write_xlsx(allCombined.markers, file_name)
  
  # Print a message to confirm completion for each resolution
  print(paste("Markers for resolution", res, "saved to", file_name))
}


#list all xlsx files in wd

xlsx_file = list.files(pattern = "\\.xlsx$")

for (file in xlsx_file){
  df = read_xlsx(file)
  dff = df[df$avg_log2FC > 0,]
  dfff = dff[dff$p_val_adj < 0.05,]
  file_new = paste0("filt",file)
  write_xlsx(dfff, file_new)
}

allCombinedclust = clustree(allCombined)

markers_resolution_0_8 <- read_excel("~/BINF/yushi scrnaseq/all six/threshold0/harmony/monocle/markers_resolution_0.8.xlsx")

markers_resolution_0_8 = markers_resolution_0_8[markers_resolution_0_8$p_val_adj < 0.05,]
write_xlsx(markers_resolution_0_8,"~/BINF/yushi scrnaseq/all six/threshold0/harmony/monocle/filtered_markers_resolution_0_8.xlsx")

paxsox <- readRDS("~/BINF/yushi scrnaseq/all six/threshold0/harmony/monocle/allcombined_with_harmony.rds")

table(paxsox$RNA_snn_res.0.9)



e105Combined <- readRDS("C:/Users/neha/Documents/BINF/yushi scrnaseq/all six/threshold0/e105/e105Combined.rds")

DimPlot(e105Combined,reduction = "umap", group.by = "orig.ident")+
  ggtitle("Combined Dataset at E10.5 by Origin")
