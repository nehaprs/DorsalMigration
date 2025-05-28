#higher esolution seem to not subcluster clusters 11 and 20. do actual subclustering
library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(hdf5r)
library(scCATCH)
library(tidyr)

Idents(dors2)
# Subset to only clusters 11, 12, and 20
setwd("~/BINF/scrnaseq general/dorsal migration/full head/subcluster/all3")
nc_sub11 = subset(dors2, idents= c("11","12","20"))

nc_sub11 = NormalizeData(nc_sub11)

nc_sub11 = FindVariableFeatures(nc_sub11,selection.method = "vst" )

top10 = head(VariableFeatures(nc_sub11),10)

#scaling data
all.genes = rownames(nc_sub11)
nc_sub11 = ScaleData(nc_sub11, features = all.genes)

#PCA
nc_sub11 = RunPCA(nc_sub11, features = VariableFeatures(object = nc_sub11))
heat = DimHeatmap(nc_sub11, dims = 1:20, cells = 500, balanced = TRUE)
#pc 5 or 6, even that is a stretch though

elbow = ElbowPlot(nc_sub11)

#all3: 5, cl12:3, cl20:5
nc_sub11 = FindNeighbors(nc_sub11, dims = 1:5)


resolution.range <- seq(from = 0, to = 1, by = 0.1)

# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  nc_sub11<- FindClusters(nc_sub11, resolution = res)
  
  # Find all markers for the clusters at this resolution
  nc_sub11.markers <- FindAllMarkers(nc_sub11, only.pos = TRUE)
  
  # Define the file name for saving the markers
  file_name <- paste0("markers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  write_xlsx(nc_sub11.markers, file_name)
  
  # Print a message to confirm completion for each resolution
  print(paste("Markers for resolution", res, "saved to", file_name))
}


#list all xlsx files in wd

xlsx_file = list.files(pattern = "\\.xlsx$")

for (file in xlsx_file){
  df = read_xlsx(file)
  dff = df[df$avg_log2FC > 1,]
  dfff = dff[dff$p_val_adj < 0.05,]
  file_new = paste0("filt",file)
  write_xlsx(dfff, file_new)
}

nc_sub11clust = clustree(nc_sub11)
#choose 0.5

nc_sub11 = RunUMAP(nc_sub11, dims = 1:5)
DimPlot(nc_sub11, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.8", pt.size = 1) + ggtitle("UMAP Plot subclustered all 3 NC, res:0.8")

saveRDS(nc_sub11,"nc_suball3.rds")
#subcluster 3 of 20 looks very promising
#subclustering of 20: res: 0.9

