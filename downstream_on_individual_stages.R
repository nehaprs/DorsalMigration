###################
##DOWNSTREAM ANALYSIS
# 10.29.2025: further subclustering of stage 19
# 12.8.2025: stage 21: count umi of dm ncc
###################

library(dplyr)
library(Seurat)
library(writexl)
library(readxl)
library(clustree)
library(hdf5r)
library(tidyr)

setwd("~/BINF/scrnaseq general/dorsal migration/full head/st19/seurat output/more_cluster")
dors <- readRDS("C:/Users/neha/Documents/BINF/scrnaseq general/dorsal migration/full head/st19/seurat output/dorsals19.rds")

resolution.range <- seq(from = 2.1, to = 3, by = 0.1)

# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  dors<- FindClusters(dors, resolution = res)
  
  # Find all markers for the clusters at this resolution
  dors.markers <- FindAllMarkers(dors, only.pos = TRUE)
  
  # Define the file name for saving the markers
  file_name <- paste0("markers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  write_xlsx(dors.markers, file_name)
  
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


#####################
setwd("~/BINF/scrnaseq general/dorsal migration/full head/st21/seurat_output")

dors <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/st21/seurat_output/dorsals21.rds")

dors$seurat_clusters = dors$RNA_snn_res.2.4
table(Idents(dors))
Idents(dors) = dors$RNA_snn_res.2.4
table(Idents(dors))

#count UMI

UMI_cluster27 = subset(dors, idents = 27)$nCount_RNA
summary(UMI_cluster27)

library(ggplot2)

df <- data.frame(UMI_cluster27)

ggplot(df, aes(x = UMI_cluster27)) +
  geom_histogram(bins = 100) +
  labs(x = "UMI count (nCount_RNA)", y = "Number of cells")

nrow(df) #2186 total
sum(df$UMI_cluster27 < 501) #218 cells lost in umi filtering
