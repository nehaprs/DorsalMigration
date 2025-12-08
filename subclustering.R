#================================================
#subclustering NCC populations of stages 21 and 24
#12.1.2025: stage 24, res 1.4 cluster 5
#================================================

library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(hdf5r)
library(tidyr)

#stage 24, res 1.4 cluster 5
#stage 21 res 1.7 clusters 4, 14
#stage 21 res 1.6 clusters 4

setwd("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage21/seurat_output/subcluster-r16-c4")

dors0 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage21/seurat_output/dorsals21.rds")


dors0$seurat_clusters = dors0$RNA_snn_res.1.6
table(Idents(dors0))
Idents(dors0) = dors0$RNA_snn_res.1.6
table(Idents(dors0))


dors = subset(dors0, idents = c("4"))
rm(dors0)

#process this cluster
vln = VlnPlot(dors, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(dors, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

dors = NormalizeData(dors)

dors = FindVariableFeatures(dors,selection.method = "vst" )

dors = ScaleData(dors)

#PCA
dors = RunPCA(dors, features = VariableFeatures(object = dors))
elbow = ElbowPlot(dors) 
#stage 24, res 1.4 cluster 5: 7
#stage 21, res 1.7 clusters 4 and 14: 5
#stage 21, res 1.6 clusters 4: 7
dors = FindNeighbors(dors, dims = 1:7)
resolution.range <- seq(from = 0, to = 0.5, by = 0.1)

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

dorsclust = clustree(dors)

dors = RunUMAP(dors, dims = 1:7)
DimPlot(dors, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.2", pt.size = 1) + ggtitle("UMAP Plot, res:0.2")

saveRDS(dors,"subcl-r16-c4.rds")
###########################################


