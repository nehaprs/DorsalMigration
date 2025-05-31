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

Idents(dorsal)
head(row.names(dorsal), n =50)
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

nc_sub12 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/subcluster/cluster12/nc_sub12.rds")
nc_sub11 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/subcluster/cluster11/nc_sub11.rds")
nc_sub20 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/subcluster/cluster20/nc_sub20.rds")
nc_suball3 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/subcluster/all3/nc_suball3.rds")

nc_sublist = c(nc_sub12, nc_sub11,nc_sub20,nc_suball3)

#convert from ensembl ids to gene names
source("~/GitHub/DorsalMigration/frog_name_ensembl2symbol.R")
nc_suball3 = frog_name_ensembl2symbol(nc_suball3)



for (nc in nc_sublist) {
  print(deparse(nc))
  obj_name <- deparse(substitute(nc))
  suffix   <- sub("^nc_sub", "", obj_name)}

#plot absolute expression of the genes of interest

nc = nc_suball3

  obj_name <- deparse(substitute(nc_suball3))
  suffix   <- sub("^nc_sub", "", obj_name)
  
  # 1) get the four feature plots as a list
  fp_list <- FeaturePlot(
    object    = nc,
    features  = c("zic2","zic3","pax3","msx1"),
    combine   = FALSE    # return a list instead of auto-combining
  )
  
  # 2) wrap into a single patchwork and 3) add a title
  p<- wrap_plots(fp_list, ncol = 2) +
    plot_annotation(
      title = paste0("Gene Expression Levels in Subcluster: ", suffix),
      theme = theme(plot.title = element_text(hjust = 0.5))
    )
  
  print(p)
rm(nc)

#get the number of cells in the subclusters


print(table(Idents(nc_suball3)))

#get the average expression level of each gene in each of the subclusters within
#these clusters.
## use the raw count so that they are comparable across different seurat objects


genes = c("zic2","zic3","pax3","msx1")
nc = nc_suball3


avg_list = AggregateExpression(nc, features = genes,
                             assays = "RNA",
                             slot = "data", #for raw expression levels
                             return.seurat = FALSE)

print(avg_list)


