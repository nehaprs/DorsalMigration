#update.packages("spatstat.utils")
library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(hdf5r)

#setwd("~/BINF/scrnaseq general/dorsal migration/CR_count/outs/filtered_feature_bc_matrix")
'
s.data =  ReadMtx(mtx = "matrix.mtx.gz",
                  cells = "barcodes.tsv.gz",
                  features = "features.tsv.gz")
' 
s.data = Read10X_h5("~/BINF/scrnaseq general/dorsal migration/full head/CR-output/filtered_feature_bc_matrix.h5")

setwd("~/BINF/scrnaseq general/dorsal migration/full head")
dors = CreateSeuratObject(counts = s.data, project = "dorsal migration")
dors[["percent.mt"]] <- PercentageFeatureSet(dors, pattern = "^MT-")

vln = VlnPlot(dors, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dors.low = subset(dors, subset = nFeature_RNA< 350)
#3520 'cells' have nFeatures < 500
#4 cells < 250
#1018 < 350
FeatureScatter(dors, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

##check for empty droplets
#BiocManager::install("DropletUtils")
library(DropletUtils)
raw <- Read10X_h5("~/BINF/scrnaseq general/dorsal migration/full head/CR-output/raw_feature_bc_matrix.h5")
e.out <- emptyDrops(raw)
# Keep barcodes with FDR < 0.01 and 0.05
keep <- e.out$FDR < 0.05
sum(keep, na.rm=TRUE)
#3573 real cells on FDR 0.01
#3777 real cells on FDR 0.05

#sox9 <- subset(sox9, subset = nFeature_RNA > 2000 & nFeature_RNA < 7000 & percent.mt < 5)
dors = NormalizeData(dors)

dors = FindVariableFeatures(dors,selection.method = "vst" )

top10 = head(VariableFeatures(dors),10)

#scaling data
all.genes = rownames(dors)
dors = ScaleData(dors, features = all.genes)

#PCA
dors = RunPCA(dors, features = VariableFeatures(object = dors))
heat = DimHeatmap(dors, dims = 1:20, cells = 500, balanced = TRUE)
#pc 5 or 6, even that is a stretch though

elbow = ElbowPlot(dors)
#choose 14
dors = FindNeighbors(dors, dims = 1:14)


resolution.range <- seq(from = 0, to = 1, by = 0.1)

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
  dff = df[df$avg_log2FC > 1,]
  file_new = paste0("filt",file)
  write_xlsx(dff, file_new)
}

dorsclust = clustree(dors)
#choose 0.5

dors = RunUMAP(dors, dims = 1:14)
DimPlot(dors, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.6", pt.size = 1) + ggtitle("UMAP Plot")

saveRDS(dors,"dorsal.rds")
