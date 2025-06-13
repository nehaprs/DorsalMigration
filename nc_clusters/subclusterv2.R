#=========================================================================
'
subclustering for v2: that is, 6000 < features <250, cell cycle regressed
nc subcluster: res 0.5, cluster 3

subclustering for v3: that is, 6000 < features < 500. cell cycle regressed
nc subcluster: res 0.9, clusters 0,4,14
'
#==========================================================================

library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(scCATCH)
library(tidyr)

setwd("~/BINF/scrnaseq general/dorsal migration/full head/version3/subcluster")
dors <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/version3/dorsal.rds")
table(Idents(dors))

nc_sub = subset(dors, idents = c(0,4,14))
rm(dors)

nc_sub[["percent.mt"]] <- PercentageFeatureSet(nc_sub, pattern = "^MT-")

vln = VlnPlot(nc_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#no filtering on nFeatures
nc_sub = NormalizeData(nc_sub)
nc_sub = FindVariableFeatures(nc_sub, selection.method = "vst")

#-----------------
#cell cycle stages
#-----------------

# A list of HUMAN cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes




# Basic function to convert human to xentrop gene names

source('~/GitHub/DorsalMigration/humanGeneName2XenTropNames.R')
s.genes_xt = humanGeneName2XenTropNames(s.genes)
g2m.genes_xt = humanGeneName2XenTropNames(g2m.genes)

#calculate cell cycle phase score

nc_sub = CellCycleScoring(object = nc_sub, s.features = s.genes_xt, g2m.features = g2m.genes_xt)
nc_sub$CC.Difference = nc_sub$S.Score - nc_sub$G2M.Score
colnames(nc_sub@meta.data)

nc_sub = ScaleData(nc_sub, vars.to.regress = c("nCount_RNA", "CC.Difference"))


nc_sub = RunPCA(nc_sub, features = VariableFeatures(object = nc_sub))
heat = DimHeatmap(nc_sub, dims = 1:20, cells = 2000, balanced = TRUE)


elbow = ElbowPlot(nc_sub)
#elbow at ~6


nc_sub = FindNeighbors(nc_sub, dims = 1:6)


resolution.range <- seq(from = 0, to = 1, by = 0.1)

# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  nc_sub<- FindClusters(nc_sub, resolution = res)
  
  # Find all markers for the clusters at this resolution
  nc_sub.markers <- FindAllMarkers(nc_sub, only.pos = TRUE)
  
  # Define the file name for saving the markers
  file_name <- paste0("markers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  write_xlsx(nc_sub.markers, file_name)
  
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

nc_subclust = clustree(nc_sub)

#res 0.7 or 0.8
#v3: start with res 0.5
Idents(nc_sub) = nc_sub$RNA_snn_res.0.5
nc_sub$seurat_clusters = nc_sub$RNA_snn_res.0.5

nc_sub = RunUMAP(nc_sub, dims = 1:6)
DimPlot(nc_sub, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.5", pt.size = 1) + ggtitle("UMAP Plot subclustered NC cluster, res:0.5")

#diff b/w res 0.6 and 0.7 is some cells shuffled between clusters 7 and 8
#use res 0.7
saveRDS(nc_sub, "nc_sub.rds")

genes = c("zic2","zic3","pax3","msx1","cdh1","cdh2",
          "foxd3", "sox10", "snai2", "snai1", "tfap2a","twist1")

avg_list = AggregateExpression(nc_sub, features = genes,
                               assays = "RNA",
                               slot = "data", #for raw expression levels
                               return.seurat = FALSE)
avg_list$genes = rownames(as.data.frame(avg_list))
write_xlsx(as.data.frame(avg_list),"aggregateExpression.xlsx")
FeaturePlot(nc_sub, features = genes)



genes2 = c("nr4a1",
           "nr4a2",
           "junb",
           "irf1",
           "fosl1",
           "jun",
           "cbx4"
)


