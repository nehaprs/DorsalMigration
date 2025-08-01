
#============================================================
#edits 6.9.2025: v2: without removing cells with < 500 features, using <250 features instead.
#including cell cycle stages: scaledata : variables to regress: cell cycle stages, nCounts
#scaling with variable features instead of all genes
#Commit named v2


#edits 6.12.2025: v2:  removing cells with < 500 features.
#including cell cycle stages: scaledata : variables to regress: cell cycle stages, nCounts
#scaling with variable features instead of all genes

#edit 6.26.2025: further study v2. find smaller clusters of nc cluster so that we can find the DM NC

#edit4: 7.17.2025: V4 what kind of filtering will bring back wnt-enriched dorsl nc cells?. 
#only nFeatures< 6000. No cell cycles
#============================================================

library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(hdf5r)
library(scCATCH)
library(tidyr)

#setwd("~/BINF/scrnaseq general/dorsal migration/CR_count/outs/filtered_feature_bc_matrix")

s.data = Read10X_h5("~/BINF/scrnaseq general/dorsal migration/full head/CR-output/filtered_feature_bc_matrix.h5")

setwd("~/BINF/scrnaseq general/dorsal migration/full head/version4/CCScoring")
dors = CreateSeuratObject(counts = s.data, project = "dorsal migration")
ncol(dors) #10882 cells without filtering
dors[["percent.mt"]] <- PercentageFeatureSet(dors, pattern = "^MT-")

vln = VlnPlot(dors, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dors.low = subset(dors, subset = nFeature_RNA< 500)

#3520 'cells' have nFeatures < 500
#4 cells < 250
#94 cells < 300
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

#dors2 <- subset(dors, subset = nFeature_RNA > 220 & nFeature_RNA < 6000 & percent.mt < 5)
dors2 <- subset(dors, subset =  nFeature_RNA < 6000 & percent.mt < 5)
ncol(dors2)
#10869 cells after nFeature_RNA < 6000: working with this on edit 4
#10881 cells after Feature_RNA > 220 & nFeature_RNA < 7000
#10875 cells after nFeature_RNA > 250 & nFeature_RNA < 7000
#10863 cells after nFeatures > 250 & nFeature_RNA < 6000
#10881 cells after Feature_RNA > 100 & nFeature_RNA < 7000
dors = dors2
rm(dors2)
dors = NormalizeData(dors)

dors = FindVariableFeatures(dors,selection.method = "vst" )

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

dors = CellCycleScoring(object = dors, s.features = s.genes_xt, g2m.features = g2m.genes_xt)
dors$CC.Difference = dors$S.Score - dors$G2M.Score
colnames(dors@meta.data)

dors = ScaleData(dors, vars.to.regress = c("nCount_RNA", "CC.Difference"))
#dors = ScaleData(dors)
top10 = head(VariableFeatures(dors),10)

#scaling data
#all.genes = rownames(dors)
#dors = ScaleData(dors, features = all.genes)

#PCA
dors = RunPCA(dors, features = VariableFeatures(object = dors))
heat = DimHeatmap(dors, dims = 1:20, cells = 2000, balanced = TRUE)
#pc 5 or 6, even that is a stretch though

elbow = ElbowPlot(dors)
#choose 8
dors = FindNeighbors(dors, dims = 1:14)


resolution.range <- seq(from = 0, to = 3, by = 0.1)

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


dors = RunUMAP(dors, dims = 1:14)
DimPlot(dors, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.1.6", pt.size = 1) + ggtitle("UMAP Plot, res:1.6")
#choose res 0.5 for subclustering
saveRDS(dors,"dorsal.rds")

dors$seurat_clusters = dors$RNA_snn_res.0.9
table(Idents(dors))
Idents(dors) = dors$RNA_snn_res.0.9
table(Idents(dors))

#in res5, cluster 3 is neural crest



#find marker genes from excel sheets for each resolution

source("~/GitHub/DorsalMigration/functions/clusters_and_markers.R")
clusters_and_markers("filt_markers")

source("~/GitHub/DorsalMigration/functions/clean_nc_clusters.R")

clean_nc_clusters("filt_markers/results/nc_Candidates.xlsx")




##finding scCatch clusters
dorsal <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/dorsal.rds")
dors = dorsal

xentro_markers <- read_excel("~/BINF/scrnaseq general/dorsal migration/ref/xentro_markers.xlsx", 
                             col_types = c("text", "text", "skip", 
                                           "skip", "skip", "skip"))

library(scAnnoX)

marker_list = split(xentro_markers$Marker_genes, xentro_markers$State)

# Run multiple annotation tools

dors2 = autoAnnoTools(dors, 
                      method = 'scCATCH',
                      marker = marker_list,
                      cluster_col = "seurat_clusters",
                      strategy = 'marker-based'
)

DimPlot(dors2, group.by = "scCATCH", pt.size = 1,label = TRUE) + ggtitle("Cell Type Annotations, software: scCATCH") + NoLegend()

#we have scCatch clusters
#try subdividing the clusters further

#use 
resolution.range <- seq(from = 0, to = 2, by = 0.1)

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

table(dors$RNA_snn_res.2)
DimPlot(dors, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.3.6", pt.size = 1) + ggtitle("UMAP Plot at resolution 3.6")


#number of cells in each cluster at desired res
dors$seurat_clusters = dors$RNA_snn_res.2.4
table(Idents(dors))
Idents(dors) = dors$RNA_snn_res.2.4
table(Idents(dors))



ncol(dors)


dors2 = autoAnnoTools(dors, 
                      method = 'scCATCH',
                      marker = marker_list,
                      cluster_col = "seurat_clusters",
                      strategy = 'marker-based'
)

DimPlot(dors2, group.by = "scCATCH", repel = TRUE, pt.size = 1,label = TRUE
    ) + 
  ggtitle("Cell Type Annotations, software: scCATCH, res:2") + NoLegend()


#print the scCatch annotations

Idents(dors2)
#res 2, non-annotated

head(dors2$seurat_clusters)

#Compute markers for every cluster

markers = FindAllMarkers(object = dors2,
                         only.pos = TRUE,
                         return.thresh = 0.05) %>% as.data.frame()

#per-cluster lookup table

cluster_ann = dors2@meta.data %>%
  as.data.frame() %>%
  select(cluster = seurat_clusters,
         cell_type = scCATCH) %>%
  distinct()

write_xlsx(cluster_ann,"scCatch_clusters.xlsx")

mergedmarkers = left_join(markers, cluster_ann, by = "cluster")
mergedmarkers = mergedmarkers[mergedmarkers$p_val_adj < 0.05,]
write_xlsx(mergedmarkers,"scCatchmarkersNoScores.xlsx")
saveRDS(dors2,"dorsal.rds")

cluster_counts = as.data.frame(table(dors2$seurat_clusters))
write_xlsx(cluster_counts,"cluster_counts.xlsx")

##############
##candidate clusters: 11, 20, 20

#find differential markers in these clusters

#cluster 20 v 11


#genes enriched in each cluster
markers11 = FindMarkers(dors2,
                        ident.1 = 11,
                        only.pos       = TRUE,
                        return.thresh = 0.05)
markers12 = FindMarkers(dors2,
                        ident.1 = 12,
                        only.pos       = TRUE,
                        return.thresh = 0.05)

markers20 = FindMarkers(dors2,
                        ident.1 = 20,
                        only.pos       = TRUE,
                        return.thresh = 0.05)


#Genes enriched in 20 but not in 11
#genes enriched in 12 not 11
in12not11genes = setdiff(rownames(markers12), rownames(markers11))
#Subset the cluster-12 markers to just those

in12not11 = markers12[in12not11genes, ,drop = FALSE] #drop FALSE prevents R from simplifying the result to a vector or matrix if we end up selecting only one column
#Rank by log2 fold-change
in12not11 = in12not11 %>% arrange(desc(avg_log2FC))
in12not11 = in12not11[in12not11$p_val_adj < 0.05,]
GeneName = rownames(in12not11)
in12not11 = cbind( GeneName, in12not11)
write_xlsx(in12not11,"enrichment/in12not11.xlsx")


#genes enriched in 11 not 20
in11not20genes = setdiff(rownames(markers11), rownames(markers20))
#Subset the cluster-12 markers to just those

in11not20 = markers11[in11not20genes, ,drop = FALSE] #drop FALSE prevents R from simplifying the result to a vector or matrix if we end up selecting only one column
#Rank by log2 fold-change
in11not20 = in11not20 %>% arrange(desc(avg_log2FC))
in11not20 = in11not20[in11not20$p_val_adj < 0.05,]
gene = rownames(in11not20)
in11not20 = cbind( gene, in11not20)
write_xlsx(in11not20,"enrichment/in11not20.xlsx")


#genes enriched in 11 not 12
in11not12genes = setdiff(rownames(markers11), rownames(markers12))
#Subset the cluster-12 markers to just those

in11not12 = markers11[in11not12genes, ,drop = FALSE] #drop FALSE prevents R from simplifying the result to a vector or matrix if we end up selecting only one column
#Rank by log2 fold-change
in11not12 = in11not12 %>% arrange(desc(avg_log2FC))
in11not12 = in11not12[in11not12$p_val_adj < 0.05,]
gene = rownames(in11not12)
in11not12 = cbind( gene, in11not12)
write_xlsx(in11not12,"enrichment/in11not12.xlsx")

#################

#subclustering beyond 2

#use 
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
  dff = df[df$avg_log2FC > 1,]
  dfff = dff[dff$p_val_adj < 0.05,]
  file_new = paste0("filt",file)
  write_xlsx(dfff, file_new)
}

dorsclust = clustree(dors)

res3nos = as.data.frame(table(dors$RNA_snn_res.3))
write_xlsx(res3nos,"cell_numberres3.xlsx")
DimPlot(dors, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.3", pt.size = 1) + ggtitle("UMAP Plot at resolution 3")
###
#continued insubclustering.R
#subclustering for version2: res0.5, cluster3: in nc_clusters/subclusterv2.R
