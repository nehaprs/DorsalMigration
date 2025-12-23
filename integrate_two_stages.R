#=================================================
#created 12.11.2025

#12.11.2025 script: integration of stage 17 with 19 and monocle 


#=================================================

library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(hdf5r)
library(scCATCH)
library(tidyr)  
library(harmony)
library(SeuratWrappers)
library(monocle3)
library(SeuratObject)
library(scAnnoX)


setwd("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/integrated/s17-19")

dorsals17 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage17/seurat_output/dorsals17.rds")
dorsals19 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage19/seurat_output/dorsals19.rds")

dorsals17 = RenameCells(dorsals17, add.cell.id = "s17")
dorsals19 = RenameCells(dorsals19, add.cell.id = "s19")

head(colnames(dorsals17))

#set pre-integration cluster ids at desired resolution
Idents(dorsals17) = dorsals17$RNA_snn_res.1.3
dorsals17$cluster_orig = Idents(dorsals17)

Idents(dorsals19) = dorsals19$RNA_snn_res.1.5
dorsals19$cluster_orig = Idents(dorsals19)

seu_list = list(dorsals17, dorsals19)

features = SelectIntegrationFeatures(seu_list, nfeatures = 3000)
seu_list <- lapply(seu_list, \(x) ScaleData(x, features = features))
seu_list <- lapply(seu_list, \(x) RunPCA(x, features = features))

I2=  merge(seu_list[[1]], seu_list[[2]])
table(I2$orig.ident)

#identify the batches
cn <- Cells(I2)
last_char <- substring(cn, 1, 3)
table(last_char)

I2$stage = last_char

I2@meta.data = I2@meta.data %>%
  mutate(orig_cluster = paste(stage, cluster_orig, sep = "_"))

table(I2$orig_cluster)
########################################
##No harmony
########################################

DefaultAssay(I2) = "RNA"

I2 = NormalizeData(I2)
I2 = FindVariableFeatures(I2, nfeatures = 3000)
I2 = ScaleData(I2)
I2 = RunPCA(I2, npcs = 30)

elbow = ElbowPlot(I2) #9
I2 = FindNeighbors(I2, dims = 1:9)
I2 = RunUMAP(I2, dims = 1:9)
I2 = FindClusters(I2)

DimPlot(I2, group.by = "batch")


FeaturePlot(I2, features = c("foxd3", "sox10", "zic2", "zic3") )
saveRDS(I2,"mergdS17-19.rds")




##############################
#Part 2: pseudotime
##############################

I2[["RNA"]] <- JoinLayers(I2[["RNA"]])
cds = as.cell_data_set(I2)


colData(cds)$orig_cluster = I2$orig_cluster
colData(cds)$stage = I2$stage

reducedDims(cds)$UMAP <- Embeddings(I2, "umap")
reducedDims(cds)$PCA  <- Embeddings(I2, "pca")

cds = cluster_cells(cds, reduction_method = "UMAP")
cds = learn_graph(cds, use_partition = TRUE)
orig_cluster <- colData(cds)$orig_cluster
root_cells = rownames(colData(cds))[colData(cds)$stage == "s17"]
cds <- order_cells(cds, root_cells = root_cells)

plot_cells(cds, color_cells_by="orig_cluster",label_roots = TRUE,
           label_leaves = FALSE, label_branch_points = FALSE)







##############################
#part3: do clustering of I2
##############################

DimPlot(I2, group.by = "orig_cluster", label = TRUE, label.size = 3, repel = TRUE
        ) + NoLegend()
#s17 NCC: 12, 18. dm = 18
#st19 ncc: 6, 10, 20. dm: 10

VlnPlot(I2, group.by = 'orig_cluster', features = 'foxd3') + NoLegend()
dors.markers <- FindAllMarkers(I2, only.pos = TRUE)
write_xlsx(dors.markers, "markers_default.xlsx")

filt.makers = dors.markers %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0)
write_xlsx(filt.makers,"filt_markers_default.xlsx")



#color nc clusters
groups_to_highlight = c("s17_12", "s17_18", "s19_6", "s19_10", "s19_20")
I3 = I2
I3$highlight <- ifelse(I3$orig_cluster %in% groups_to_highlight, I3$orig_cluster, 
                       "other")

DimPlot(
  I3,
  reduction = "umap",
  group.by = "highlight",
  cols = c(
    "gray80",   
    "red",      
    "blue",  
    "green",
    "yellow", 
    "pink"
  ))

#variable resolutions
resolution.range <- seq(from =0, to = 3, by = 0.1)
resolution.range = 2.9
# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  I2<- FindClusters(I2, resolution = res)
  
  # Find all markers for the clusters at this resolution
  #I2.markers <- FindAllMarkers(I2, only.pos = TRUE)
  
  # Define the file name for saving the markers
  #file_name <- paste0("markers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  
  #write_xlsx(I2.markers, file_name)
  
  # Print a message to confirm completion for each resolution
  #print(paste("Markers for resolution", res, "saved to", file_name))
}


#list all xlsx files in wd

xlsx_file = list.files(pattern = "\\.xlsx$")
#xlsx_file = "markers_resolution_2.4.xlsx"
for (file in xlsx_file){
  df = read_xlsx(file)
  dff = df[df$avg_log2FC > 0,]
  dfff = dff[dff$p_val_adj < 0.05,]
  file_new = paste0("filt",file)
  write_xlsx(dfff, file_new)
}
DimPlot(
  I2, label = TRUE,
  reduction = "umap")


FeaturePlot(I2, features = c("foxd3", "sox10", "zic2", "zic3") )

#subcluster ncc and see what else is in there

table(I2$seurat_clusters)

I2 = autoAnnoTools(I2, 
                      method = 'scCATCH',
                      marker = marker_list,
                      cluster_col = "seurat_clusters",
                      strategy = 'marker-based'
)

DimPlot(I2, group.by = "scCATCH", pt.size = 1,label = TRUE) + ggtitle("Stage 24 res 4, software: scCATCH") + NoLegend()


#we have scCatch clusters. Save 'em

scCatch_clusters = I2@meta.data %>%
  select(seurat_clusters, scCATCH) %>%
  distinct()


write_xlsx(scCatch_clusters,"scCatchs17_19_res3.xlsx")
saveRDS(I2,"s17_19_scCatch_res3.rds")



I2$cluster40 = ifelse(I2$seurat_clusters == "40", "40", "other")
# Plot with cluster 40 in color and others in grey


DimPlot(I2, group.by = "cluster40", cols = c( "red", "grey"))

###################
#across clusters
#at lower res, what do dm nccs cluster with? which stage are they from?
###################



DimPlot(I2, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.2.9", pt.size = 1) + ggtitle("UMAP Plot, res:2.9")


#color nc clusters

I2$seurat_clusters = I2$RNA_snn_res.2.9
Idents(I2) = I2$RNA_snn_res.2.9
table(Idents(I2))

I2$cluster = ifelse(I2$seurat_clusters == "9", "9", "other")
# Plot with cluster 40 in color and others in grey


DimPlot(I2, group.by = "cluster", cols = c( "blue", "grey")) + ggtitle("cluster 9 at res 2.9")
table(I2$orig_cluster)

#how many s17_18 cells in each seurat_cluster?

s19_10 <- subset(I2, subset = orig_cluster == "s19_10")
table(s19_10$seurat_clusters)

write.table(
  table(s19_10$seurat_clusters),
  file = "clipboard",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


#scCatch annnotation


xentro_markers <- read_excel("~/BINF/scrnaseq general/dorsal migration/ref/xentro_markers.xlsx", 
                             col_types = c("text", "text", "skip", 
                                           "skip", "skip", "skip"))


marker_list = split(xentro_markers$Marker_genes, xentro_markers$State)

I2 = autoAnnoTools(I2, 
                   method = 'scCATCH',
                   marker = marker_list,
                   cluster_col = "seurat_clusters",
                   strategy = 'marker-based'
)

DimPlot(I2, group.by = "scCATCH", pt.size = 1,label = TRUE) + ggtitle("Stage 17+-19 res 2.9, software: scCATCH") + NoLegend()

scCatch_clusters = I2@meta.data %>%
  select(seurat_clusters, scCATCH) %>%
  distinct()


write_xlsx(scCatch_clusters,"scCatchs17-19_res4.xlsx")

############

I2$seurat_clusters = I2$RNA_snn_res.2.9
Idents(I2) = I2$RNA_snn_res.2.9
table(Idents(I2))


cells_s19_backtrace <- WhichCells(I2, idents = c(9,8,0,20,5,22,25,1,27))

# Build table of orig_cluster counts within these seurat clusters
orig_cluster_table <- I2@meta.data[cells_s19_backtrace, ] %>%
  dplyr::count(seurat_clusters, orig_cluster, name = "n_cells") %>%
  dplyr::arrange(seurat_clusters, dplyr::desc(n_cells))

orig_cluster_table

write_xlsx(orig_cluster_table,"orig_cluster_s19_backtrace-res2.9.xlsx")
