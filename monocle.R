#=====================================
#monocle pseudotime
#created 10.16.2025
#############################

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

setwd("~/BINF/scrnaseq general/dorsal migration/full head/merged/noHarmPT")
m <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/merged/noHarmPT/ready_for_monocle.rds")
#DefaultAssay(m)
#Layers(m)


m[["RNA"]] <- JoinLayers(m[["RNA"]])
#Layers(m)
cds = as.cell_data_set(m) #m is the output from merge_preserve)cluster.R

colData(cds)$orig_cluster = m$orig_cluster
colData(cds)$timepoint = m$batch

reducedDims(cds)$UMAP <- Embeddings(m, "umap")
reducedDims(cds)$PCA  <- Embeddings(m, "pca")


cds = cluster_cells(cds, reduction_method = "UMAP")
cds = learn_graph(cds, use_partition = TRUE)

     