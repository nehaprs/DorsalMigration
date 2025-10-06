#=================================================
#created 10.1.2025
#DM NCC Pseudotime starting with stage 17 dm ncc
#10.1.2025: integrated the 4 datasets. with merge here, integratedata() in biomix
#10.2.2025: quantify batch effect correction if we need to do harmony.

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

setwd("~/BINF/scrnaseq general/dorsal migration/full head/merged")
#load seurat objects

dorsals17 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/st17force cells40k/seurat_output/dorsals17.rds")
dorsals19 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/st19/seurat output/dorsals19.rds")
dorsals21 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/st21/seurat_output/dorsals21.rds")
dorsals24 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/st24-force-cells-100k/seurat_output/dorsals24.rds")


#rm( dorsals19, dorsals21, dorsals24)

#before integration, create metadata of original clusters

#make cell names unique across samples
s17 = RenameCells(dorsals17, add.cell.id = "s17")
s19 = RenameCells(dorsals19, add.cell.id = "s19")
s21 = RenameCells(dorsals21, add.cell.id = "s21")
s24 = RenameCells(dorsals24, add.cell.id = "s24")

#set pre-integration cluster ids at desired resolution for future root cluster

Idents(dorsals17) = dorsals17$RNA_snn_res.0.6 #dm ncc cluster 10
dorsals17$cluster_pre_s1 = Idents(dorsals17)

seu_list = list(dorsals17, dorsals19, dorsals21, dorsals24)

#variable features were selected with the same parameters. go directly to select integration features

#first do RPCA based integration
#better than CCA for similar datasets. CCA can overcorrect.

features = SelectIntegrationFeatures(seu_list, nfeatures = 3000) #2000 is default, but 3k is recommended.


seu_list <- lapply(seu_list, \(x) ScaleData(x, features = features))
seu_list <- lapply(seu_list, \(x) RunPCA(x, features = features))
#anchors <- FindIntegrationAnchors(seu_list, anchor.features = features,
                                  reduction = "rpca", dims = 1:30)
#I1 <- IntegrateData(anchors, dims = 1:30)
#saveRDS(I1, "integrated.R")

#do merge

I2=  merge(seu_list[[1]], c(seu_list[[2]], seu_list[[3]], seu_list[[4]]))
saveRDS(I2,"mergedDM.rds")

rm(dorsals17, dorsals19, dorsals21, dorsals24)
rm(seu_list)



table(I2$orig.ident)

cn <- Cells(I2)
last_char <- substring(cn, nchar(cn), nchar(cn))
table(last_char)  # should show counts for 1,2,3,4 only
#topifnot(all(last_char %in% c("1","2","3","4")))  # fail early if not true



I2$batch = last_char


########################################
##try first without harmony
########################################

DefaultAssay(I2) #gives RNA. If not set RNA

I2 = NormalizeData(I2)
I2 = FindVariableFeatures(I2, nfeatures = 3000)
I2 = ScaleData(I2)
I2 = RunPCA(I2, npcs = 30)

elbow = ElbowPlot(I2) #choose 6 pcs
I2 = FindNeighbors(I2, dims = 1:6)
I2 = RunUMAP(I2, dims = 1:6)
I2 = FindClusters(I2)

DimPlot(I2, group.by = "batch")
