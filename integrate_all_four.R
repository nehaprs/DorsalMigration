#integrate_all_four

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

dorsals17 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage17/seurat_output/dorsals17.rds")
dorsals19 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage19/seurat_output/dorsals19.rds")
dorsals21 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage21/seurat_output/dors21_annot.rds")
dorsals24 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage24/seurat_output/dors24_annot.rds")


dorsals17 = RenameCells(dorsals17, add.cell.id = "s17")
dorsals19 = RenameCells(dorsals19, add.cell.id = "s19")
dorsals21 = RenameCells(dorsals21, add.cell.id = "s21")
dorsals24 = RenameCells(dorsals24, add.cell.id = "s24")


#set pre-integration cluster ids at desired resolution
Idents(dorsals17) = dorsals17$RNA_snn_res.1.3
dorsals17$cluster_orig = Idents(dorsals17)

Idents(dorsals19) = dorsals19$RNA_snn_res.1.5
dorsals19$cluster_orig = Idents(dorsals19)


Idents(dorsals21) = dorsals21$RNA_snn_res.2.6
dorsals21$cluster_orig = Idents(dorsals21)

Idents(dorsals24) = dorsals24$RNA_snn_res.3
dorsals24$cluster_orig = Idents(dorsals24)

seu_list = list(dorsals17, dorsals19, dorsals21, dorsals24)

features = SelectIntegrationFeatures(seu_list, nfeatures = 3000)
seu_list <- lapply(seu_list, \(x) ScaleData(x, features = features))
seu_list <- lapply(seu_list, \(x) RunPCA(x, features = features))


setwd("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/integrated/all4")

I2=  merge(seu_list[[1]], list(seu_list[[2]], seu_list[[3]], seu_list[[4]]))
table(I2$orig.ident)

#identify the batches
cn <- Cells(I2)
head(cn)
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

elbow = ElbowPlot(I2) #9 for 17-19 5 for 19-21 and 21-24. 10 for all 4
I2 = FindNeighbors(I2, dims = 1:10)
I2 = RunUMAP(I2, dims = 1:10)
I2 = FindClusters(I2)

p = DimPlot(I2, group.by = "orig_cluster", label = TRUE, label.size = 6, repel = TRUE)+ NoLegend()
ggsave("umapall4.png", plot = p, width = 17, height = 12, dpi = 300)

table(I2$stage)
DimPlot(I2, group.by = "stage", label = FALSE, label.size = 6, repel = TRUE)
saveRDS(I2,"all4.rds")

#color the clusters ide-d


clusters_keep <- c("s17_18", "s19_10", "s17_15", "s19_17", "s21_26", "s21_32",
                   "s21_33", "s24_33", "s24_28", "s24_29")


I2$highlight <- ifelse(
  I2$orig_cluster %in% clusters_keep,
  I2$orig_cluster,
  "Other"
)

# plot
DimPlot(
  I2,
  group.by = "highlight",
  cols = c(
    "s17_18" = "red",
    "s19_10" = "blue",
    "s17_15" = "green", "s19_17"= "yellow", "s21_26" = "pink" , "s21_32"="purple",
    "s21_33", "s24_33"="orange", "s24_28"="brown", "s24_29" = "cyan",
    "Other" = "grey80"
  )
)

###subset into 4 constituent clusters
s17 = subset(I2, subset = stage == "s17" )
s19 = subset(I2, subset = stage == "s19")
s21 = subset(I2, subset = stage == "s21")
s24 = subset(I2, subset = stage == "s24")

DimPlot(ss, group.by = "orig_cluster", label = TRUE,label.size = 5, repel = TRUE) + NoLegend() 

saveRDS(s17,"s17.rds")
saveRDS(s19,"s19.rds")
saveRDS(s21,"s21.rds")
saveRDS(s24,"s24.rds")


s17_s19 <- subset(I2, subset = stage %in% c("s17", "s19"))
s19_s21 = subset(I2, subset = stage %in% c("s19", "s21"))
s21_s24 = subset(I2, subset = stage %in% c("s21", "s24"))

saveRDS(s17_s19,"s17_s19.rds")
saveRDS(s19_s21, "s19_s21.rds")
saveRDS(s21_s24,"s21_s24.rds")
