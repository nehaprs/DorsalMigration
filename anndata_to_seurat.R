library(Seurat)

#library(SeuratData)
library(ggplot2)
library(readxl)
library(dplyr)
library(ggrepel)
library(scuttle)
s.query = dors

library(zellkonverter)
sce = readH5AD("~/BINF/scrnaseq general/dorsal migration/ref/sparsed_s18ref.h5ad")
assayNames(sce)       # X
colnames(sce)         # Cell names
rownames(sce)         # gene names in format Xetrov
assay(sce, "counts") <- assay(sce, "X")

sce <- logNormCounts(sce)
ref = as.Seurat(sce)


'library(SeuratDisk)

Convert("~/BINF/scrnaseq general/dorsal migration/ref/sparsed_s18ref.h5ad", dest = "h5seurat", overwrite = TRUE)
s.ref <- LoadH5Seurat("dataset.h5seurat")
'




