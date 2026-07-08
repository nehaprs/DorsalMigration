########seurat to anndata

library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(zellkonverter)
setwd("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/integrated/all4/high_res")


obj <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/integrated/all4/s21_s24.rds")
#obj <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/integrated/all4/high_res/s17_s19.rds")
obj@meta.data <- obj@meta.data[, 
                               !grepl("^RNA_snn_res\\.", colnames(obj@meta.data)) & 
                                 colnames(obj@meta.data) != "scCATCH"
]
obj <- JoinLayers(obj, assay = "RNA") #layers to be joined
DefaultAssay(obj) <- "RNA"
sce <- as.SingleCellExperiment(obj)
assayNames(sce)

writeH5AD(sce, file = "s21_24.h5ad", X_name = "counts")
