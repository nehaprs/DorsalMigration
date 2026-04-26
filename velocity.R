########
#convert seurat to AnnData

library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

library(SeuratDisk)

setwd("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/velocity")

#convert seurat v5 to h5ad object
#obj <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/integrated/s21-24/round2/mergdS21-24.rds")
obj = s21_s24
obj@meta.data <- obj@meta.data[, 
                               !grepl("^RNA_snn_res\\.", colnames(obj@meta.data)) & 
                                 colnames(obj@meta.data) != "scCATCH"
]

obj <- JoinLayers(obj, assay = "RNA")
DefaultAssay(obj) <- "RNA"
sce <- as.SingleCellExperiment(obj)
assayNames(sce)
# "counts"    "logcounts"
writeH5AD(sce, file = "s21_24.h5ad", X_name = "counts")


