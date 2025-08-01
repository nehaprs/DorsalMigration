
#################
#annotate the existing object
#################

library(GPTCelltype)

dors = dorsal
rm(dorsal)


DimPlot(dors, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.5", pt.size = 1) + ggtitle("UMAP Plot at resolution 0.5")
dors$seurat_clusters = dors$RNA_snn_res.0.5
table(Idents(dors))
Idents(dors) = dors$RNA_snn_res.0.5
table(Idents(dors))


dors$cell_type = as.character(Idents(dors))

#annotate neural crest cell types

seurat_obj$celltype[which(Idents(seurat_obj) == 3)] <- "Neural crest"
seurat_obj$celltype[which(Idents(seurat_obj) == 7)] <- "Placodal cells"

#run a clustify annotation using that file
