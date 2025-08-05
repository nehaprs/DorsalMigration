
#################
#annotate the existing object
#################
setwd("~/BINF/scrnaseq general/dorsal migration/full head/version4/no CCScoring/annotation")
library(GPTCelltype)
library(write)
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



#run a clustify annotation using that file

#find which clusters are what under scCatch
dors = dors2
table(dors$type)
scCatch_annots = as.data.frame(table(dors2$scCATCH))
write_xlsx(scCatch_annots,"scCatch_annots_res0.5.xlsx")


# Step 2: Make a clean data frame
scCatch_table <- data.frame(
  Cluster = names(scCatch_annotations),
  CellType = unname(scCatch_annotations),
  stringsAsFactors = FALSE
)

#markergenes vlns
VlnPlot(dors, features = c("fezf1", "arx", "gdf10", "nkx2-1"))+plot_layout(ncol = 2) +  
  plot_annotation(title = "pre-chordal neural plate markers") 

VlnPlot(dors, features = c("rax")) + ggtitle("optic vesicle markers")

VlnPlot(dors, features = c("tbx1", "s1pr5"))+plot_layout(ncol = 2) +  
  plot_annotation(title = "otic placode markers") 

VlnPlot(dors, features = c("tp63", "foxi1e", "grhl1", "ctbs")) +
  plot_annotation(title = "Epidermal Progenitor Markers")
        

VlnPlot(dors, features = c("ocm3", "mybpc3", "ckb", "gapdh", "ankrd2", "myod1", "cap2", "myf6",  "actc1")) +
  plot_annotation(title = "Somite Markers")


VlnPlot(dors, features = c("slc5a8", "rbp4", "hebp2", "sox17a", "ca7", "poldip2", "gata4", "frzb2", "frzb", "sox17b")) +
  plot_annotation(title = "Endoderm Markers")


VlnPlot(dors, features = c("sox2", "sox3", "zeb2", "otx2", "zic1", "zic2", "zic3", "ncam1" )) +
  plot_annotation(title = "neural plate Markers")

VlnPlot(dors, features = c("aqp3", "pax3", "gbx2", "tfap2a", "zic1", "msx1")) +
  plot_annotation(title = "neural plate border Markers")


VlnPlot(dors, features = c("dlx1", "tfap2b", "dlx2", "rpe65", "sox10", "sox8", "twist1", "itga5",  "egr2")) +
  plot_annotation(title = "cranial neural crest Markers")

VlnPlot(dors, features = c("lhx9", "hoxd3", "mst1", "dhh",  "hoxd1")) +
  plot_annotation(title = "Hindbrain Markers")





VlnPlot(dors, features = c("foxc1", "lhx1", "pkd2", "osr2")) +
  plot_annotation(title = "Interediate Mesoderm Markers")

VlnPlot(dors, features = c("foxc1", "lhx1", "pkd2", "osr2")) +
  plot_annotation(title = "Interediate Mesoderm Markers")

VlnPlot(dors, features = c("foxd1", "fezf1", "otx2", "foxg1", "arx", "s1pr1","six3", "six6")) +
  plot_annotation(title = "Anterior Neural Tube Markers")

VlnPlot(dors, features = c("six3", "six6", "otx2", "pax6", "hhex")) +
  plot_annotation(title = "")

###########################################

#annotate using Anndata

##########################################

library(SeuratDisk)
#ref from anndata_to_seurat.R

table(ref$Celltype_reannotated)
Idents(ref) = "Celltype_reannotated"


# Find anchors
anchors = FindTransferAnchors(reference = ref, query = dors, dims = 1:20)
# Transfer labels
predictions <- TransferData(anchorset = anchors, refdata = Idents(ref), dims = 1:30)

# Add predictions to query
query <- AddMetaData(query, metadata = predictions)