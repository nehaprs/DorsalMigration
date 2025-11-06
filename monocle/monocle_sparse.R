library(Seurat)
library(SeuratWrappers)   # as.cell_data_set
library(Matrix)
library(SingleCellExperiment)
library(monocle3)

# 1) Convert from Seurat
cds <- as.cell_data_set(m, assay = DefaultAssay(m))
rm(m)
# 2) Force sparse storage
assays(cds)$counts <- as(assays(cds)$counts, "dgCMatrix")
if ("logcounts" %in% assayNames(cds))
  assays(cds)$logcounts <- as(assays(cds)$logcounts, "dgCMatrix")

# 3) Remove any dense assay Monocle won't need
if ("exprs" %in% assayNames(cds)) assays(cds)[["exprs"]] <- NULL
if ("data"  %in% assayNames(cds)) assays(cds)[["data"]]  <- NULL

# 4) Optionally restrict to HVGs to cut memory
hvg <- VariableFeatures(m)
if (length(hvg) > 0) {
  keep <- intersect(rownames(cds), hvg)
  if (length(keep) > 0) cds <- cds[keep, ]
}

# 5) Proceed with lightweight settings
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)
cds <- learn_graph(
  cds,
  use_partition = TRUE,
  learn_graph_control = list(ncenter = 100)
)
saveRDS(cds,"cds_learned_graph.rds")
save_monocle_objects(cds,"~/binf/dorsal migration/dorsal migration")
