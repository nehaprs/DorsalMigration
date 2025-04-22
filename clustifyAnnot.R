BiocManager::install("clustifyr")
library(clustifyr)
library(dplyr)
library(Seurat)
library(writexl)
library(readxl)
library(tibble)
dors = readRDS("dorsal.rds")
xentro_markers <- read_excel("~/BINF/scrnaseq general/dorsal migration/ref/xentro_markers.xlsx", 
                               col_types = c("text", "text", "skip", 
                                                      "skip", "skip", "skip"))

# Create a named list from marker table
marker_list = split(xentro_markers$Marker_genes, xentro_markers$State)

# Assume your Seurat object is named seurat_obj
expr_mat <- GetAssayData(seurat_obj, slot = "data", assay = "RNA")
cluster_ids <- Idents(seurat_obj)

expr_mat = GetAssayData(dors, layer = "data", assay = "RNA")
cluster_ids = Idents(dors)
cluster_ids_fixed <- setNames(as.character(cluster_ids), names(cluster_ids))



metadata_df <- data.frame(
  cell_barcode = names(cluster_ids_fixed),
  cluster = cluster_ids_fixed,
  stringsAsFactors = FALSE
)

markerdf = as.data.frame(marker_list)
marker_scores <- clustify_lists(
  input = expr_mat,
  metadata = metadata_df,
  cluster_col = "cluster",     # name of the column in metadata
  marker = markerdf
)
marker_scores = clustify_lists(input = expr_mat,
                               cluster_col = cluster_ids,
                               markers = marker_list)
#####################



marker_matrix <- xentro_markers %>%
  group_by(State) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = State, values_from = Marker_genes) %>%
  select(-row) %>%
  as.data.frame() 


dors <- clustify_lists(
  input = dors,          
  marker = marker_matrix,      # Prepared marker matrix
  cluster_col = "seurat_clusters",  # Name of cluster column
  metric = "hyper",            # Hypergeometric test
  genome_n = 21000             # Approximate  genome size
)

DimPlot(dors, group.by = "type") + ggtitle("Cell Type Annotations")
