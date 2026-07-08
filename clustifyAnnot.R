#BiocManager::install("clustifyr")
library(clustifyr)
library(dplyr)
library(Seurat)
library(writexl)
library(readxl)
library(tibble)
library(tidyr)
library(readxl)



dors <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage17/seurat_output/stage17-res22.rds")

table(dors$RNA_snn_res.2.2)
dors$seurat_clusters = dors$RNA_snn_res.2.2

xentro_markers <- read_excel("~/BINF/scrnaseq general/dorsal migration/ref/xentro_markers.xlsx", 
                               col_types = c("text", "text", "skip", 
                                                      "skip", "skip", "skip"))
Idents(dors) = dors$RNA_snn_res.2.2


#####################



marker_matrix <- xentro_markers %>%
  group_by(State) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = State, values_from = Marker_genes) %>%
  select(-row) %>%
  as.data.frame() 

dors$seurat_clusters

dors <- clustify_lists(
  input = dors,          
  marker = marker_matrix,      # Prepared marker matrix
  cluster_col = "seurat_clusters",  # Name of cluster column
  metric = "hyper",            # Hypergeometric test
  genome_n = 21000             # Approximate  genome size
)

DimPlot(dors, group.by = "type") + ggtitle("Cell Type Annotations S17 at 2.2")

######################################

#annotate with the same reference, use different annotating software

#devtools::install_github('XQ-hub/scAnnoX')
#BiocManager::install('CHETAH')
#BiocManager::install("CelliD")

#install.packages("RcppEigen")
#install.packages("Rcpp")
#install.packages("ggsci")
#install.packages("viridis")
#install.packages("tidyverse")
#devtools::install_github("PaulingLiu/scibet")
#devtools::install_github("zwj-tina/scibetR")


library(scAnnoX)

marker_list = split(xentro_markers$Marker_genes, xentro_markers$State)

# Run multiple annotation tools

dors2 = autoAnnoTools(dors, 
                     method = 'scCATCH',
                     marker = marker_list,
                     cluster_col = "seurat_clusters",
                     strategy = 'marker-based'
                     )

DimPlot(dors2, group.by = "scCATCH", pt.size = 1,label = TRUE, repel = TRUE, label.size = 8) + ggtitle("Cell Type Annotations, software: scCATCH") + NoLegend()
table(dors2$scCATCH)


cluster_celltype_table = dors2@meta.data %>%
  select(seurat_clusters, scCATCH) %>%
  distinct() %>%
  arrange(seurat_clusters)

dors2$RNA_snn_res.2.2
cluster_to_celltype <- dors2@meta.data %>%
  dplyr::count(RNA_snn_res.2.2, scCATCH, name = "n") %>%
  group_by(RNA_snn_res.2.2) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(RNA_snn_res.2.2)

unique(dors2$scCATCH)
##############

#look for markers
neuralCrest.markers = c("Foxd3", "Snail2", "Sox9","Zic5", "Tfap2a")
FeaturePlot(dors, reduction = "umap", pt.size = 0.5,
            features = neuralCrest.markers)

saveRDS(dors2,"s19-23-annotated.rds")

