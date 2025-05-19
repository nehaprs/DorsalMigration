#BiocManager::install("clustifyr")
library(clustifyr)
library(dplyr)
library(Seurat)
library(writexl)
library(readxl)
library(tibble)
library(tidyr)
dors = readRDS("dorsal.rds")
xentro_markers <- read_excel("~/BINF/scrnaseq general/dorsal migration/ref/xentro_markers.xlsx", 
                               col_types = c("text", "text", "skip", 
                                                      "skip", "skip", "skip"))


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

DimPlot(dors2, group.by = "scCATCH", pt.size = 1,label = TRUE) + ggtitle("Cell Type Annotations, software: scCATCH") + NoLegend()


##############

#look for markers
neuralCrest.markers = c("Foxd3", "Snail2", "Sox9","Zic5", "Tfap2a")
FeaturePlot(dors, reduction = "umap", pt.size = 0.5,
            features = neuralCrest.markers)
