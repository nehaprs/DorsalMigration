
#================================================================================
#analysis of the three clusters of interest within the larger seurat object 'dorsal'
#================================================================================

#cluster 11, 12, 20 seem nc-like
#suspects 12 and 20 includes dorsally migrating cells
#analysis to check this

library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(hdf5r)
library(scCATCH)
library(tidyr)

dorsal <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/dorsal.rds")

setwd("~/BINF/scrnaseq general/dorsal migration/full head/three clusters")

#feature plots of the genes of interest
FeaturePlot(object    = dorsal,
  features  = c("zic2","zic3","pax3","msx1","twist1"))

#find aggregate expression level of these genes 

genes = c("zic2","zic3","pax3","msx1","twist1")



avg_list = AggregateExpression(dorsal, features = genes,
                               assays = "RNA",
                               slot = "data", #for raw expression levels
                               return.seurat = FALSE)

print(avg_list)
write_xlsx(as.data.frame(avg_list),"aggregate_gene_expression.xlsx")
row.names(as.data.frame(avg_list))


########################
#paired analysis to find the markers:
# gene significantly higher in cluster 20/12 than in cluster 11
#find markers that are enriched in 20 or 12 compared to 11
pairedmarkers20vs11 = FindMarkers(dorsal,
                        ident.1 = 20,
                        ident.2 = 11,
                        only.pos       = TRUE,
                        return.thresh = 0.05)
pairedmarkers20vs11$gene = row.names(pairedmarkers20vs11)
write_xlsx(pairedmarkers20vs11, "pairedmarkers20vs11.xlsx")


pairedmarkers12vs11 = FindMarkers(dorsal,
                                  ident.1 = 12,
                                  ident.2 = 11,
                                  only.pos       = TRUE,
                                  return.thresh = 0.05)
pairedmarkers12vs11$gene = row.names(pairedmarkers12vs11)
write_xlsx(pairedmarkers12vs11, "pairedmarkers12vs11.xlsx")

###################

