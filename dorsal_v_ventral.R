library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(hdf5r)
library(scCATCH)
library(tidyr)
library(scAnnoX)

dors <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/integrated/all4/s17_s19.rds")
Idents(dors) = dors$orig_cluster
dors = JoinLayers(dors)
markers <- FindMarkers(
  dors,
  ident.1 = c("s17_18"),
  ident.2 = c("s19_20", "s19_6","s17_12", "s17_6", "s17_7"),
  logfc.threshold = 0.25,
  min.pct = 0.1
)

markers$gene = rownames(markers)

write_xlsx(markers,"s17d_v_allv.xlsx")


#filter for av_log_fx

xlsx_file = list.files(pattern = "\\.xlsx$")
#xlsx_file = "markers_resolution_2.4.xlsx"
for (file in xlsx_file){
  df = read_xlsx(file)
  dff = df[df$avg_log2FC > 0,]
  dfff = dff[dff$p_val_adj < 0.05,]
  file_new = paste0("filt",file)
  write_xlsx(dfff, file_new)
}

#