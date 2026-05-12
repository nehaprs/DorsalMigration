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
  #ident.1 = c("s17_12", "s17_6", "s17_7"),
  #ident.2 = c("s19_20", "s19_6"),
  ident.1 = "s17_17",
  ident.2 = c("s19_12"),
  only.pos = FALSE,
  min.pct = 0.1
)

markers$gene = rownames(markers)

write_xlsx(markers,"s17_vs_19-dorsal lateral plate region.xlsx")

df = markers
df= df[abs(df$avg_log2FC) > 0.25,]
#df = dff[df$pct.1 - df$pct.2 > 0.1,]
df = df[df$p_val_adj < 0.05,]
write_xlsx(df,"filt_r1_s17_vs_19-dorsal lateral plate region.xlsx")

#filter for av_log_fx

xlsx_file = list.files(pattern = "\\.xlsx$")
#xlsx_file = "markers_resolution_2.4.xlsx"
for (file in xlsx_file){
  df = read_xlsx(file)
  dff = df[abs(df$avg_log2FC) > 0.25,]
  dfff = dff[dff$p_val_adj < 0.05,]
  file_new = paste0("filt",file)
  write_xlsx(dfff, file_new)
}

#

xlsx_file = list.files(pattern = "\\.xlsx$")

for (file in xlsx_file){
  df = read_xlsx(file)
  dff = df[abs(df$avg_log2FC) > 1,]
  dfff = dff[dff$pct.1 - dff$pct.2 > 0.1,]
  dfff = dff[dff$p_val_adj < 0.05,]
  file_new = paste0("filt",file)
  write_xlsx(dfff, file_new)
}