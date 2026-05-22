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

dors <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/integrated/all4/high_res/s17_s19.rds")
Idents(dors) = dors$orig_cluster
dors = JoinLayers(dors)
markers <- FindMarkers(
  dors,
  #ident.1 = c("s17_11"),
  ident.2 = c("s19_4", "s19_22"),
  #ident.1 = "s17_16",
  ident.1 = c("s19_10"),
  only.pos = FALSE,
  min.pct = 0.1
)

markers$gene = rownames(markers)

write_xlsx(markers,"s19d_vs_19v.xlsx")

df = markers
df= df[abs(df$avg_log2FC) > 0.1,]
#df = dff[df$pct.1 - df$pct.2 > 0.1,]
df = df[df$p_val_adj < 0.05,]
write_xlsx(df,"filt_r1_s19d_vs_19v.xlsx")


df = markers
df= df[abs(df$avg_log2FC) > 0.1,]
df = df[df$pct.1 - df$pct.2 > 0.1,]
df = df[df$p_val_adj < 0.05,]
write_xlsx(df,"filt_r2_s19d_vs_19v.xlsx")

