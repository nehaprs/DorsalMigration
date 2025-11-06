#=================================================
# run merge() without harmony.
# preserve the original standalone clusters
#started 10.16.2025
#=================================================


library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(hdf5r)
library(scCATCH)
library(tidyr)  
library(harmony)
setwd("~/BINF/scrnaseq general/dorsal migration/full head/merged/noHarmPT")
s1 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/st17force cells40k/seurat_output/dorsals17.rds")
s2 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/st19/seurat output/dorsals19.rds")
s3 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/st21/seurat_output/dorsals21.rds")
s4 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/st24-force-cells-100k/seurat_output/dorsals24.rds")

#make Idents = desired clusters for each of the object
Idents(s1) = s1$RNA_snn_res.0.6
Idents(s2) = s2$RNA_snn_res.2.9
Idents(s3) = s3$RNA_snn_res.2.4
Idents(s4) = s4$RNA_snn_res.3

lab <- list(
  s1 = setNames(paste0("s1_", as.character(Idents(s1))), Cells(s1)),
  s2 = setNames(paste0("s2_", as.character(Idents(s2))), Cells(s2)),
  s3 = setNames(paste0("s3_", as.character(Idents(s3))), Cells(s3)),
  s4 = setNames(paste0("s4_", as.character(Idents(s4))), Cells(s4))
)
lab_vec <- do.call(c, lab)

#cell names from the merged object
m <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/merged/mergNoHarm.rds")

if (grepl("[1-4]$", Cells(m)[1])) {
  nm <- names(lab_vec)
  
  append_src <- function(nm, pool, suf) {
    i <- nm %in% pool & !grepl("[1-4]$", nm)  # only unsuffixed names
    nm[i] <- paste0(nm[i], suf)
    nm
  }
  
  nm <- append_src(nm, Cells(s1), "1")
  nm <- append_src(nm, Cells(s2), "2")
  nm <- append_src(nm, Cells(s3), "3")
  nm <- append_src(nm, Cells(s4), "4")
  
  names(lab_vec) <- nm
}

names(lab_vec) <- sub("^(s[1-4])\\.(.*)", "\\2_\\1", names(lab_vec))
names(lab_vec) <- sub("_s", "_", names(lab_vec)) 
#names(lab_vec) <- sub("^s[1-4]\\.", "", names(lab_vec))

head(names(lab_vec))

m$orig_cluster <- factor(lab_vec[Cells(m)])
head(Cells(m))
table(m$orig_cluster)
head(lab_vec)

p = DimPlot(m, group.by = "orig_cluster") + NoLegend()

#look at s1_10


DimPlot(m, group.by = "orig_cluster",
        cells.highlight = which(m$orig_cluster == "s1_10"),
        cols.highlight = "red",
        cols = "grey") + NoLegend()


saveRDS(m, "ready_for_monocle.rds")
