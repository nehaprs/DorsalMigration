#################
#Differentiate between rohon-beard and trigeminal markers. 
#################

library(dplyr)
library(readxl)
setwd("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/integrated/terminals")
' 
df = read_excel("BINF/scrnaseq general/dorsal migration/full head/highUMI/stage17/seurat_output/filt_markers/stage17_resolution_1.3.xlsx")
df <- read_excel("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage19/seurat_output/filt_markers/stage19_resolution_1.5.xlsx")
df <- read_excel("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage24/seurat_output/filtmarkers_resolution_3.xlsx")
df <- read_excel("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage21/seurat_output/filt_markers/stage21_highUMI-resolution_2.6.xlsx")
' 


df <- read_csv("s21-24/S21_26par_endpoint_summary.csv")
df <- read_csv("s21-24/S21_2par_endpoint_summary.csv")


rbgenes = c("cyp26b1", "gm53", "hoxa6", "hoxa7", "hoxa9", "hoxa10", "hoxb3", "hoxb5", "hoxb6",
            "hoxb7", "hoxc6", "hoxc8", "hoxc9", "hoxc10", "hoxd8", "hoxd9", "hoxd10", "kcnq5")

tggenes = c("avpr1a", "gabrd", "prlr")


genes_cl = df$gene[df$cluster == 2]


genes_in_rb = intersect(genes_cl, rbgenes)

write.table(
  genes_in_rb,
  file = "rb_st2426.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

genes_in_tg = intersect(genes_cl, tggenes)

write.table(genes_in_tg,
            file = "tg_st2426.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE
            )