#check mitochondrial contamination

dors <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage24/seurat_output/dorsals24-2.rds")

dors[["percent.mt"]] <- PercentageFeatureSet(dors, pattern = "^MT-|^mt|^cox|^nd|^cytb|^atp")


VlnPlot(
  dors,
  features =  "percent.mt"
)

