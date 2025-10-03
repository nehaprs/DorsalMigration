#################
#Downstream analysis of the obtained results
#created 25 sep 2025
##################

library(dplyr)

#dm NCC st 17 and 19: common genes

s17c10 <- read_excel("~/BINF/scrnaseq general/dorsal migration/full head/st17force cells40k/seurat_output/filt_markers/stage17_resolution_0.6.xlsx", 
                     +     sheet = "cluster 10")

s19c43 <- read_excel("~/BINF/scrnaseq general/dorsal migration/full head/st19/seurat output/filt_markers/stage19_resolution_2.9.xlsx", 
                     +     sheet = "cluster 43")

in17n19 = inner_join(s17c10, s19c43, by = "gene")
setwd("~/BINF/scrnaseq general/dorsal migration/full head/across_times")
writexl::write_xlsx(in17n19,"dmNCC_17n19.xlsx")


#######################################################
#Prioritizing markers
########################################################