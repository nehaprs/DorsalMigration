####
#does stage 17  and 19 require higher resolutions?
###

#try dotplots at current res
#choose higher res
#try diotplots again. is there a difference?

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(clustree)
library(writexl)

s17 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage17/seurat_output/dors17_annot.rds")
s19 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage19/seurat_output/dors19_annot.rds")


genes <- c("foxd3", "sox10", "sox9",
           "tfap2b",	"tfap2a",	"tfap2c",	"tfap2e",	"snai2",	"snai1",
           "sox8",	"twist1",	"adam33",	"mmp14",	"cdh2",	"itga5",	"mafb",
           "msx1",	"pax3",
             "zic2", "zic3")


p = DotPlot(s17, features = genes, dot.scale = 5)
p = p + 
  theme(
    axis.text.y = element_text(size = 10,angle = 45 ),
    axis.text.x = element_text(size = 15,angle = 45, hjust = 1 )
    #axis.title.y = element_text(size = 10)
  )

p = DotPlot(s19, features = genes, dot.scale = 5)
p = p + 
  theme(
    axis.text.y = element_text(size = 10,angle = 45 ),
    axis.text.x = element_text(size = 15,angle = 45, hjust = 1 )
    #axis.title.y = element_text(size = 10)
  )

clustree(s17) #try 1.8 for s17

Idents(s17) = s17$RNA_snn_res.1.8

p = DotPlot(s17, features = genes, dot.scale = 5)
p = p + 
  theme(
    axis.text.y = element_text(size = 10,angle = 45 ),
    axis.text.x = element_text(size = 15,angle = 45, hjust = 1 )
    #axis.title.y = element_text(size = 10)
  )


####
#higher resolutions
dors = s17

#use 
resolution.range <- seq(from = 2.1, to = 3, by = 0.1)

# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  dors<- FindClusters(dors, resolution = res)
  
  # Find all markers for the clusters at this resolution
  dors.markers <- FindAllMarkers(dors, only.pos = TRUE)
  
  # Define the file name for saving the markers
  file_name <- paste0("markers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  write_xlsx(dors.markers, file_name)
  
  # Print a message to confirm completion for each resolution
  print(paste("Markers for resolution", res, "saved to", file_name))
}


#list all xlsx files in wd

xlsx_file = list.files(pattern = "\\.xlsx$")

for (file in xlsx_file){
  df = read_xlsx(file)
  dff = df[df$avg_log2FC > 1,]
  dfff = dff[dff$p_val_adj < 0.05,]
  file_new = paste0("filt",file)
  write_xlsx(dfff, file_new)
}

source("~/GitHub/DorsalMigration/functions/clusters_and_markers.R")
clusters_and_markers("filt_markers")

source("~/GitHub/DorsalMigration/functions/clean_nc_clusters.R")

clean_nc_clusters("filt_markers/nc_Candidates.xlsx")

source("~/GitHub/DorsalMigration/functions/count_nc_clusters.R")
make_input_counted("nc_Candidates_reshaped.xlsx")

clustree(dors)

#try dotplot for 2.1, 2.2, 2.8

Idents(dors) = dors$RNA_snn_res.2.2 #fix 2.2 for further consideration

DotPlot(dors, features = genes)
saveRDS(dors,"stage17-res22.rds")

#stricter filetring for 2.2
filt1 <- read_excel("filt_markers/filtmarkers_resolution_2.2.xlsx")

filt2 = filt1[
  filt1$avg_log2FC > 0.5 &
  filt1$pct.1 - filt1$pct.2 > 0.1,
]

write_xlsx(filt2, "morefilt_st17_res2.2.xlsx")

########
#stage 19
setwd("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage19/seurat_output")

clustree(s19)
dors = s19


#use 
resolution.range <- seq(from = 2.1, to = 3, by = 0.1)

# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  dors<- FindClusters(dors, resolution = res)
  
  # Find all markers for the clusters at this resolution
  dors.markers <- FindAllMarkers(dors, only.pos = TRUE)
  
  # Define the file name for saving the markers
  file_name <- paste0("markers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  write_xlsx(dors.markers, file_name)
  
  # Print a message to confirm completion for each resolution
  print(paste("Markers for resolution", res, "saved to", file_name))
}


#list all xlsx files in wd

xlsx_file = list.files(pattern = "\\.xlsx$")

for (file in xlsx_file){
  df = read_xlsx(file)
  dff = df[df$avg_log2FC > 1,]
  dfff = dff[dff$p_val_adj < 0.05,]
  file_new = paste0("filt",file)
  write_xlsx(dfff, file_new)
}

source("~/GitHub/DorsalMigration/functions/clusters_and_markers.R")
clusters_and_markers("filt_markers")

source("~/GitHub/DorsalMigration/functions/clean_nc_clusters.R")

clean_nc_clusters("filt_markers/nc_Candidates.xlsx")

source("~/GitHub/DorsalMigration/functions/count_nc_clusters.R")
make_input_counted("nc_Candidates_reshaped.xlsx")

clustree(dors)

#try 2.3 and 2,4
Idents(dors) = dors$RNA_snn_res.2.3 #fix 2.3 for further consideration

DotPlot(dors, features = genes)

filt1 <- read_excel("filt_markers/filtmarkers_resolution_2.3.xlsx")
filt2 = filt1[
  filt1$avg_log2FC > 0.5 &
    filt1$pct.1 - filt1$pct.2 > 0.1,
]

write_xlsx(filt2, "morefilt_st19_res2.3.xlsx")
saveRDS(dors, "stage19-res23.rds")
