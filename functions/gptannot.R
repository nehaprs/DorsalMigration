

library(GPTCelltype)
library(openai)

all.markers = filtmarkers_resolution_0_5
cl0.markers = all.markers[all.markers$cluster == 0,]
class(all.markers)
marker_list <- split(all.markers$gene, all.markers$cluster)
markerscl0 = split(cl0.markers$gene, cl0.markers$cluster)
head(markerscl0)

markerscl0.1 = markerscl0[1:1000]
res <- gptcelltype(markerscl0.1, 
                   tissuename = 'frog embryo', 
                   model = 'gpt-4'
)