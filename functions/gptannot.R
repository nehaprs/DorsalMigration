Sys.setenv(OPENAI_API_KEY = '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')

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