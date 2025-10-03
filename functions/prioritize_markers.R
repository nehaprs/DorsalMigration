###
#from the results of findallmarkers, find the markers we need to prioritize
###

library(dplyr)
library(readxl)
library(writexl)

#input: output file from findAllMarkers

#filter for biological relevance

markers = read_excel("~/BINF/scrnaseq general/dorsal migration/full head/st19/seurat output/filt_markers/stage19_resolution_2.9.xlsx", 
                    )

keepers = markers %>%
  filter(p_val_adj < 0.05) %>%
  mutate(
    spec = pct.1 - pct.2
  ) %>%
  filter(avg_log2FC >= 0.2,
         pct.1 >= 0.01,
        spec > 0.01) %>%
  #filter low-informative genes
  filter(!(grepl("^(RPL|RPS|MT-)|^MALAT1$|^XIST$", gene)))
#write_xlsx(keepers,"s17_prioritized.xlsx")

#foxd3: finding nc clusters

#5,7,13



exclude = keepers %>% filter(cluster %in% c(30,37
))
dmncc = keepers %>% filter(cluster == 43) %>%
  anti_join(exclude, by = "gene")

write_xlsx(dmncc,"dmncc_unique_st19.xlsx")

################################
##common 17_unique and 19_unique
###############################

dmncc_unique_st17 <- read_excel("dmncc_unique_st17.xlsx")
dmncc_unique_st19 <- read_excel("dmncc_unique_st19.xlsx")

dmncc_unique_17n19 = inner_join(dmncc_unique_st17, dmncc_unique_st19, by = "gene")
write_xlsx(dmncc_unique_17n19, "dmncc_unique_17n19.xlsx")
