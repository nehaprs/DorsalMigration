#finding common genes
setwd("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/integrated")

library(readxl)
library(dplyr)
library(purrr)
library(tidyr)

xlsx_path <- "clusters-of-interest.xlsx"

sheets <- excel_sheets(xlsx_path)

df_list <- setNames(
  map(sheets, ~ read_excel(xlsx_path, sheet = .x)),
  sheets
)

df_long <- imap_dfr(df_list, function(dat, sh) {
  dat %>%
    select(gene, p_val_adj, avg_log2FC) %>%
    mutate(
      gene = as.character(gene),
      sheet = sh
    ) %>%
    select(gene, sheet, p_val_adj, avg_log2FC)
})

df_long = df_long %>% filter(avg_log2FC > 1)

genes_all <- df_long %>%
  distinct(sheet, gene) %>%                
  count(gene, name = "n_sheets") %>%        
  arrange(desc(n_sheets)) 

genes_all2 <- df_long %>%
  distinct(sheet, gene) %>% 
  group_by(gene) %>%
  summarise(
    n_sheets = n(),
    sheets = paste(sort(sheet), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sheets))

writexl::write_xlsx(genes_all2, "top-genes-in-many-clusters.xlsx")

################
#common genes
setwd("~/BINF/scrnaseq general/dorsal migration/full head/highUMI")
s17 <- read_excel("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage17/seurat_output/filt_markers/stage17_resolution_1.3.xlsx")

s19 <- read_excel("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage19/seurat_output/filt_markers/stage19_resolution_1.5.xlsx")

s17_18 = s17[s17$cluster == 18,]
s19_10 = s19[s19$cluster == 10,]
dm_commons = inner_join(s17_18, s19_10, by = "gene")
writexl::write_xlsx(dm_commons,"comm_s1718-s1910.xlsx")
dm_comm_high = dm_commons[dm_commons$avg_log2FC.x > 1,]
dm_comm_high = dm_comm_high[dm_comm_high$avg_log2FC.y > 1,]
writexl::write_xlsx(dm_comm_high,"comm_s1718-s1910high.xlsx")

s17_12 = s17[s17$cluster == 12,]
s17_dmOnly = anti_join(s17_18, s17_12, by = "gene")
writexl::write_xlsx(s17_dmOnly,"s17_dmOnly.xlsx")

s19_6 = s19[s19$cluster == 6,]
s19_20 = s19[s19$cluster == 20,]
s19_dmonly =  s19_10 %>%
  anti_join(s19_6, by = "gene") %>%
  anti_join(s19_20, by = "gene")

writexl::write_xlsx(s19_dmonly,"s19_dmOnly.xlsx")

dm_only_comms = inner_join(s17_dmOnly, s19_dmonly, by = "gene")
writexl::write_xlsx(dm_only_comms,"dm_only_commsS1719.xlsx")
