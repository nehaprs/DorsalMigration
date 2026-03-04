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
