library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

s17 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage17/seurat_output/dors17_annot.rds")
s19 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage19/seurat_output/dors19_annot.rds")
s21 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage21/seurat_output/dors21_annot.rds")
s24 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage24/seurat_output/dorsals24-2.rds")
obj_list <- list(
  s17 = s17,
  s19 = s19,
  s21 = s21,
  s24 = s24
)

# Genes to plot

genes <- c("foxd3", "sox10", "sox9", "zic1", "zic2", "zic3", "zic4")


cluster_map <- data.frame(
  object_name = c("s17", "s19", "s21", "s24"),
  cluster_id  = c(18, 10),
  row_label   = c("s17-18", "s17-15", "s19-10", "s19-17" ),
  stringsAsFactors = FALSE
)

# -----------------------------
# Function to calculate dotplot values
# -----------------------------
get_dotplot_data <- function(seu, cluster_id, genes, row_label) {
  
  cells_use <- WhichCells(seu, idents = cluster_id)
  
  if (length(cells_use) == 0) {
    stop(paste("No cells found for cluster", cluster_id, "in", row_label))
  }
  
  expr <- FetchData(seu, vars = genes, cells = cells_use)
  
  data.frame(
    gene = genes,
    avg_exp = colMeans(expr),
    pct_exp = colMeans(expr > 0) * 100,
    cluster = row_label
  )
}

# -----------------------------
# Build plotting dataframe
# -----------------------------
plot_df <- lapply(seq_len(nrow(cluster_map)), function(i) {
  obj_name <- cluster_map$object_name[i]
  clust    <- cluster_map$cluster_id[i]
  label    <- cluster_map$row_label[i]
  
  get_dotplot_data(
    seu = obj_list[[obj_name]],
    cluster_id = clust,
    genes = genes,
    row_label = label
  )
}) %>%
  bind_rows()

# scale average expression by gene for better visual comparison
plot_df <- plot_df %>%
  group_by(gene) %>%
  mutate(avg_exp_scaled = scale(avg_exp)[, 1]) %>%
  ungroup()

# Keep row/column order exactly as desired
plot_df$cluster <- factor(plot_df$cluster, levels = rev(cluster_map$row_label))
plot_df$gene <- factor(plot_df$gene, levels = genes)

# -----------------------------
# Dot plot
# -----------------------------
ggplot(plot_df, aes(x = gene, y = cluster)) +
  geom_point(aes(size = pct_exp, color = avg_exp_scaled)) +
  scale_size(range = c(0, 10), name = "% expressing") +
  scale_color_gradient(low = "lightgrey", high = "blue", name = "Avg expression\n(scaled)") +
  theme_bw() +
  labs(x = "Genes", y = "Clusters") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

FeaturePlot(s24, features = c("foxd3", "sox10", "sox9", "zic1", "zic2", "zic3", "zic4"))


dm_genes <- c("foxd3", "sox10", "sox9", "zic2", "zic3")
pn_genes = c("isl1","isl2","ebf2","olfm1","pou4f4","tubb2b","runx3",
             "runx1","stmn2","cbfb","prph")
dm_genes <- c("foxd3", "sox10", "sox9","tfap2b", "tfap2a","tfap2c","tfap2e",
              "snai1","snai2","sox8","twist1","adam33","mmp14","cdh2","itga5","mafb",
              "msx1",
              "zic2", "zic3")
#FeaturePlot(s19, features = pn_genes)
Idents(s19) = s19$RNA_snn_res.2.3
table(Idents(s19))
p = DotPlot(s19, features = dm_genes, dot.scale = 5)
p = p + 
  theme(
    axis.text.y = element_text(size = 10,angle = 45 ),
    axis.text.x = element_text(size = 15,angle = 45, hjust = 1 )
    #axis.title.y = element_text(size = 10)
  ) + ggtitle("stage 19 at resolution 3")

ggsave("st19res3.png",p, width = 13, height = 8)
########################################

df17 <- read_excel("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage17/seurat_output/filt_markers/stage17_resolution_1.3.xlsx")
df19 <- read_excel("BINF/scrnaseq general/dorsal migration/full head/highUMI/stage19/seurat_output/filt_markers/stage19_resolution_1.5.xlsx")
df21 <- read_excel("BINF/scrnaseq general/dorsal migration/full head/highUMI/stage21/seurat_output/filt_markers/stage21_highUMI-resolution_2.6.xlsx")
df24 <- read_excel("BINF/scrnaseq general/dorsal migration/full head/highUMI/stage24/seurat_output/filtmarkers_resolution_3.xlsx")

#markers w/o filtering
df17 <- read_excel("BINF/scrnaseq general/dorsal migration/full head/highUMI/stage17/seurat_output/markers/markers_resolution_1.3.xlsx")
df19 <- read_excel("BINF/scrnaseq general/dorsal migration/full head/highUMI/stage19/seurat_output/markers_resolution_1.5.xlsx")
df21 <- read_excel("BINF/scrnaseq general/dorsal migration/full head/highUMI/stage21/seurat_output/markers_resolution_2.6.xlsx")
df24 <- read_excel("BINF/scrnaseq general/dorsal migration/full head/highUMI/stage24/seurat_output/markers_resolution_3.xlsx")



library(dplyr)
library(tidyr)

# -----------------------------
# INPUTS
# -----------------------------
dm_genes <- c("foxd3", "sox10", "sox9", "zic1", "zic2", "zic3", "zic4")
pn_genes = c("isl1","isl2","ebf2","olfm1","pou4f4","tubb2b",
             "stmn2","cbfb","prph")


# clusters of interest
cluster_map <- data.frame(
  dataset = c("df17", "df19","df21","df21","df21","df24","df24","df24"),
  cluster = c(15, 17,26,32,33,33,28,29),
  row_name = c("st17-cl15",  "st19-cl17","st21-cl26","st21-cl32","st21-cl33","st24-cl33","st24-cl28","st24-cl29"),
  stringsAsFactors = FALSE
)


# clusters of interest
cluster_map <- data.frame(
  dataset = c("df17", "df19"),
  cluster = c(18,10),
  row_name = c("st17-cl18",  "st19-cl10"),
  stringsAsFactors = FALSE
)

# put dfs in a list
df_list <- list(df17 = df17, df19 = df19, df21= df21, df24=df24)
df_list <- list(df17 = df17, df19 = df19)

# -----------------------------
# BUILD MATRIX
# -----------------------------
result <- lapply(seq_len(nrow(cluster_map)), function(i) {
  
  df_name <- cluster_map$dataset[i]
  clust   <- cluster_map$cluster[i]
  row_lab <- cluster_map$row_name[i]
  
  df <- df_list[[df_name]]
  
  # genes present in this cluster
  genes_in_cluster <- df %>%
    filter(cluster == clust) %>%
    pull(gene) %>%
    unique()
  
  # create row with checkmarks
  data.frame(
    row = row_lab,
    gene = pn_genes,
    value = ifelse(pn_genes %in% genes_in_cluster, "✓", ""),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  pivot_wider(names_from = gene, values_from = value)

# -----------------------------
# FINAL MATRIX
# -----------------------------
result

###############
#HEATMAP
###############

# -----------------------------
# BUILD HEATMAP DATA
# -----------------------------
plot_df <- lapply(seq_len(nrow(cluster_map)), function(i) {
  
  df_name <- cluster_map$dataset[i]
  clust   <- cluster_map$cluster[i]
  row_lab <- cluster_map$row_name[i]
  
  df <- df_list[[df_name]]
  
  # keep only requested cluster and genes
  sub_df <- df %>%
    filter(cluster == clust, gene %in% dm_genes) %>%
    select(gene, avg_log2FC) %>%
    distinct(gene, .keep_all = TRUE)
  
  # complete all genes; absent genes will have NA avg_log2FC
  data.frame(gene = dm_genes, stringsAsFactors = FALSE) %>%
    left_join(sub_df, by = "gene") %>%
    mutate(row_name = row_lab)
  
}) %>%
  bind_rows()

# order rows/columns
plot_df$row_name <- factor(plot_df$row_name, levels = rev(cluster_map$row_name))
plot_df$gene <- factor(plot_df$gene, levels = dm_genes)

# -----------------------------
# HEATMAP
# -----------------------------
ggplot(plot_df, aes(x = gene, y = row_name, fill = avg_log2FC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey",
    name = "avg_log2FC",
    
  ) +
  theme_bw() +
  labs(x = "", y = "") +
  theme(
    panel.grid = element_blank(),
    # tick labels
    axis.text.x = element_text(size = 30, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30)
    #axis.text.x = element_text(angle = 45, hjust = 1)
  )


###################################################

s17 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage17/seurat_output/stage17-res22.rds")
s19 <- readRDS("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/stage19/seurat_output/stage19-res23.rds")


Idents(s19) = s19$RNA_snn_res.2.3
table(Idents(s19))

dm_genes =  c("sox2","sox3", "nes",  "epcam", "cdh1",
           "snai1", "snai2", "foxd3", "twist1", "sox9", "sox10", "tfap2a", "tfap2b",
           "adam33"," mmp14", "pax3", "msx1", "zic2", "zic3")  

dm_genes1 =  c("tfap2a", "tfap2b","sox9", "sox10", "snai1", "snai2",
               "twist1", "pax3", "msx1","zic1", "zic2", "zic3") 



pn_genes = c("isl1", "isl2", "prph","calca", "tac1", "trpv1", "ret", "ngfr", "ntrk1", "ntrk2", "ntrk3", 
             "nefl", "nefm", "nefh",  "mbp", "pvalb", "slc17a7", "chat", "onecut1", "onecut2", "mnx1")


wnt_genes = c("wnt1", "wnt3a", "wnt8a", "wnt8b", "wnt2b", "wnt4", "wnt5b", "wnt7a", "wnt7c", 
              "wnt9a", "wnt9b", "wnt10b", "wnt11", "wnt11b" )

epi_genes = c("mcidas" ,  "myb", "foxj1",
              "ccno", "deup1",  "cep152",
              "tuba1a",  "tekt2",  "dynlrb2", "dnah5", "dnah9"
)
epi2 =  c( "foxj1", "deup1")

  
p = DotPlot(s24, features = epi2,
            #idents = c(3,9,10,20),
            dot.scale = 5)
p = p + 
  theme(
    axis.text.y = element_text(size = 15,angle = 45 ),
    axis.text.x = element_text(size = 20,angle = 45, hjust = 1 )
    #axis.title.y = element_text(size = 10)
  ) + ggtitle("Stage 24: ciliated epidermal progenitor markers")

############
table(Idents(s17))
