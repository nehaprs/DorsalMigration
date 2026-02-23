
#barplot

library(dplyr)
library(ggplot2)

setwd("~/BINF/scrnaseq general/dorsal migration/full head/highUMI/integrated")

# Summarize total counts per cluster

'
df <- read_csv("s17-19/sam17-19/S17_18_endpoint_summary.csv")
df <- read_csv("s19-21/S19_10v2par_endpoint_summary.csv")
df <- read_csv("s19-21/S19_17v2par_endpoint_summary.csv")
df <- read_csv("s21-24/S21_10par_endpoint_summary.csv")
df <- read_csv("s21-24/S21_32par_endpoint_summary.csv")
df <- read_csv("s21-24/S21_33par_endpoint_summary.csv")
'
df <- read_csv("s21-24/S21_26par_endpoint_summary.csv")
df <- read_csv("s21-24/S21_2par_endpoint_summary.csv")
df <- read_csv("s21-24/S21_16par_endpoint_summary.csv")
colnames(df)

df_sum = df %>% group_by(endpoint_orig_cluster) %>% summarise(total_counts = sum(count))

# Bar plot

ggplot(df_sum, aes(x = endpoint_orig_cluster, y = total_counts)) +
  geom_bar(stat = "identity") +
  xlab("terminal cluster") +
  ylab("Total hits") +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

##label top clusters only
# Identify top 5 clusters
top5_clusters <- df_sum %>%
  slice_max(total_counts, n = 7, with_ties = FALSE) %>%
  pull(endpoint_orig_cluster)

# Bar plot with labels only for top 5 clusters
ggplot(df_sum, aes(x = endpoint_orig_cluster, y = total_counts)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(
    labels = function(x) ifelse(x %in% top5_clusters, x, "")
  ) +
  xlab("Terminal cluster") +
  ylab("Total hits") +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )