# Load required packages
library(Seurat)
library(monocle3)
library(SeuratWrappers)  # required for as.cell_data_set()

seurat_obj = allCombined
# STEP 1: Convert Seurat object to cds
cds <- as.cell_data_set(seurat_obj)

# STEP 2: Transfer cluster and timepoint metadata
# (If not already present)
colData(cds)$cluster <- Idents(seurat_obj)
colData(cds)$timepoint <- seurat_obj$timepoint  # replace with correct column name

# STEP 3: Preprocess CDS
cds <- preprocess_cds(cds, num_dim = 50)

# STEP 4: Reduce dimension (UMAP)
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# STEP 5: Cluster the cells (needed before learning graph, but ignored later)
cds <- cluster_cells(cds, reduction_method = "UMAP")

# STEP 6: Learn trajectory graph, but constrain to 1 partition
cds <- learn_graph(cds, use_partition = FALSE)

# STEP 7: Select root cells from E9.5 (earliest timepoint)
# Identify cells annotated as E9.5
root_cells <- colnames(cds)[colData(cds)$timepoint == "E9.5"]

# Pick one or more root cells (best to use one near the center of E9.5 cluster)
# If too many, optionally downsample:
# root_cells <- sample(root_cells, 5)

# Order cells in pseudotime from E9.5 root(s)
cds <- order_cells(cds, root_cells = root_cells)

# STEP 8: Visualize pseudotime
plot_cells(
  cds,
  color_cells_by = "timepoint",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_groups_by_cluster = FALSE,
  label_roots = FALSE
)




##########

library(dplyr)
library(monocle3)

# 1. Assign pseudotime, cluster, and timepoint info


pseudotime_df <- as.data.frame(colData(cds))
pseudotime_df$cell_id <- rownames(pseudotime_df)
colnames(pseudotime_df)

ptime = pseudotime(cds)
# Add pseudotime as a new column by matching rownames
pseudotime_df$pseudotime <- ptime[rownames(pseudotime_df)]



# Reorder columns if desired
pseudotime_df <- pseudotime_df %>%
  select(cell_id, pseudotime, cluster, timepoint)
# 2. For each cluster, get its dominant timepoint
cluster_timepoint_map <- pseudotime_df %>%
  group_by(cluster, timepoint) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = n, n = 1) %>%
  ungroup()

# 3. Compute average pseudotime per cluster
cluster_ptime <- pseudotime_df %>%
  group_by(cluster) %>%
  summarise(mean_pseudotime = mean(pseudotime, na.rm = TRUE))

# 4. Combine both
cluster_info <- inner_join(cluster_timepoint_map, cluster_ptime, by = "cluster")









# Compare each E11.5 cluster to all E10.5 clusters before it
e11_clusters <- cluster_info %>% filter(timepoint == "E11.5")
e10_clusters <- cluster_info %>% filter(timepoint == "E10.5")

for (i in seq_len(nrow(e11_clusters))) {
  target <- e11_clusters[i, ]
  cat("\nCluster", target$cluster, "(E11.5):\n  Likely parents from E10.5:\n")
  likely_parents <- e10_clusters %>%
    filter(mean_pseudotime < target$mean_pseudotime) %>%
    arrange(desc(mean_pseudotime))
  print(likely_parents)
}



# Show segment-level transitions
pr_graph <- principal_graph(cds)[["UMAP"]]
edges <- as.data.frame(igraph::as_data_frame(pr_graph))

# You can map cells to nearest graph segment (optional)
# And use igraph functions like igraph::shortest_paths() or distances()





# Get metadata
df <- as.data.frame(colData(cds))
df$cell_id <- rownames(df)
df$pseudotime <- pseudotime(cds)[rownames(df)]

# Get UMAP coordinates
umap_coords <- reducedDims(cds)$UMAP
df$UMAP_1 <- umap_coords[, 1]
df$UMAP_2 <- umap_coords[, 2]

# Check if 'cluster' and 'timepoint' are present
stopifnot(all(c("cluster", "timepoint") %in% colnames(df)))



library(dplyr)

centroids <- df %>%
  group_by(cluster, timepoint) %>%
  summarise(
    x = mean(UMAP_1),
    y = mean(UMAP_2),
    pt = mean(pseudotime, na.rm = TRUE),
    .groups = "drop"
  )






get_parent_links <- function(from_time, to_time) {
  from <- centroids %>% filter(timepoint == from_time)
  to <- centroids %>% filter(timepoint == to_time)
  
  links <- to %>%
    rowwise() %>%
    mutate(
      parent_cluster = from$cluster[which.min(abs(pt - from$pt))],
      parent_x = from$x[which.min(abs(pt - from$pt))],
      parent_y = from$y[which.min(abs(pt - from$pt))]
    ) %>%
    select(parent_cluster, parent_x, parent_y, child_cluster = cluster, child_x = x, child_y = y)
  
  return(links)
}

# Get links from E9.5 → E10.5 and E10.5 → E11.5
links_9to10 <- get_parent_links("E9.5", "E10.5")
links_10to11 <- get_parent_links("E10.5", "E11.5")
all_links <- bind_rows(links_9to10, links_10to11)





library(ggplot2)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = timepoint)) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_segment(data = all_links,
               aes(x = parent_x, y = parent_y, xend = child_x, yend = child_y),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = centroids, aes(x = x, y = y, label = cluster), size = 3, color = "black") +
  theme_minimal() +
  labs(title = "Cluster transitions across developmental timepoints")

############################
#Cluster Transition Network
##############################


library(dplyr)

# Define a clean cluster label: "timepoint_cluster"
centroids <- df %>%
  group_by(cluster, timepoint) %>%
  summarise(
    x = mean(UMAP_1),
    y = mean(UMAP_2),
    pt = mean(pseudotime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(node = paste0(timepoint, "_", cluster))


# Helper function to find nearest earlier cluster
get_edges <- function(from_time, to_time) {
  from <- centroids %>% filter(timepoint == from_time)
  to <- centroids %>% filter(timepoint == to_time)
  
  to %>%
    rowwise() %>%
    mutate(
      parent_idx = which.min(abs(pt - from$pt)),
      from_node = from$node[parent_idx],
      to_node = node
    ) %>%
    select(from_node, to_node)
}

edges_9to10 <- get_edges("E9.5", "E10.5")
edges_10to11 <- get_edges("E10.5", "E11.5")
edges <- bind_rows(edges_9to10, edges_10to11)


# Helper function to find nearest earlier cluster
get_edges <- function(from_time, to_time) {
  from <- centroids %>% filter(timepoint == from_time)
  to <- centroids %>% filter(timepoint == to_time)
  
  to %>%
    rowwise() %>%
    mutate(
      parent_idx = which.min(abs(pt - from$pt)),
      from_node = from$node[parent_idx],
      to_node = node
    ) %>%
    select(from_node, to_node)
}

edges_9to10 <- get_edges("E9.5", "E10.5")
edges_10to11 <- get_edges("E10.5", "E11.5")
edges <- bind_rows(edges_9to10, edges_10to11)



library(igraph)
library(ggraph)
library(ggplot2)

# Create graph
g <- graph_from_data_frame(edges, vertices = centroids %>% select(node, timepoint), directed = TRUE)

# Layout: align by timepoint
time_levels <- c("E9.5", "E10.5", "E11.5")
V(g)$x <- as.numeric(factor(V(g)$timepoint, levels = time_levels))
V(g)$y <- ave(runif(length(V(g))), V(g)$x, FUN = seq_along)  # spread nodes vertically

# Plot
ggraph(g, layout = "manual", x = V(g)$x, y = V(g)$y) +
  geom_edge_link(arrow = arrow(length = unit(4, "mm")), end_cap = circle(3, "mm")) +
  geom_node_circle(aes(r = 0.1), fill = "lightblue") +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  theme_void() +
  ggtitle("Inferred Cluster Transition Network from Pseudotime")

