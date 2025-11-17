library(dplyr)
library(igraph)
library(tidygraph)
library(ggraph)
library(grid)

## 0) Graph + cell→vertex map from desc_cds
g  <- principal_graph(desc_cds)$UMAP
vn <- V(g)$name

head(root_cells)

closest <- desc_cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
closest <- as.matrix(closest[colnames(desc_cds), , drop = FALSE])
vmap <- setNames(if (is.numeric(closest[,1])) vn[closest[,1]] else as.character(closest[,1]),
                 rownames(closest))

## 1) Root at S95_4
root_cells <- colnames(desc_cds)[tolower(colData(desc_cds)$orig_cluster) == "s95_4"]
stopifnot(length(root_cells) > 0)
root_v <- names(sort(table(vmap[root_cells]), decreasing = TRUE))[1]

## 2) Target leaves
#get the list of leaves from lin_flow.png
#
target_leaves <- c("Y_121","Y_134","Y_19","Y_194","Y_195",
                   "Y_25", "Y_257", "Y_262","Y_276","Y_304",
                   "Y_328","Y_334","Y_380","Y_45","Y_515","Y_524",
                   "Y_525", "Y_531", "Y_578","Y_583",
                   "Y_6","Y_93","Y_95") 

#target_leaves <- c("Y_121","Y_194")
target_leaves <- target_leaves[target_leaves %in% vn]
stopifnot(length(target_leaves) > 0)

## 3) Paths root→each leaf (vertex-name sequences)
path_vs <- lapply(target_leaves, function(L)
  names(shortest_paths(g, from = root_v, to = L)$vpath[[1]])
)
names(path_vs) <- target_leaves

## 4) Trunk (pre-bifurcation): vertices that lie on ≥2 of the selected paths
all_vs <- unlist(path_vs, use.names = FALSE)
pre_bifurcation <- names(which(table(all_vs) >= 2))

## 5) Modal cluster per vertex (for node labels)
modal <- function(x) if (length(x)) names(sort(table(x), decreasing = TRUE))[1] else NA_character_
v2cl <- tapply(as.character(colData(desc_cds)$orig_cluster),
               INDEX = factor(vmap, levels = vn),
               FUN = modal)
v2cl <- unlist(v2cl[!is.na(names(v2cl))], use.names = TRUE)

## 6) Build union subgraph (edges along each path)
edge_df <- do.call(rbind, lapply(names(path_vs), function(leaf){
  pv <- path_vs[[leaf]]
  if (length(pv) < 2) return(NULL)
  data.frame(leaf = leaf, from = pv[-length(pv)], to = pv[-1], stringsAsFactors = FALSE)
})) %>% distinct()

node_df <- data.frame(name = unique(c(edge_df$from, edge_df$to)), stringsAsFactors = FALSE) %>%
  mutate(role = case_when(
    name == root_v ~ "root",
    name %in% pre_bifurcation ~ "pre_bifurcation",
    TRUE ~ "branch"
  ),
  cluster = v2cl[name])

#####################
#find all orig_clusters that can be mapped to these nodes
##########################
library(dplyr)

# node_df already has all vertices used in the lineage tree
# vmap: named vector of cell → vertex
# colData(desc_cds)$orig_cluster: cluster label for each cell

vertex_cluster_tbl <- tibble(
  cell = names(vmap),
  vertex = vmap,
  cluster = as.character(colData(desc_cds)$orig_cluster)
) %>%
  # keep only vertices that are in node_df
  filter(vertex %in% node_df$name) %>%
  group_by(vertex) %>%
  summarise(
    clusters = paste(unique(cluster), collapse = ", "),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(match(vertex, node_df$name))  # optional: keep same order as node_df

vertex_cluster_tbl

# starting from vertex_cluster_tbl created earlier
# Build full mapping first
vertex_cluster_tbl_top5 <- tibble(
  cell = names(vmap),
  vertex = vmap,
  cluster = as.character(colData(desc_cds)$orig_cluster)
) %>%
  filter(vertex %in% node_df$name) %>%
  group_by(vertex, cluster) %>%
  summarise(n_cells = n(), .groups = "drop_last") %>%
  arrange(vertex, desc(n_cells)) %>%
  slice_head(n = 5) %>%                                  # keep top 5 clusters per vertex
  summarise(
    top5_clusters = paste(cluster, collapse = ", "),
    total_cells = sum(n_cells),
    .groups = "drop"
  ) %>%
  arrange(match(vertex, node_df$name))  


#download vertex_cluster_tbl_top5 and manually select the earliest clusters.
#because code is fucked.

#write_xlsx(vertex_cluster_tbl_top5,"vertex_cluster_tbl_top5.xlsx")

node_df1 <- read_excel("vertex_cluster_tbl_top5.xlsx")

role = node_df[,-3]
#node_df2 = inner_join(node_df1, role, by = c("vertex" = "name"))
node_df2 = node_df
################################################

## 7) Map endpoints to indices and build tbl_graph
edges_idx <- edge_df %>%
  transmute(leaf, from = match(from, node_df2$name), to = match(to, node_df2$name)) %>%
  filter(!is.na(from), !is.na(to))

tg <- tbl_graph(nodes = node_df2, edges = edges_idx, directed = TRUE)

## 8) Layout from the root; plot a single tree with branches colored by leaf
lay <- create_layout(tg, layout = "sugiyama")

ggraph(lay) +
  geom_edge_link(aes(colour = leaf),
                 arrow = arrow(type = "open", length = unit(5, "mm")),
                 show.legend = TRUE) +
  geom_node_label(aes(label = ifelse(!is.na(cluster), cluster, name),
                      fill = role),
                  label.size = 0.2, size = 3, show.legend = TRUE) +
  scale_fill_manual(values = c(root = "tomato", pre_bifurcation = "khaki", branch = "lightblue")) +
  labs(colour = "Lineage leaf", fill = "Node role") +
  theme_minimal() +
  theme(panel.grid = element_blank())


ggraph(lay) +
  geom_edge_link(aes(colour = leaf),
                 arrow = arrow(type = "open", length = unit(3, "mm")),
                 show.legend = TRUE) +
  geom_node_label(aes(label = name,
                      fill = role),
                  label.size = 0.2, size = 3, show.legend = TRUE) +
  scale_fill_manual(values = c(root = "tomato", pre_bifurcation = "khaki", branch = "lightblue")) +
  labs(colour = "Lineage leaf", fill = "Node role") +
  theme_minimal() +
  theme(panel.grid = element_blank())

