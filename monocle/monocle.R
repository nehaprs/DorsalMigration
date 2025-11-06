#=====================================
#monocle pseudotime
#created 10.16.2025
#############################
install.packages("R.utils")
remotes::install_github('satijalab/seurat-wrappers')
library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(hdf5r)
library(scCATCH)
library(tidyr)  
library(harmony)
library(SeuratWrappers)
library(monocle3)
library(SeuratObject)

#setwd("~/BINF/scrnaseq general/dorsal migration/full head/merged/noHarmPT")
setwd("~/binf/dorsal migration/dorsal migration")
m <- readRDS("ready_for_monocle.rds")
#DefaultAssay(m)
#Layers(m)


m[["RNA"]] <- JoinLayers(m[["RNA"]])
#Layers(m)
cds = as.cell_data_set(m) #m is the output from merge_preserve)cluster.R

colData(cds)$orig_cluster = m$orig_cluster
colData(cds)$timepoint = m$batch

reducedDims(cds)$UMAP <- Embeddings(m, "umap")
reducedDims(cds)$PCA  <- Embeddings(m, "pca")


cds = cluster_cells(cds, reduction_method = "UMAP")
cds = learn_graph(cds, use_partition = TRUE)
#graph learned using monocle_sparse.R
#name root cells

orig_cluster <- colData(cds)$orig_cluster
root_cells = rownames(colData(cds))[orig_cluster == "s1_10"]
cds <- order_cells(cds, root_cells = root_cells)
# sanity plots
plot_cells(cds, color_cells_by="pseudotime", label_groups_by_cluster=FALSE)
saveRDS(cds,"completedCDS.RDS")

#assign branch lineages of s1_10
library(igraph)
head(g)
# principal graph on UMAP
g <- principal_graph(cds)$UMAP 

class(g) #makes sure g is an igraph object
head(V(g))

comp <- components(g)  




# matrix 'closest' linking nearest graph vertex for each cell
#rows(closest)== cells
#clms(closest)==nearest vertex
#principal_graph_aux: list of auxillary data encapsulating Spatial relationship between cells and the graph
closest <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
closest <- as.matrix(closest[colnames(cds), , drop=FALSE])
head(closest)
# Replace numeric indices of closest with the corresponding vertex names
if (is.numeric(closest[,1])) {
  closest[,1] <- vn[closest[,1]]
}
head(closest)



#Create a named vector mapping each cell barcode to its closest graph vertex
vmap <- setNames(closest[,1], rownames(closest))
head(vmap)



# root vertex:retrieve the nearest vertex (node) for root cells.
root_vs <- vmap[root_cells]
'
table(root_vs)
root_vs
  1   2   3   4   5  10  12  15  16  17  20  22  24  25  27  28  31  36 
  3 267   6   2   1   3  26   8   6  32   2 239  11 131   1   1   4  12

'



#select the most common graph vertex among the root cells
root_v  <- names(sort(table(root_vs), decreasing=TRUE))[1]
#root_v = Y_2


# leaves and partitions
'
definitions:
degree(v)=number of edges incident to v
Interpretation in Monocle3

In the trajectory graph:

Degree = 1 → a leaf node (end of a lineage).

Degree = 2 → an intermediate node along a continuous branch.

Degree ≥ 3 → a branch point (where a lineage splits).

'
deg <- degree(g)
leaf_vs <- names(deg[deg==1])

# precompute shortest paths from root to each leaf
paths <- lapply(leaf_vs, function(L) shortest_paths(g, from=root_v, to=L)$vpath[[1]])
#convert each graph path into a vector of vertex names
#each element of path_set = character vector of node names along one path
path_sets <- lapply(paths, function(p) names(p))



# assign each vertex to the furthest leaf whose path contains it
# vertices shared by multiple leaf paths are "pre-bifurcation"
assign_leaf <- function(v){
  idx <- which(vapply(path_sets, function(S) v %in% S, logical(1)))
  if(length(idx)==0) return(NA_character_)
  if(length(idx)==1) return(leaf_vs[idx])
  return("pre_bifurcation")
}
vertex_to_leaf <- setNames(vapply(V(g)$name, assign_leaf, character(1)), V(g)$name)




# cell lineage label
cell_lineage <- setNames(vertex_to_leaf[vmap], names(vmap))
colData(cds)$lineage_leaf <- cell_lineage
head(colData(cds)$lineage_leaf)



# restrict to descendants reachable from cluster 10 root
descendant_cells <- names(cell_lineage)[!is.na(cell_lineage)]
desc_cds <- cds[, descendant_cells]

#check what label was used for time point
table(cds@colData$batch)
table(cds@colData$orig_cluster)



# summary: where do stage 17 cluster 10 cells go by timepoint?

lin_flow <- as.data.frame(colData(desc_cds)) |>
  mutate(leaf = lineage_leaf,
         tp = batch,
         is_root = (orig_cluster == "s1_10")) |>
  group_by(tp, leaf) |>
  summarise(n=n(), .groups="drop") |>
  arrange(tp, desc(n))
lin_flow
#lineages in lin_flow are Y_1, Y_43, Y_5
write_xlsx(lin_flow,"lin_flow.xlsx")
ggplot(lin_flow, aes(tp, n, fill = leaf)) + geom_bar(stat = "identity", position = "fill")


#through which clusters do they flow into which lineages?”

#install.packages("ggalluvial")
library(ggalluvial)
df <- as.data.frame(colData(desc_cds)) |>
  select(batch, orig_cluster, lineage_leaf) |>
  mutate(n=1) |>
  group_by(batch, orig_cluster, lineage_leaf) |>
  summarise(n=n(), .groups="drop")


ggplot(df,
       aes(axis1=batch, axis2=orig_cluster, axis3=lineage_leaf, y=n)) +
  geom_alluvium() + geom_stratum() + geom_text(stat="stratum", infer.label=TRUE) +
  scale_x_discrete(limits=c("batch","cluster","leaf"), expand=c(.1,.1)) +
  ylab("cells")
table(df)

df1 <- df |>
  dplyr::group_by(orig_cluster) |>
  dplyr::summarise(n = sum(n), .groups = "drop") |>
  dplyr::arrange(desc(n))

df1 = df %>%
  group_by(orig_cluster) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  arrange(desc(n))
  
write_xlsx(df1,"desc_clustersR1.xlsx")

###checking why we see biologically impossible lineages
igraph::components(principal_graph(cds)$UMAP)$no


leaf_tab <- as.data.frame(colData(desc_cds)) |>
  filter(batch== 4) |>
  dplyr::count(lineage_leaf, orig_cluster, name="n") |>
  group_by(lineage_leaf) |>
  mutate(frac = n/sum(n)) |>
  arrange(lineage_leaf, desc(frac))
leaf_tab


names(colData(desc_cds))



end_alloc <- as.data.frame(colData(desc_cds)) |>
  group_by(batch, lineage_leaf) |>
  summarise(cells=n(), .groups="drop")
write.csv(end_alloc, "descendants_endpoint_allocation.csv", row.names=FALSE)



























#plot cells
plot_cells(cds, color_cells_by="lineage_leaf", label_groups_by_cluster=FALSE,
           label_leaves=TRUE, label_branch_points=TRUE)

plot_cells(cds, color_cells_by="orig_cluster", label_groups_by_cluster=FALSE,
           label_leaves=TRUE, label_branch_points=TRUE)


#highlight root cells

# 1. Flag root cells
colData(cds)$is_root <- ifelse(colnames(cds) %in% root_cells, "root", "other")

# 2. Plot with custom colors
plot_cells(
  cds,
  color_cells_by = "is_root",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  show_trajectory_graph = TRUE
) + scale_color_manual(values = c("root" = "red", "other" = "grey"))








