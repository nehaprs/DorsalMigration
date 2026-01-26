import scanpy as sc
import scvelo as scv
import numpy as np
from pathlib import Path

scv.settings.verbosity = 3
sc.set_figure_params(dpi=120)

IN_ADATA = Path("/home/neha/dorsalmigration/velocity/velocity/results/adata_with_velocity.h5ad")

OUT_ADATA = Path("/home/neha/dorsalmigration/velocity/velocity/results/adata_velocity_dynamical.h5ad")
FIGDIR = Path("/home/neha/dorsalmigration/velocity/velocity/results/figures")
#FIGDIR.mkdir(parents=True, exist_ok=True)
scv.settings.figdir = str(FIGDIR)

adata = sc.read_h5ad(IN_ADATA)

print("loaded files")

#filter and normalize
scv.pp.filter_and_normalize(
    adata,
    min_shared_counts=20,
    n_top_genes=2000
)

#compute PCA, or use PCA from seurat
if "X_pca" not in adata.obsm:
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")

#compute neighbors on pca space
sc.pp.neighbors(adata,use_rep="X_pca",  n_neighbors=15, n_pcs=40)

#compute moments 
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

print("PCA and moments computed")

#run rna velocity using dynamic model

scv.tl.recover_dynamics(adata, n_jobs = -1)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

#Latent time
scv.tl.latent_time(adata)

#visualization
if "X_umap" not in adata.obsm:
    sc.tl.umap(adata, min_dist=0.3)

#plot velocity stream
scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color=["stage"],
    legend_loc="right",
    save="_stream_stage.png"
)

# Latent time
scv.pl.scatter(
    adata,
    basis="umap",
    color="latent_time",
    color_map="viridis",
    save="_latent_time.png"
)

OUT_ADATA.parent.mkdir(parents=True, exist_ok=True)
adata.write(OUT_ADATA)
print(f"Saved: {OUT_ADATA}")

