import scanpy as sc
import scvelo as scv
import numpy as np

adata = sc.read_h5ad('/home/neha/dorsalmigration/velocity/velocity/results/adatasam_with_velocity.h5ad')

scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata)

scv.tl.velocity_graph(adata, n_jobs=1)
scv.set_figure_params()

X = adata.obsm["UMAP"]
if hasattr(X, "loc"):  # pandas DataFrame
    X = X.loc[adata.obs_names].to_numpy()

adata.obsm["X_umap"] = np.asarray(X)

scv.pl.velocity_embedding(adata, basis="umap", color = "orig_cluster",
 save="/home/neha/dorsalmigration/velocity/velocity/results/sam/stochastic_velocity.png")

scv.pl.velocity_embedding_stream(adata, basis="umap", color = "orig_cluster",
    legend_loc="none",
    size=20,
    alpha=0.6, save="/home/neha/dorsalmigration/velocity/velocity/results/sam/stochastic_velocity_stream.png")

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c = keys, perc = [5,95], save = "/home/neha/dorsalmigration/velocity/velocity/results/sam/stochastic_velocity_confidence.png")

scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color = 'velocity_pseudotime', save = "/home/neha/dorsalmigration/velocity/velocity/results/sam/stochastic_velocity_ptime.png")

#dynamic model

scv.tl.recover_dynamics(
    adata,
    n_top_genes=2000,
    n_jobs=-1,
    show_progress_bar=False
)

scv.tl.velocity(adata, mode = 'dynamical')


scv.pl.velocity_embedding_stream(adata, basis = 'umap', color = "orig_cluster",
    legend_loc="none",  save = "/home/neha/dorsalmigration/velocity/velocity/results/sam/dynamic_velocity_stream.png")


scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80,  save = "/home/neha/dorsalmigration/velocity/velocity/results/sam/dynamic_velocity_latenttime.png")

adata.write_h5ad("/home/neha/dorsalmigration/velocity/velocity/results/sam/adata_scvelo_computed.h5ad")
