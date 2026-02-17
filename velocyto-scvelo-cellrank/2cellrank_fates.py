
import multiprocessing as mp


def main():
    import scanpy as sc
    import scvelo as scv
    import cellrank as cr
    import numpy as np
    import pandas as pd
    from pathlib import Path

    IN_ADATA = Path("/home/neha/dorsalmigration/velocity/velocity/results/sam/adata_scvelo_computed.h5ad")
    OUT_ADATA = Path("/home/neha/dorsalmigration/velocity/velocity/results/sam/adatasam_velocity_cellrank.h5ad")

    FIGDIR = Path("/home/neha/dorsalmigration/velocity/velocity/results/sam")

    adata = sc.read_h5ad(IN_ADATA)
    print("loaded files")

    #choose cluster annotation columns from seurat
    cluster_key = "orig_cluster"

    #build cellrank kernel from velocity
    vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix(n_jobs=1, 
    backend="threading",show_progress_bar=False,)

    print("velocity kernel built")

    #combine with connectivity kernel

    ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
    kernel = 0.6 * vk + 0.4 * ck

    print("velocity and connectivity kernels combined")

    #estimate macrostate

    g = cr.estimators.GPCCA(kernel)
    g.compute_schur(n_components=20)
    g.compute_macrostates(n_states=12)

   #init_mask = adata.obs["stage"].isin(["s17"]).to_numpy()
   # init_ix = adata.obs_names[init_mask].tolist()

    g.predict_terminal_states()
    g.predict_initial_states(allow_overlap=True)
    
   

    print("computed initial and terminal states")

    #fate probabilities of each cell towards terminal state
    g.compute_fate_probabilities()
    print("fate probabilities computed")


    #save results 
    adata = g.adata

    #visualize terminal states and fate probabilities
    #terminal states over umap
    g.plot_macrostates(which="terminal",basis="umap",save=str(FIGDIR / "terminal_states.png"))

    #fate probabilities in embedding
    g.plot_fate_probabilities(basis="umap", save=str(FIGDIR / "fate_probabilities_umap.png"))

    print("images saved")

    # identifier for stages

    stage1_mask = (adata.obs["stage"] == "s17")
    stage2_mask = (adata.obs["stage"] == "s19")

    #define origin cells of interest: DM cells of the earlier cluster
    stage1_interest = stage1_mask & (adata.obs["orig_cluster"] == "s17_18")

    print("defined cells of interest")

    fates = adata.obsm["lineages_fwd"] #fate probailities are stored here
    fate_names = list(fates.names)         


    fate_df = pd.DataFrame(
        np.asarray(fates),                 
        index=adata.obs_names,
        columns=fate_names,
    )

    print("fates identified")

    summary = pd.DataFrame({
        "terminal_state": fate_names,
        "mean_fate_prob_stage1_interest": np.asarray(
            fate_df.loc[stage1_interest].mean(axis=0)
        ).ravel(),
        "mean_fate_prob_all_stage1": np.asarray(
            fate_df.loc[stage1_mask].mean(axis=0)
        ).ravel()
    })

    summary = summary.sort_values("mean_fate_prob_stage1_interest", ascending=False)

    print("created summary dataframe of fates")

    summary.to_csv("/home/neha/dorsalmigration/velocity/velocity/results/sam/stage17_dm_fate_summary.csv",
                index=False)

    #compute cluster-level fate mapping
    clust_summary = (
        fate_df.loc[stage1_mask]
        .groupby(adata.obs.loc[stage1_mask, "orig_cluster"].astype(str))
        .mean()
    )


    print(clust_summary.shape)
    print(clust_summary.index)
    print(clust_summary.columns)

    clust_summary.to_csv("/home/neha/dorsalmigration/velocity/velocity/results/sam/stage17_cluster_to_terminalstates.csv")

    adata.write(OUT_ADATA)
    print("Wrote:", OUT_ADATA)
    print("Top predicted terminal states for stage1 interest:")
    print(summary.head(10))


if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)
    main()
































































