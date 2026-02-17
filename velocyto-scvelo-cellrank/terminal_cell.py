import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv


H5AD_PATH ="/home/neha/dorsalmigration/velocity/velocity/results/sam/adata_scvelo_computed.h5ad" 
CLUSTER_COL = "orig_cluster"
CLUSTER_VALUE = "s17_18"

#simulation controls
N_STEPS=200
N_REPEATS_PER_CELL=200
RANDOM_SEED_BASE=0

#limit the number of start cells if you want
MAX_START_CELLS = None

#output
OUT_ENDPOINTS_SUMMARY_CSV = "/home/neha/dorsalmigration/velocity/velocity/results/sam/S17_18_endpoint_summary.csv"
OUT_PER_START_LONG_CSV = "/home/neha/dorsalmigration/velocity/velocity/results/sam/s17_18_endpoints_per_start_long.csv"
OUT_PER_START_TOPK_CSV = "/home/neha/dorsalmigration/velocity/velocity/results/sam/s17_18_endpoints_per_start_topk.csv"  # per-start topK endpoints
TOPK = 10

adata = sc.read_h5ad(H5AD_PATH)

# Basic checks
if CLUSTER_COL not in adata.obs.columns:
    raise KeyError(f"{CLUSTER_COL=} not found in adata.obs. Available columns: {list(adata.obs.columns)[:20]} ...")

#ensure velocity graph exists

if "velocity_graph" not in adata.uns:
    raise KeyError(
        "adata.uns['velocity_graph'] not found. "
        "You need to run scv.tl.velocity_graph(adata) (or load an object that already has it)."
    )

print("adata loaded, checks done")

#select start cells

start_mask = adata.obs[CLUSTER_COL].astype(str).values == str(CLUSTER_VALUE)
start_indices = np.where(start_mask)[0]


if start_indices.size == 0:
    raise ValueError(f"No cells found with adata.obs[{CLUSTER_COL!r}] == {CLUSTER_VALUE!r}.")


start_barcodes = adata.obs_names[start_indices].to_numpy()

print(f"Loaded:{adata.n_obs} cells, {adata.n_vars} genes")
print(f"Starting cells in {CLUSTER_VALUE}: {len(start_indices)}")

#transition matrix options
BACKWARD = False
N_NEIGHBORS = None
SCALE = 10
SELF_TRANSITIONS=True


#simulate one trajectory and return endpoint index
#If a cell follows RNA velocity for n steps, where does it end up?
def simulate_endpoint_idx(
    adata,
    starting_cell,
    n_steps,
    random_state,
    backward=False,
    n_neighbors=None,
    scale=10,
    #self_transitions=True,
):
    """
    Uses scv.utils.get_cell_transitions.
    If basis is None, returns indices of visited cells. :contentReference[oaicite:6]{index=6}
    Endpoint = last visited index.
    """

    kwargs= dict(n_steps=n_steps,
            backward=backward,
            random_state=random_state,
            scale=scale,
            #self_transitions=self_transitions,
    )
    if n_neighbors is not None:
        kwargs["n_neighbors"] = n_Neighbors



    path = scv.utils.get_cell_transitions(
    adata,
    starting_cell=starting_cell,
    basis=None,
    **kwargs,
    )

    endpoint = int(np.asarray(path)[-1])
    return endpoint

#run simulation

overall_endpoint_counts = {}
rows_long = []

for i, (s_idx, s_bc) in enumerate(zip(start_indices, start_barcodes), start=1):
    endpoints = np.empty(N_REPEATS_PER_CELL, dtype = int)

    for r in range(N_REPEATS_PER_CELL):
        seed = RANDOM_SEED_BASE + (i* 1_000_000) + r
        endpoints[r] = simulate_endpoint_idx(
            adata=adata,
            starting_cell= int(s_idx),
            n_steps = N_STEPS,
            random_state=int(seed),
            backward=BACKWARD,
            n_neighbors=N_NEIGHBORS,
            scale=SCALE,
            #self_transitions=SELF_TRANSITIONS,
            )
        

    # Update overall counts
    ep_unique, ep_counts = np.unique(endpoints, return_counts=True)
    for ep, c in zip(ep_unique, ep_counts):
        overall_endpoint_counts[int(ep)] = overall_endpoint_counts.get(int(ep), 0) + int(c)

    # Store per-start distribution (long format)
    total = float(N_REPEATS_PER_CELL)
    for ep, c in zip(ep_unique, ep_counts):
        rows_long.append(
            {
                "start_index": int(s_idx),
                "start_barcode": str(s_bc),
                "endpoint_index": int(ep),
                "endpoint_barcode": str(adata.obs_names[int(ep)]),
                "count": int(c),
                "prob": float(c) / total,
            }
        )

    if i % 25 == 0 or i == len(start_indices):
        print(f"Processed {i}/{len(start_indices)} start cells")


# Overall endpoint summary
overall = (
    pd.DataFrame(
        [{"endpoint_index": k, "count": v} for k, v in overall_endpoint_counts.items()]
    )
    .sort_values("count", ascending=False)
    .reset_index(drop=True)
)

overall["endpoint_barcode"] = overall["endpoint_index"].map(lambda ix: str(adata.obs_names[int(ix)]))
overall["prob_overall"] = overall["count"] / overall["count"].sum()



for col in [CLUSTER_COL]:
    if col in adata.obs.columns:
        overall[f"endpoint_{col}"] = overall["endpoint_index"].map(
            lambda ix: str(adata.obs.iloc[int(ix)][col])
        )

# Per-start long table
df_long = pd.DataFrame(rows_long).sort_values(
    ["start_barcode", "prob"], ascending=[True, False]
)

# Per-start topK
df_topk = (
    df_long.groupby("start_barcode", sort=False, as_index=False)
    .apply(lambda g: g.nlargest(TOPK, "prob"))
    .reset_index(drop=True)
)


for col in [CLUSTER_COL]:
    if col in adata.obs.columns:
        df_long[f"endpoint_{col}"] = df_long["endpoint_index"].map(
            lambda ix: str(adata.obs.iloc[int(ix)][col])
        )
        df_topk[f"endpoint_{col}"] = df_topk["endpoint_index"].map(
            lambda ix: str(adata.obs.iloc[int(ix)][col])
        )


# -------------------------
# SAVE
# -------------------------
overall.to_csv(OUT_ENDPOINTS_SUMMARY_CSV, index=False)
df_long.to_csv(OUT_PER_START_LONG_CSV, index=False)
df_topk.to_csv(OUT_PER_START_TOPK_CSV, index=False)

print("\nSaved:")
print(f"  Overall endpoint summary: {OUT_ENDPOINTS_SUMMARY_CSV}")
print(f"  Per-start endpoints (long): {OUT_PER_START_LONG_CSV}")
print(f"  Per-start endpoints (top{TOPK}): {OUT_PER_START_TOPK_CSV}")


















































    
















