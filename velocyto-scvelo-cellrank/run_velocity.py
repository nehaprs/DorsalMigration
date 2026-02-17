import scanpy as sc
import scvelo as scv
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
import re

scv.settings.verbosity = 3
scv.settings.presenter_view = True
sc.set_figure_params(dpi=120)

ROOT = Path(__file__).resolve().parent

H5AD = ROOT / "s17_19.h5ad"
LOOM_S1 = ROOT.parent /"stage17" / "data" / "cellranger" / "velocyto" \
         /"stage17.loom"
LOOM_S2 = ROOT.parent /"stage19" / "data" / "cellranger" / "velocyto" \
          / "stage19.loom"

OUT_ADATA = ROOT / "results" / "adatasam_with_velocity.h5ad"

ROOT.joinpath("results").mkdir(parents=True, exist_ok=True)

#load anndata
adata = sc.read_h5ad(H5AD)

print("loaded anndata")

#load loom files
ldata1 = sc.read_loom(LOOM_S1)
ldata2 = sc.read_loom(LOOM_S2)
print("loaded looms")


#make functions to harmonize cell barcodes


def std_barcode_integrated(index):
    out = []
    for x in index.astype(str):
        # keep stage/sample prefix
        prefix, bc = x.split("_", 1)
        bc = re.sub(r"-\d+$", "", bc)   # remove -1
        out.append(f"{prefix}_{bc}")
    return pd.Index(out)

def std_barcode_loom(index):
    out = []
    for x in index.astype(str):
        x = x.replace("cellranger:", "")
        bc = re.sub(r"[xX]$", "", x)
        
        out.append(bc)
    return pd.Index(out)




"""
adata_join = standardize_barcode(adata.obs_names)
ldata1.obs_names = standardize_barcode(ldata1.obs_names)
ldata2.obs_names = standardize_barcode(ldata2.obs_names)
"""
# Tag stage in loom objects 
ldata1.obs["stage"] = "st17"
ldata2.obs["stage"] = "st19"

# prepend stage after 'cellranger:' in the barcode
ldata1.obs_names = pd.Index(
    [bc.replace("cellranger:", "cellranger:s17_") for bc in ldata1.obs_names]
)
ldata2.obs_names = pd.Index(
    [bc.replace("cellranger:", "cellranger:s19_") for bc in ldata2.obs_names]
)


print("ldata1 var_names unique:", ldata1.var_names.is_unique)
print("ldata2 var_names unique:", ldata2.var_names.is_unique)

ldata1.var_names_make_unique()
ldata2.var_names_make_unique()

#concatenate loom objects
ldata = ad.concat(
    [ldata1, ldata2],
    join="outer",
    merge="first",
    label="loom_batch",
    fill_value=0
)

print("ldata concatenated")

#standardize names


adata.obs_names = std_barcode_integrated(adata.obs_names)
ldata.obs_names = std_barcode_loom(ldata.obs_names)


print("adata barcodes:")
print(adata.obs_names[:10])

print("\nldata barcodes:")
print(ldata.obs_names[:10])

common = adata.obs_names.intersection(ldata.obs_names)
print("common:", len(common), "of", adata.n_obs)

#subset to common cells, and ensure same order
adata = adata[common].copy()
ldata = ldata[common].copy()

#verfiy
print((adata.obs_names == ldata.obs_names).all()) #should return True
print("ldata layers:", ldata.layers.keys())
#bring spliced/unspliced layer to integrated data
'''
for layer in ["spliced", "unspliced", "ambiguous"]:
    if layer in ldata.layers:
        adata.layers[layer] = ldata.layers[layer]
    else:
        print(f"Layer {layer} not found in loom; continuing.")

print("spliced and unspliced brought in")
print(adata)
print("Layers:", list(adata.layers.keys()))
'''
adatasam = scv.utils.merge(adata, ldata)
adatasam.write(OUT_ADATA)
print(f"Saved: {OUT_ADATA}")








''' 
#intersect cell from integrated AnnData and loom-derived AnnData
common = adata_join.intersection(ldata.obs_names)
if len(common) < 0.8 * adata.n_obs:
    print(f"Warning: only {len(common)} / {adata.n_obs} integrated cells matched loom cells.")
    print("This usually means barcode naming mismatch. Fix standardize_barcode().")


print("common cells checked")
'''
