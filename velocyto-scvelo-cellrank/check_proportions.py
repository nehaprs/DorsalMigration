#check the proportion of spliced vs unspliced

import scvelo as scv
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
adata = sc.read_h5ad('/home/neha/dorsalmigration/velocity/velocity/results/adata_with_velocity.h5ad')

#scv.pl.proportions(adata)
scv.pl.proportions(adata, save='_proportions.png')
#print(scv.get_df(adata, 'proportions'))
#scv.pl.proportions(adata, groupby='clusters')

