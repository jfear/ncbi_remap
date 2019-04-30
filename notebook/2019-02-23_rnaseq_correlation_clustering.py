# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.0.0
#   kernelspec:
#     display_name: Python [conda env:ncbi_remap]
#     language: python
#     name: conda-env-ncbi_remap-py
# ---

# %%
from pickle import dump

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist


# %%
X = pd.read_parquet('/home/fearjm/scratch/ncbi_remap/rnaseq_correlations.parquet').values

# %%
d = pdist(X)
del X

# %%
L = linkage(d, method='average')

# %%
tree = dendrogram(L, no_plot=True)

# %%
with open('../output/notebook/2019-02-23_rnaseq_correlation_tree.pkl', 'wb') as fh:
    dump(tree, fh)

# %%
leaves = tree['leaves']

# %%
df = pd.read_parquet('/home/fearjm/scratch/ncbi_remap/rnaseq_correlations.parquet')

# %%
idx = df.index[leaves]

# %%
df = df.reindex(index=idx, columns=idx)

# %%
df.to_parquet('/home/fearjm/scratch/ncbi_remap/rnaseq_correlations_clustered.parquet')

# %%
