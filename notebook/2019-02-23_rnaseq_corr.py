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
import numpy as np
import pandas as pd

# %%
X = pd.read_parquet('/home/fearjm/scratch/ncbi_remap/rnaseq_wide_ranks.parquet').T.values

# %%
X.shape

# %%
corr = np.corrcoef(X)

# %%
index = pd.read_parquet('/home/fearjm/scratch/ncbi_remap/rnaseq_wide_ranks.parquet').columns

# %%
corr_df = pd.DataFrame(corr, index=index, columns=index)

# %%
corr_df.to_parquet('/home/fearjm/scratch/ncbi_remap/rnaseq_correlations.parquet')

# %%
