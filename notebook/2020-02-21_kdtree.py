# %%
import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.normalization import tpm

# %%
gene_lengths = pd.read_csv("../output/gene_ts_lengths.tsv", sep="\t", index_col=0).squeeze()

# %%
df = tpm(
    pd.read_csv("../output/agg-rnaseq-wf/gene_counts.tsv", sep="\t", index_col=0).T,
    gene_length=gene_lengths
).dropna().T

# %% 
print(df.shape)
df.head()

# %%
df.T.plot(kind="kde")

# %%
highly_var = df.apply(lambda x: x.std() / x.shape[0]).sort_values().tail(2000).index

# %%
binned = df.loc[:, highly_var].head().apply(pd.cut, bins=100, labels=[x for x in range(100)], axis=1)

#%%
fig = plt.figure(figsize=(20, 8))
sns.heatmap(binned)


# %%
tree = BallTree(binned, metric="hamming")
tree.query(binned.values[:1], k=5)

# %%
