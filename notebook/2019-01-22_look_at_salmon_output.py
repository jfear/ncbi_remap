# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 0.8.6
#   kernelspec:
#     display_name: Python [conda env:ncbi_remap]
#     language: python
#     name: conda-env-ncbi_remap-py
# ---

# %% [markdown]
# # modENCODE Salmon Output

# %% [markdown]
# Chris and I are trying to figure out optimized Salmon parameters. I went ahead and ran Salmon the 703 modENCODE RNA-seq like samples:
#
# ```                
# # PE Libraries
# salmon quant --index {index_dir} --output {outdir} --threads {threads} --libType=A --gcBias --seqBias -1 {fastqs[0]} -2 {fastqs[1]}
#
# # SE Libraries
# salmon quant --index {index_dir} --output {outdir} --threads {threads} --libType=A --gcBias --seqBias -r {input.fastq}
# ```
#
# There were 13 samples that would not complete the workflow. 
#
# Here I look at how well the samples correlate using the isoforms profiles.

# %%
import os
import sys
from pathlib import Path

from IPython.display import display, HTML, Markdown
import numpy as np
from scipy.stats import spearmanr
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

from dask.delayed import delayed
from dask.distributed import Client

# Project level imports
sys.path.insert(0, '../lib')
from ncbi_remap.notebook import Nb
from ncbi_remap.plotting import make_figs

# Setup notebook
nbconfig = Nb.setup_notebook()

# Connect to data store
store = pd.HDFStore('../output/sra.h5', mode='r')

# %% [markdown]
# ## Read in all the needed data

# %%
# Set-up minicluster
client = Client()
client

# %% [markdown]
# ### modENCODE data
#
# Here I am using my modENCODE sample table and metadata I have generated as part of the SRA project. Note I only have metadata generated for 239 samples.

# %%
# get list of modENCODE srxs from RNA-Seq experiments
mod = pd.read_csv('../output/modENCODE_sampletable.tsv', sep='\t')
modENCODE_srxs = mod[mod.modENCODE_type.str.contains('RNA-seq')].srx.unique().tolist()

# metadata
metadata = pd.read_csv('../output/geo-wf/rnaseq_metadata.tsv', sep='\t', index_col=0)
metadata = metadata.loc[metadata.index.isin(modENCODE_srxs), ['sex', 'developmental stage', 'tissue', 'cell type']].copy()

# %% [markdown]
# ### Read in Salmon Quantifications
#
# Read in each of the salmon results files keeping the TPM value.

# %%
# Read in salmon counts keeping only TPM information, use minicluster to do this is parallel
@delayed
def read_salmon(srx):
    try:
        fname = f'../output/rnaseq-analysis-wf/salmon_counts/{srx}.salmon/quant.sf'
        return pd.read_csv(fname, sep='\t', names=['FBtr', 'Length', 'EffectiveLength', srx, 'NumReads'], header=0, index_col=0)[srx]
    except FileNotFoundError:
        return

work = [read_salmon(srx) for srx in modENCODE_srxs]
futures = client.compute(work)
data = client.gather(futures)
df = pd.concat([sr for sr in data if sr is not None], sort=True, axis=1)

# %% [markdown]
# ## Munge Data for Clustering

# %%
# scale data using mean std
X = df.T
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# %%
# Fit PCA to figure out how many components to use
pca = PCA()
pca.fit(X_scaled)
y = pca.explained_variance_ratio_

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=plt.figaspect(1/2), sharey=True)
ax1.scatter(range(1, len(y) + 1), y, s=2)
ax1.set_ylabel('Proportion Explained Variance')
ax1.set_xlabel('PCs')
ax1.set_title('All PCs')

ax2.scatter(range(1, len(y) + 1), y, s=2)
ax2.set_xlabel('PCs')
ax2.set_xlim(None, 100)
ax2.set_title('All PCs (Zoomed)')

# %%
# Reduce dimensions for 40 PCs
pca = PCA(n_components=40)
pca.fit(X_scaled)
X_pca = pca.transform(X_scaled)

# %%
plt.scatter(X_pca[:, 0], X_pca[:, 1])
plt.xlabel('PC1: Variance Explained ({})'.format(np.round(pca.explained_variance_ratio_[0], 4)))
plt.ylabel('PC2: Variance Explained ({})'.format(np.round(pca.explained_variance_ratio_[1], 4)));

# %% [markdown]
# ## KMeans clustering

# %% [markdown]
# Run a quick KMeans clustering with 15 clusters. Also plot a correlation heatmap ordered by hierarchical clustering. This will give a quick look at how well samples group together using two clustering methods. I also plot the sample in tSNE space to see how well they separate.

# %%
# Run KMeans with 15 clusters
kmeans = KMeans(n_clusters=15).fit(X_pca)
labels = kmeans.labels_

# Make colors for the KMeans clusters
cmap = {k: v for k, v in enumerate(sns.color_palette(n_colors=15))}
colors = pd.Series(list(map(lambda x: cmap[x], labels)), index=df.columns, name='KMeans')

# %%
# Plot correlation heatmap and add row colors for KMeans
corr = pd.DataFrame(spearmanr(X_pca, axis=1)[0], index=df.columns, columns=df.columns)
sns.clustermap(corr, row_colors=colors, xticklabels=False, yticklabels=False);

# %%
# Plot tSNE colored by KMeans
from sklearn.manifold import TSNE
tsne = TSNE()
X_tsne = tsne.fit_transform(X_pca)
df_tsne = pd.DataFrame(X_tsne, columns=['tsne1', 'tsne2'], index=df.columns)
df_tsne['KMean_labels'] = labels
sns.lmplot('tsne1', 'tsne2', hue='KMean_labels', data=df_tsne, fit_reg=False, palette=sns.color_palette(n_colors=15));

# %% [markdown]
# There are some outliers that behave strange, I am guessin these may be the small-RNAs.

# %%
assay = mod[['srx', 'modENCODE_assay']].drop_duplicates('srx').set_index('srx')
assay = assay.reindex(df_tsne.index).dropna().copy()
df_tsne['flag_small_rna'] = (assay.modENCODE_assay == 'small-RNA').values
sns.lmplot('tsne1', 'tsne2', hue='flag_small_rna', data=df_tsne, fit_reg=False, palette=sns.color_palette(n_colors=2));

# %% [markdown]
# ## Compare clustering with biological metadata

# %% [markdown]
# Using my metadata for 239 sample, look and see how well they cluster.

# %%
# Pull out the samples I have metadata for
metadata.fillna('None', inplace=True)
corr_short = corr.reindex(metadata.index).reindex(metadata.index, axis=1)

# %%
# Create colormaps for different attributes
sex_cmap = dict(zip(metadata.sex.unique(), sns.color_palette(n_colors=4)))
dev_cmap = dict(zip(metadata['developmental stage'].unique(), sns.color_palette(n_colors=19)))
tissue_cmap = dict(zip(metadata.tissue.unique(), sns.color_palette(n_colors=13)))
cell_cmap = dict(zip(metadata['cell type'].unique(), sns.color_palette(n_colors=10)))

# Make color vectors
row_colors = pd.concat([
    metadata.sex.map(sex_cmap), 
    metadata['developmental stage'].map(dev_cmap), 
    metadata.tissue.map(tissue_cmap), 
    metadata['cell type'].map(cell_cmap), 
], axis=1, sort=True)

# %%
# Plot cluster map
g = sns.clustermap(corr_short, row_colors=row_colors, xticklabels=False, yticklabels=False, figsize=(12, 12))
ax = g.ax_heatmap

markers = [plt.Line2D([], [], color=color, label=label, marker='o', linestyle='') for label, color in sex_cmap.items()]
sex_legend = ax.legend(handles=markers, title='sex', loc='upper left', bbox_to_anchor=[1.2, 1])

markers = [plt.Line2D([], [], color=color, label=label, marker='o', linestyle='') for label, color in tissue_cmap.items()]
tissue_legend = ax.legend(handles=markers, title='tissue', loc='upper left', bbox_to_anchor=[1.4, .6])

markers = [plt.Line2D([], [], color=color, label=label, marker='o', linestyle='') for label, color in cell_cmap.items()]
cell_legend = ax.legend(handles=markers, title='Cell Type', loc='upper left', bbox_to_anchor=[1, 1])

markers = [plt.Line2D([], [], color=color, label=label, marker='o', linestyle='') for label, color in dev_cmap.items()]
dev_legend = ax.legend(handles=markers, title='Dev Stage', loc='upper left', bbox_to_anchor=[1, .6])

ax.add_artist(sex_legend)
ax.add_artist(tissue_legend)
ax.add_artist(cell_legend)
