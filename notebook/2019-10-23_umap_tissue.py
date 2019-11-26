#%%
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP

import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use(["seaborn-poster", "seaborn-white"])
#%%
import os

try:
    os.chdir(os.path.join(os.getcwd(), "notebook"))
    print(os.getcwd())
except:
    pass

#%%
def tpm(df, gene_length, scale_library=1e6, scale_length=1e3, log=None):
    """Transcripts Per Killobase Million.

    Calcualtes TPM which normalizes by library size and gene length.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame with genes as rows and samples as columns.
    gene_length : pd.Series
        Series were index matches df.index and values are gene lengths.
    scale_library : int or float
        Scaling factor to scale the library size.
    scale_length : int or float
        Scaling factor to scale gene model lengths.
    log : None or function or str
        If a function is giving will apply this to the data before returning.
        If a string is given then it uses numpy log functions that correspond
        to the provided string.

    """

    rpk = (df.T / (gene_length / scale_length)).T
    totals = rpk.sum()

    if log is None:
        log = lambda x: x - 1
    elif log == "log2":
        log = np.log2
    elif log == "log10":
        log = np.log10
    elif log == "ln":
        log = np.log

    return log((rpk / (totals / scale_library)) + 1)


def zscore(df):
    mu = df.mean(axis=1)
    sigma = df.std(axis=1)
    return df.subtract(mu, axis="rows").div(sigma, axis="rows")

#%%
gene_lengths = pd.read_feather("../../larval_gonad/references/gene_annotation_dmel_r6-26.feather", columns=["FBgn", "length"]).set_index("FBgn").squeeze()

#%%
metadata = (
    pd.read_csv(
        "../output/geo-wf/rnaseq_metadata.tsv",
        sep="\t",
        index_col=0,
        usecols=["sample_name", "tissue"],
    )
    .squeeze()
    .dropna()
)

#%%
counts = pd.read_parquet("../output/aln-downstream-wf/gene_counts_wide_tpm.parquet", columns=metadata.index)

#%%
counts = counts.reindex(columns=metadata.sort_values().index)

#%%
norm = tpm(counts, gene_lengths).fillna(0)


#%%

sns.clustermap(
    norm.iloc[:, :1000],
    xticklabels=False,
    yticklabels=False,
    cmap="RdBu_r",
    rasterized=True,
    col_cluster=False,
    figsize=(8, 8)
)




















#%%
counts.head()

#%%
umap = UMAP(n_components=2, random_state=42,)
umap_embeded = umap.fit_transform(counts.T).T

#%%
tsne = TSNE(n_components=2, init='pca', random_state=42)
tsne_embeded = tsne.fit_transform(counts.T).T

#%%
print(umap_embeded.shape)

#%%
metadata = metadata.reindex(counts.columns)

#%%
df = pd.DataFrame(tsne_embeded.T, index=counts.columns, columns=['tSNE 1', 'tSNE 2']).join(metadata)

#%%
fig, ax = plt.subplots(figsize=(20, 20))

for (tissue, dd), color in zip(df.sort_values("tissue").groupby("tissue"), sns.color_palette("husl", n_colors=df.tissue.unique().shape[0])):
    plt.scatter(dd['tSNE 1'], dd['tSNE 2'], c=[color], label=tissue, s=20, linewidths=.1, edgecolors='k')

x, y = df.query("tissue == 'fat body'").mean()
ax.text(x, y, "Fat Body")

#%%

fig, ax = plt.subplots(figsize=(20, 20))

defaults = dict(s=20, linewidths=.1, edgecolors='k')
plt.scatter(df['tSNE 1'], df['tSNE 2'], c='lightgrey', **defaults)
defaults['s'] = 80
plt.scatter(df.query("tissue == 'ovary'")['tSNE 1'], df.query("tissue == 'ovary'")['tSNE 2'], c='red', label='ovary', **defaults)
plt.scatter(df.query("tissue == 'testis'")['tSNE 1'], df.query("tissue == 'testis'")['tSNE 2'], c='blue', label='testis', **defaults)
plt.scatter(df.query("tissue == 'head'")['tSNE 1'], df.query("tissue == 'head'")['tSNE 2'], c='green', label='head', **defaults)
plt.scatter(df.query("tissue == 'fat body'")['tSNE 1'], df.query("tissue == 'fat body'")['tSNE 2'], c='black', label='fat body', **defaults)

plt.legend()