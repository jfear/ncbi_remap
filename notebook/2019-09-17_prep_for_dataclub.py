#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pymongo import MongoClient
from sklearn.decomposition.pca import PCA

plt.style.use("seaborn-talk")
#%%
# Connect to databases
client = MongoClient()
db = client["sramongo"]
ncbi = db["ncbi"]

store = pd.HDFStore("output/sra.h5")

#%%[markdown]
# # Peak at new samles
#%%
new_srxs = store["prealn/queue"].srx.tolist()
new_samples = (
    pd.DataFrame(
        list(ncbi.find({"srx": {"$in": new_srxs}}, {"_id": 0, "srx": 1, "library_strategy": 1}))
    )
    .set_index("srx")
    .squeeze()
)

new_samples.value_counts().to_frame()

#%%[markdown]
# # Summary what is in the SRA
#%%
# library strategy
top_ten_data_types = (
    pd.read_parquet("output/metadata-wf/select_library_strategy.parquet")
    .squeeze()
    .value_counts()
    .head(10)
    .rename("count")
    .to_frame()
    .rename_axis("strategy")
    .reset_index()
)

ax = sns.barplot("strategy", "count", data=top_ten_data_types)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
ax.set(title="Library Strategy", ylabel="Number of Samples", xlabel="")

#%%
# Sample type

#%%
def _value_counts(df, column):
    _df = biometa[column].fillna("not annotated").value_counts().reset_index()
    _df.loc[_df[column] < 100, "index"] = "Other"
    return _df.groupby("index").sum().sort_values(column, ascending=False).squeeze()


def annot_func(pct, allvals):
    absolute = int(pct / 100.0 * np.sum(allvals))
    if pct > 3:
        return f"{pct:.1f}\n({absolute:,})".format(pct, absolute)
    return ""


def make_pie(df, ax):
    ax.pie(df, labels=df.index, autopct=lambda pct: func(pct, df), wedgeprops=dict(width=0.5))


biometa = pd.read_csv(
    "output/geo-wf/rnaseq_metadata.tsv",
    sep="\t",
    index_col=0,
    usecols=["sample_name", "sex", "developmental stage", "tissue", "cell type"],
)

_, axes = plt.subplots(2, 2)
make_pie(_value_counts(biometa, "sex"), axes.flat[0])
make_pie(_value_counts(biometa, "tissue"), axes.flat[1])
make_pie(_value_counts(biometa, "developmental stage"), axes.flat[2])
axes.flat[-1].set_visible(False)

#%%[markdown]
# # Mappability
#%%
rnaseq_srxs = (
    pd.read_parquet("output/metadata-wf/select_library_strategy.parquet")
    .squeeze()
    .pipe(lambda x: x[x == "RNA-Seq"])
    .index.to_list()
)

aln = store["/aln/workflow/hisat2"].query(f"srx == {rnaseq_srxs}")[["num_reads", "per_alignment"]]
dates = pd.DataFrame(
    list(ncbi.find({"srx": {"$in": rnaseq_srxs}}, {"_id": 0, "srx": 1, "sra_create_date": 1}))
).set_index("srx")
df_map = aln.join(dates).sort_values("sra_create_date").set_index("sra_create_date")

_, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw=dict(wspace=0.1))
df_map.num_reads.plot(style="k.", ax=ax1, rasterized=True)
df_map.per_alignment.plot(style="k.", ax=ax2, rasterized=True)
ax1.set(title="Number of Reads", xlabel="Date")
ax2.set(title="Percent Alignment", xlabel="Date")
sns.despine(ax=ax1)
sns.despine(ax=ax2, left=True)

#%%[markdown]
# # Testis

#%%
testis_srx = pd.read_csv("output/geo-wf/rnaseq_metadata.tsv", sep="\t", index_col=0).query("tissue == 'testis'").index.to_list()
gene_counts = pd.read_parquet("output/aln-downstream-wf/gene_counts_wide_tpm.parquet", columns=testis_srx)
corr = gene_counts.corr()

#%%
sns.clustermap(corr, xticklabels=False, yticklabels=False, rasterized=True)

#%%
pca = PCA(n_components=2).fit(gene_counts)
plt.scatter(pca.components_[0, :], pca.components_[1, :], color='k', rasterized=True)
ax = plt.gca()
ax.set(xlabel="PC1", ylabel="PC2")
sns.despine(ax=ax)

#%%
