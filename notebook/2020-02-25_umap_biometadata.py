# %%
from pickle import load

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting.sra_umap import Plot
from ncbi_remap.plotting import umap_panel

# %%
plt.style.use(("sra_talk", "sra"))


# %%
metadata = pd.read_csv("../output/geo-wf/rnaseq_metadata.tsv", sep="\t", index_col=0).rename_axis("srx")

# %%
embeddings = pd.read_parquet("../output/agg-rnaseq-wf/umap_embedding_gene_counts.parquet").reindex(metadata.index)


# %%
ovary = metadata.tissue == "ovary"
mask = (metadata.tissue != "ovary") & (metadata.tissue != "testis") & (metadata.tissue != "head")
# %%
metadata["tissue2"] = metadata.tissue
metadata.loc[mask, "tissue2"] = "Other"


# %%
umap_panel.Plot(embeddings=embeddings, labels=metadata.tissue2, plot_kwargs=dict(s=2))

# %%
(metadata.tissue == "testis").sum()

# %%
