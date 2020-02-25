# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap import plotting

# %%
plt.style.use(("sra_talk", "sra"))
plt.style.use("sra")

# %%
gene_metadata = (
    pd.read_feather(
        "../../larval_gonad/references/gene_annotation_dmel_r6-26.feather",
        columns=["FBgn", "gene_symbol"],
    )
    .set_index("FBgn")
    .squeeze()
)
symbol2fbgn = {v: k for k, v in gene_metadata.items()}

# %%
gene_expression = pd.read_csv("../output/agg-rnaseq-wf/tpm_gene_counts.tsv", sep="\t", index_col=0)


# %%
def zscore(x):
    return (x - x.mean()) / x.std()

# %%
dsx = zscore(gene_expression[symbol2fbgn["dsx"]]).rename("dsx")
pd.cut(dsx, bins=4, labels=["low", "low-mid", "mid-high", "high"]).value_counts().map(lambda x: f"{x:,}").to_frame()

# %%
pd.cut(dsx, bins=4, labels=["low", "low-mid", "mid-high", "high"]).pipe(lambda x: x[x == "high"]).index.tolist()


# %%
ax = sns.kdeplot(dsx)
ax.legend_.remove()
ax.set(ylabel="Density", xlabel="Z-score (TPM)")
sns.despine(ax=ax, left=True, right=True)
ax.set_title("dsx", fontstyle="italic")