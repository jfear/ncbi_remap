# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import umap_panel
from ncbi_remap.plotting import sra_umap

plt.style.use(("sra_talk", "sra"))
# %%
embeddings = pd.read_parquet("../output/library_strategy-wf/umap_prealn_features_embeddings.parquet")
metadata = pd.read_parquet("../output/library_strategy-wf/prealn_outliers.parquet")
sra_rnaseq = metadata.library_strategy.str.startswith("RNA-Seq")

# %%
umap_panel.Plot(
    embeddings=embeddings,
    labels=metadata.library_strategy
)

# %%
outliers = (
    pd.read_parquet("../output/library_strategy-wf/prealn_outliers.parquet")
    .pipe(lambda x: (x.library_strategy == "RNA-Seq") & (x.library_strategy_flag_outlier))
)
outliers.sum()

# %%
metadata2 = pd.read_parquet("../output/library_strategy-wf/summarized_metadata.parquet")
final_rnaseq = metadata2.Fear_et_al_library_strategy == "RNA-Seq"
metadata2.loc[final_rnaseq, "library_strategy"] = "RNA-Seq"
metadata2.loc[~final_rnaseq, "library_strategy"] = "Other"


# %%
sra_umap.Plot(
    embeddings=embeddings,
    labels=metadata2.library_strategy,
    outlier_flag=outliers,
    plot_kwargs=dict(s=20)

)

# %%
