# %%
import pandas as pd
from joblib import load
from sklearn.neighbors import BallTree


# %%
metadata = pd.read_feather("../../larval_gonad/references/gene_annotation_dmel_r6-26.feather", columns=["FBgn", "gene_symbol"]).set_index("FBgn").squeeze()

# %%
model = load("../output/agg-rnaseq-wf/feature_balltree.pkl") # type: BallTree
reduced = pd.read_parquet("../output/agg-rnaseq-wf/sample_pca_select_components.parquet")
labels = reduced.index.values
# %%
fbgn = metadata[metadata == "dsx"].index[0]
X = reduced.query(f"srx == '{fbgn}'").values

# %%
dist, idx = model.query(X, k=100)
genes = labels[idx[0]]

# %%
res = metadata.reindex(genes).to_frame()
res["Distance"] = dist[0]
res

# %%
pd.options.display.max_rows = 100


# %%
res.to_csv("/home/fearjm/Documents/dsx_neighbors.csv")

# %%
