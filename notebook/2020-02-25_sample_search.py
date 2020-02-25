# %%
import pandas as pd
from joblib import load
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import BallTree


# %%
metadata = pd.read_csv("../output/geo-wf/rnaseq_metadata.tsv", sep="\t", index_col=0, usecols=["sample_name", "title", "study", "tissue", "description"])

# %%
model = load("../output/agg-rnaseq-wf/sample_balltree.pkl") # type: BallTree
reduced = pd.read_parquet("../output/agg-rnaseq-wf/feature_pca_select_components.parquet")
labels = reduced.index.values
# %%
X = reduced.query("srx == 'SRX2877981'").values
dist, idx = model.query(X, k=10)
srxs = labels[idx[0]]

# %%
res = metadata.reindex(srxs)
res["Distance"] = dist[0]
res

# %%
