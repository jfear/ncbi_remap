# %% [markdown]
# # RNA-Seq Library Sizes
#
# I am having some issues with running jobs efficiently on Biowulf. I want to
# take a peak at library sizes and see if there is anything I can do to
# simplify my life.

# %%
import pandas as pd
import joblib
from pymongo import MongoClient
from IPython.display import display, Markdown


# %% [markdown]
# ## Comparison of Counts with SRA Values

# %%
# SRA Counts
client = MongoClient()
db = client["sramongo"]["ncbi"]

sra = pd.DataFrame(
    db.aggregate(
        [
            {"$unwind": {"path": "$runs"}},
            {"$match": {"runs.nspots": {"$ne": ""}}},
            {"$project": {"_id": 0, "srx": 1, "srr": "$runs.srr", "libsize": "$runs.nspots"}},
        ]
    )
)
client.close()

# My Counts
libsize = pd.read_parquet("../output/fastq-wf/libsize").libsize.rename("read_counts").astype(int)

# Merged Counts
df = pd.merge(sra, libsize, on="srr")

# %%
num_equal = (df.libsize ==  df.read_counts).sum()
prop_equal = num_equal / df.shape[0]

# %%
display(Markdown(f"Almost all samples ({prop_equal * 100:.2f}%) are the same."))

# %% [markdown]
# ## RNA-Seq Data Sets

# %%
rnaseq = joblib.load("../output/library_strategy-wf/rnaseq_inliers.pkl")
rnaseq_counts = df.query(f"srx == {rnaseq}")
lt_1mil = rnaseq_counts.query("libsize <= 1000000").shape[0]
lt_100k = rnaseq_counts.query("libsize <= 100000").shape[0]

display(rnaseq_counts.libsize.describe().map(lambda x: f"{x:,.0f}").to_frame().T)
print(f"# <= 1,000,000: {lt_1mil:,}")
print(f"# <= 100,000: {lt_100k:,}")
