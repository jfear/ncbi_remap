# %%
import os
from yaml import full_load
from dateutil.parser import parse as date_parse

import pandas as pd

from pymongo import MongoClient

try:
    os.chdir(os.path.join(os.getcwd(), "notebook"))
    print(os.getcwd())
except:
    pass

pd.options.display.max_rows = 100

# %%
config = full_load(open("../config/common.yaml"))
SRA_TERMS = full_load(open("../config/sra_terms.yaml"))

# %%
# Get list of library strategy in SRA.
client = MongoClient()
db = client["sramongo"]
ncbi = db["ncbi"]

sra_strategies = pd.DataFrame(
    ncbi.aggregate(
        [
            {
                "$match": {
                    "sra_create_date": {"$lte": date_parse(config["download_date"])}
                }
            },
            {
                "$project": {
                    "_id": False,
                    "samples": "$srx",
                    "studies": "$study.accn",
                    "library_strategy": True,
                }
            }
        ]
    )
)

client.close()

# %%
# Number of Samples/Studies by SRA Library Strategy
sra_strategies.replace(dict(other="OTHER")).groupby("library_strategy").agg("nunique").drop(
    "library_strategy", axis=1
).reindex(SRA_TERMS["library_strategy"].keys()).sort_values("samples", ascending=False)


# %%
# Split up free text term(s)
df = pd.read_parquet("../output/library_strategy-wf/free_text_library_strategy.parquet").squeeze()
res = []
for srx, value in df.items():
    if value is None:
        continue
    for v in value.split("|"):
        if v in SRA_TERMS["library_strategy"].keys():
            continue
        res.append((srx, v))

free = (
    pd.DataFrame(res, columns=["srx", "library_strategy"])
    .set_index("srx")
    .squeeze()
    .rename_axis("samples")
)

# %%
pd.merge(sra_strategies[["samples", "studies"]], free, on="samples").groupby(
    "library_strategy"
).agg("nunique").drop("library_strategy", axis=1).sort_values("samples", ascending=False)


# %%
[
    term
    for term in full_load(open("../src/lookups/library_strategy.yaml"))["terms"]
    if term in SRA_TERMS["library_strategy"].keys()
]

# %%
