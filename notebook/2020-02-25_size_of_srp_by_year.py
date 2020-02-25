# %%
import pandas as pd
from pymongo import MongoClient

from ncbi_remap.plotting.sample_submission import get_data

# %%
client = MongoClient()
db = client["sramongo"]
ncbi = db["ncbi"]

# %%
project = (
    pd.DataFrame(ncbi.aggregate([{"$project": {"_id": 0, "srx": 1, "srp": "$study.accn"}}]))
    .set_index("srx")
    .squeeze()
)

# %%
df = (
    pd.read_csv(
        "../output/paper-wf/srx_metadata.tsv",
        sep="\t",
        usecols=["srx", "date_created"],
        index_col="srx",
        parse_dates=["date_created"],
    )
    .assign(year=lambda x: x.date_created.dt.year)
    .join(project)
    .reset_index()
)

# %%
cnts = df.groupby("srp").aggregate({"year": "first", "srx": "size"})

# %%
cnts.sort_values("srx", ascending=False)

# %%
