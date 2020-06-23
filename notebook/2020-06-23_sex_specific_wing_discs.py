# %% [markdown]
# # Sex Specific Wing Discs
#
# Brian forwarded an email inquery if there was any RNA-Seq data on
# sex-specific wing discs. This notebook is a quick poke to see what is there.

# %%
import sqlite3
from pathlib import Path

import pandas as pd
from pymongo import MongoClient

# %%
# Get biometa annotations
con = sqlite3.connect("../data/biometa.db")
df = pd.read_sql("SELECT * FROM biometa", con)
con.close()

# %%
# Pull out wing disc data with sex annotations
wing_with_sex = (
    df.query("tissue == 'wing disc' and (sex == 'male' or sex == 'female')")
    .drop(columns=["cell_line", "diet", "chemical", "radiation", "temperature", "other"])
    .rename(columns={"notes": "pubmed"})
)
# %%
client = MongoClient()
db = client["sramongo"]["ncbi"]
sra = pd.DataFrame(
    db.aggregate(
        [
            {"$match": {"BioSample.accn": {"$in": wing_with_sex.biosample.to_list()}}},
            {
                "$project": {
                    "_id": 0,
                    "srx": 1,
                    "srp": "$study.accn",
                    "bioproject": "$BioProject.accn",
                    "biosample": "$BioSample.accn",
                    "geo": "$study.geo",
                    "title": "$BioSample.title",
                    "attributes": "$BioSample.attributes",
                }
            },
        ]
    )
)
client.close()

# %%
def format_author_attributes(sr: pd.Series):
    attributes = "\n".join([f"{attr['name']}: {attr['value']}" for attr in sr.attributes])
    sr["author attributes"] = f"{sr.title}\n{attributes}"
    del sr["title"]
    del sr["attributes"]
    return sr


sample_table = sra.apply(format_author_attributes, axis=1).merge(
    wing_with_sex, left_on="biosample", right_on="biosample"
)
sample_table.to_excel("../output/notebook/2020-06-23_sex_specific_wing_discs.xlsx", index=False)


# %%
srxs = sample_table.srx.unique()
gene_counts = pd.pivot(
    pd.concat(
        [
            pd.read_parquet(f"../output/rnaseq-wf/gene_counts/{srx}.parquet")
            for srx in srxs
            if Path(f"../output/rnaseq-wf/gene_counts/{srx}.parquet").exists()
        ]
    ).reset_index(),
    index="FBgn",
    columns="srx",
    values="count",
)
gene_counts.to_csv("../output/notebook/2020-06-23_sex_specific_wing_discs_counts.tsv", sep="\t")

