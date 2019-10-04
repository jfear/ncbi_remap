"""ChIP-Seq for Claire

> email: 10-03-2019
> It would be great for us to get the list of accessions and the metadata. I
> am sending you the samples covered in modERN and there are only a few samples
> from cell lines so it would be great to include what you have already done
> with S2 cells.
"""
#%%
from bson import json_util
from json import dump
import pandas as pd
from pymongo import MongoClient


#%%
client = MongoClient()
db = client["sramongo"]
ncbi = db["ncbi"]

#%%
print(dumps(json_util.dumps(ncbi.find_one())))

#%%
# Get From Database
srx2srr = pd.DataFrame(
    ncbi.aggregate(
        [
            {"$match": {"library_strategy": "ChIP-Seq"}},
            {
                "$project": {
                    "_id": 0,
                    "srx": 1,
                    "srrs": {
                        "$reduce": {
                            "input": "$runs.srr",
                            "initialValue": "",
                            "in": {
                                "$concat": [
                                    "$$value",
                                    {"$cond": [{"$eq": ["$$value", ""]}, "", "|"]},
                                    "$$this",
                                ]
                            },
                        }
                    },
                }
            },
        ]
    )
).set_index("srx")

srx2attrs = pd.DataFrame(
    ncbi.aggregate(
        [
            {"$match": {"library_strategy": "ChIP-Seq"}},
            {"$unwind": {"path": "$sample.attributes"}},
            {
                "$project": {
                    "_id": 0,
                    "srx": 1,
                    "name": "$sample.attributes.name",
                    "value": "$sample.attributes.value",
                }
            },
        ]
    )
).pivot_table(index="srx", columns="name", values="value", aggfunc="first")

srx2srr.join(srx2attrs).to_csv("output/notebook/2019-10-04_chipseq_for_claire.tsv", sep='\t')

with open("output/notebook/2019-10-04_chipseq_for_claire.json", "w") as fh:
    fh.write(
        json_util.dumps(
            ncbi.find({"library_strategy": "ChIP-Seq"}),
            json_options=json_util.RELAXED_JSON_OPTIONS,
            indent=2,
        )
    )
    fh.write("\n")

#%%
df = pd.read_csv("../s2cell-prior/chipseq-wf/config/sampletable_all.tsv", sep="\t").drop(
    ["orig_filename", "biological_material"], axis=1
)
df.columns = ["srr", "antibody", "replicate", "srx"]
df.set_index("srx").to_csv("output/notebook/2019-10-04_chipseq_for_claire_s2_samples.tsv", sep="\t")

#%%
