"""Get various counts from the pre-alignment workflow."""
#%%
import os
from datetime import datetime

import pandas as pd
from pymongo import MongoClient

try:
    os.chdir(os.path.join(os.getcwd(), "docs"))
    print(os.getcwd())
except:
    pass

#%%[markdown]
# # What we did not download

#%%
client = MongoClient()
db = client["sramongo"]
ncbi = db["ncbi"]

total_samples_in_db = pd.DataFrame(
    ncbi.aggregate(
        [
            {"$match": {"sra_create_date": {"$lte": datetime(2019, 3, 31)}}},
            {"$unwind": {"path": "$runs"}},
            {"$match": {"runs.srr": {"$exists": True}}},
            {"$project": {"_id": 0, "srx": "$srx", "srr": "$runs.srr"}},
        ]
    )
)

print(f"SRXs in DataBase: {total_samples_in_db.srx.unique().shape[0]:,}")
print(f"SRRs in DataBase: {total_samples_in_db.srr.unique().shape[0]:,}")

too_few_reads = pd.DataFrame(ncbi.aggregate(
    [
        {"$match": {"sra_create_date": {"$lte": datetime(2019, 3, 31)}}},
        {"$unwind": {"path": "$runs"}},
        {"$match": {"runs.srr": {"$exists": True}}},
        {
            "$project": {
                "_id": 0,
                "srx": 1,
                "srr": "$runs.srr",
                "nspots": {
                    "$cond": [{"$eq": ["$runs.nspots", ""]}, 0, {"$ifNull": ["$runs.nspots", 0]}]
                },
                "read_count_r1": {"$ifNull": ["$runs.read_count_r1", 0]},
                "read_count_r2": {"$ifNull": ["$runs.read_count_r2", 0]},
            }
        },
        {
            # Keep samples with > 1,000 reads
            "$match": {
                "$and": [
                    {
                        # Keep samples with all 0 b/c we just don't know
                        "$or": [
                            {"nspots": {"$gt": 0}},
                            {"read_count_r1": {"$gt": 0}},
                            {"read_count_r2": {"$gt": 0}},
                        ]
                    },
                    {
                        "$and": [
                            {"nspots": {"$lte": 1000}},
                            {"read_count_r1": {"$lte": 1000}},
                            {"read_count_r2": {"$lte": 1000}},
                        ]
                    },
                ]
            }
        },
        {"$project": {"_id": 0, "srx": "$srx", "srr": "$srr"}},
    ]
))

print(f"SRRs not downloaded b/c < 1,000 reads: {too_few_reads.srr.unique().shape[0]:,}")

too_short_reads = pd.DataFrame(ncbi.aggregate(
    [
        {"$match": {"sra_create_date": {"$lte": datetime(2019, 3, 31)}}},
        {"$unwind": {"path": "$runs"}},
        {"$match": {"runs.srr": {"$exists": True}}},
        {
            "$project": {
                "_id": 0,
                "srx": 1,
                "srr": "$runs.srr",
                "read_len_r1": {"$ifNull": ["$runs.read_len_r1", 0]},
                "read_len_r2": {"$ifNull": ["$runs.read_len_r2", 0]},
            }
        },
        {
            # Keep samples with > 25 bp reads lens
            "$match": {
                "$and": [
                    {
                        # Keep samples with all 0 b/c we just don't know
                        "$or": [{"read_len_r1": {"$gt": 0}}, {"read_len_r2": {"$gt": 0}}]
                    },
                    {"$and": [{"read_len_r1": {"$lt": 25}}, {"read_len_r2": {"$lt": 25}}]},
                ]
            }
        },
        {"$project": {"_id": 0, "srx": "$srx", "srr": "$srr"}},
    ]
))

print(f"SRRs not downloaded b/c < 25 bp long: {too_short_reads.srr.unique().shape[0]:,}")

#%%[markdown]
# # Total Samples Processed

#%%
store = pd.HDFStore("../output/sra.h5")

#%%
# Ignore samples still in the alignment queue
srx_aln_queue = store["aln/queue"].srx.tolist()
srx_aln_complete_cnt = len(store["aln/complete"].srx.unique().tolist())
srr_aln_complete_cnt = len(store["aln/complete"].srr.unique().tolist())

# Complete prealn-workflow
srr_cnt = len(store.select("/prealn/complete", where="srx != srx_aln_queue").srr.tolist())
srx_cnt = len(store.select("/prealn/complete", where="srx != srx_aln_queue").srx.unique().tolist())

# Download Bad
srr_cnt += len(store.select("/prealn/download_bad", where="srx != srx_aln_queue").srr.tolist())
srx_cnt += len(
    store.select("/prealn/download_bad", where="srx != srx_aln_queue").srx.unique().tolist()
)

# Quality Score Bad
srr_cnt += len(
    store.select("/prealn/quality_scores_bad", where="srx != srx_aln_queue").srr.tolist()
)
srx_cnt += len(
    store.select("/prealn/quality_scores_bad", where="srx != srx_aln_queue").srx.unique().tolist()
)

# Alignment Bad
srr_cnt += len(store.select("/prealn/alignment_bad", where="srx != srx_aln_queue").srr.tolist())
srx_cnt += len(
    store.select("/prealn/alignment_bad", where="srx != srx_aln_queue").srx.unique().tolist()
)

# Abi Solid
srr_cnt += len(store.select("/prealn/abi_solid", where="srx != srx_aln_queue").srr.tolist())
srx_cnt += len(
    store.select("/prealn/abi_solid", where="srx != srx_aln_queue").srx.unique().tolist()
)

print(f"SRX Processed: {srx_cnt:,}")
print(f"SRR Processed: {srr_cnt:,}")

print(f"SRX Complete Both Workflows: {srx_aln_complete_cnt:,}")
print(f"SRR Complete Both Workflows: {srr_aln_complete_cnt:,}")

#%%[markdown]
# # Details about problems

#%%[markdown]
# ## Download Bad

#%%
# Get Data
download_bad = store.select("/prealn/download_bad", where="srx != srx_aln_queue")


class downloadBad:
    missing = []
    df = pd.DataFrame(
        columns=[
            "srx",
            "srr",
            "md5_R1",
            "libsize_R1",
            "avgLen_R1",
            "md5_R2",
            "libsize_R2",
            "avgLen_R2",
        ]
    )

    def load(self, srx, srr):
        try:
            sr = pd.read_csv(
                f"../output/prealn-wf/samples/{srx}/{srr}/{srr}.fastq.tsv", sep="\t"
            ).T.squeeze()
            sr["srx"] = srx
            sr["srr"] = srr
            self.df = self.df.append(sr, ignore_index=True)
        except FileNotFoundError:
            self.missing.append([srx, srr])


db = downloadBad()
for i, (srx, srr) in download_bad.iterrows():
    db.load(srx, srr)

print(f"SRRs could not be downloaded: {len(db.missing):,}")

#%%
db.df

#%%
fastq = store["/prealn/workflow/fastq"]

#%%
srr_complete = store["/prealn/complete"].srr.tolist()
srr_alignment_bad = store["/prealn/alignment_bad"].srr.tolist()
srr_quality_scores_bad = store["/prealn/quality_scores_bad"].srr.tolist()
srr_abi = store["/prealn/abi_solid"].srr.tolist()

#%%[markdown]
# Download Bad

#%%
srr_download_bad = store["/prealn/download_bad"].srr.tolist()

#%%
fastq.query(f"srr == {srr_complete}")

#%%
