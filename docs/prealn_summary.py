"""Get various counts from the pre-alignment workflow."""
#%%
import pandas as pd
import os

try:
    os.chdir(os.path.join(os.getcwd(), "docs"))
    print(os.getcwd())
except:
    pass

#%%
store = pd.HDFStore("../output/sra.h5")

#%%[markdown]
# # Total Samples Processed

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
            sr['srx'] = srx
            sr['srr'] = srr
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
