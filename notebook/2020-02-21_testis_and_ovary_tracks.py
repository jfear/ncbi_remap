# %%
from pathlib import Path
from subprocess import check_call
import shlex
import pickle

import numpy as np
import pandas as pd

# %%
# Get a list of RNA-Seq
rnaseq = pickle.load(open("../output/library_strategy-wf/rnaseq.pkl", "rb"))

# %%
# Get current metadata and make sure only in current RNA-Seq
bio = (
    pd.read_csv("../output/geo-wf/rnaseq_metadata.tsv", sep="\t", index_col=0)
    .pipe(lambda x: x[x.index.isin(rnaseq)])
)

# %%
# Pull out testis and ovary
testis = bio[bio.tissue.fillna("").str.contains("testis")].index.tolist()
ovary = bio[bio.tissue.fillna("").str.contains("ovary")].index.tolist()

# %%
# Save SRXs that are testis
with open("../output/notebook/2020-02-21_testis.txt", "w") as fh:
    fh.write("\n".join(testis))
    fh.write("\n")

# %%
# Save SRXs that are ovary
with open("../output/notebook/2020-02-21_ovary.txt", "w") as fh:
    fh.write("\n".join(ovary))
    fh.write("\n")

# %%
# Download BigWigs for a sample
tout = Path("/home/fearjm/tmp/testis")
tout.mkdir(exist_ok=True)

oout = Path("/home/fearjm/tmp/ovary")
oout.mkdir(exist_ok=True)

np.random.seed(42)
idx = np.random.randint(0, 100, size=10).tolist()

for i, srx in enumerate(testis):
    if i not in idx:
        continue
    first = f"/data/MiegNCBI/ncbi_remap/output/rnaseq-wf/samples/{srx}/{srx}.flybase.first.bw"
    first_out = Path(tout, f"{srx}.testis.first.bw").as_posix()
    check_call(shlex.split(f"scp helix:{first} {first_out}"))

    second = f"/data/MiegNCBI/ncbi_remap/output/rnaseq-wf/samples/{srx}/{srx}.flybase.second.bw"
    second_out = Path(tout, f"{srx}.testis.second.bw").as_posix()
    check_call(shlex.split(f"scp helix:{second} {second_out}"))

for i, srx in enumerate(ovary):
    if i not in idx:
        continue
    first = f"/data/MiegNCBI/ncbi_remap/output/rnaseq-wf/samples/{srx}/{srx}.flybase.first.bw"
    first_out = Path(oout, f"{srx}.ovary.first.bw").as_posix()
    check_call(shlex.split(f"scp helix:{first} {first_out}"))

    second = f"/data/MiegNCBI/ncbi_remap/output/rnaseq-wf/samples/{srx}/{srx}.flybase.second.bw"
    second_out = Path(oout, f"{srx}.ovary.second.bw").as_posix()
    check_call(shlex.split(f"scp helix:{second} {second_out}"))


# %%
# paper
pout = Path("/home/fearjm/tmp/GSE93699")
pout.mkdir(exist_ok=True)

paper_srxs = bio[(bio["GEO Experiment"] == "GSE93699")].index.tolist()
for srx in paper_srxs:
    first = f"/data/MiegNCBI/ncbi_remap/output/rnaseq-wf/samples/{srx}/{srx}.flybase.first.bw"
    first_out = Path(pout, f"{srx}.first.bw").as_posix()
    check_call(shlex.split(f"scp helix:{first} {first_out}"))

    second = f"/data/MiegNCBI/ncbi_remap/output/rnaseq-wf/samples/{srx}/{srx}.flybase.second.bw"
    second_out = Path(pout, f"{srx}.second.bw").as_posix()
    check_call(shlex.split(f"scp helix:{second} {second_out}"))



# %%
print("\n".join(bio[(bio["GEO Experiment"] == "GSE93699")].sort_index().title))

# %%
