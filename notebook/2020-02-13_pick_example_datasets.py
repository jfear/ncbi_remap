# %%
import pickle

import numpy as np
import pandas as pd

from ncbi_remap.queue import Queue

# %%
rnaseq = Queue("../output/library_strategy-wf/rnaseq.pkl").sample_table

# %%
# Below is how I found candidate samples.
# store = pd.HDFStore("../output/sra.h5", mode="r")
# layout = store["layout"]
# libsize = (
#     store.select("prealn/workflow/fastq", columns=["libsize_R1", "libsize_R2"])
#     .max(axis=1)
#     .rename("libsize")
# )
# data = pd.concat([libsize, layout], axis=1, sort=False).reset_index()
# store.close()


# df = rnaseq.merge(data, on=["srx", "srr"], how="inner")

# df.groupby("layout").apply(lambda g: g.query("libsize < 6_000_000").sample(n=2))

# %%
se = "SRX014472"
pe = "SRX014459"
srxs = [se, pe]

# %%
srrs = rnaseq.query(f"srx == {srxs}").srr.unique().tolist()

# %%
pickle.dump(srxs, open("../output/example_srxs.pkl", "wb"))
pickle.dump(srrs, open("../output/example_srrs.pkl", "wb"))


# %%
