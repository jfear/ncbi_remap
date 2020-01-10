# %%
import pandas as pd
from yaml import full_load


# %%
SRA_TERMS = full_load(open("../config/sra_terms.yaml"))

# %%
free_text = pd.read_parquet("../output/library_strategy-wf/free_text_library_strategy.parquet").squeeze()

# %%
pd.DataFrame([
    (i, x)
    for i, x in free_text.dropna().items()
    if x not in SRA_TERMS["library_strategy"]
], columns=["srx", "library_strategy"]).to_csv("/home/fearjm/bob.tsv", sep="\t", index=False)

# %%
free_text.to_csv("/home/fearjm/bob2.tsv", sep="\t")