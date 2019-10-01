# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.0.0
#   kernelspec:
#     display_name: Python [conda env:ncbi_remap]
#     language: python
#     name: conda-env-ncbi_remap-py
# ---

# %%
import os
import sys
from pathlib import Path
from more_itertools import chunked

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from dask.distributed import Client, as_completed
from dask import dataframe as dd
from dask import delayed

# Project level imports
from ncbi_remap.notebook import Nb

# %%
# Setup notebook
nbconfig = Nb.setup_notebook()

# %%
# Start dask cluster
dask_client = Client()
dask_client


# %%
def read_gene_counts(srx):
    return  pd.read_parquet(f'/home/fearjm/scratch/ncbi_remap/output/aln-wf/gene_counts/{srx}.parquet', columns=['count'])['count'].rename(srx)


def chunk_and_run(iterable, chunk_size=1_000):
    chunks = chunked(iterable, chunk_size)
    for chunk in chunks:
        futures = dask_client.map(read_gene_counts, chunk)
        yield pd.concat(dask_client.gather(futures), axis=1)
        
def build_table(iterable, chunk_size=1000):
    dfs = chunk_and_run(iterable, chunk_size)
    df = next(dfs)
    for _df in dfs:
        df = df.merge(_df, on='FBgn')
    return df


# %%
rnaseq = (
    pd.read_parquet('../output/metadata-wf/select_library_strategy.parquet')
    .rename(columns={'Fear_et_al_library_strategy': 'strategy'})
    .strategy
    .pipe(lambda x: x[x == 'RNA-Seq'])
    .index.unique().tolist()
)

# %%
df = build_table(rnaseq)
df.to_parquet('/home/fearjm/scratch/ncbi_remap/rnaseq_wide.parquet')

# %%
df.rank().to_parquet('/home/fearjm/scratch/ncbi_remap/rnaseq_wide_ranks.parquet')

# %%
