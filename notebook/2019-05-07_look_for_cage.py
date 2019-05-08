# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.0.3
#   kernelspec:
#     display_name: Python [conda env:ncbi_remap]
#     language: python
#     name: conda-env-ncbi_remap-py
# ---

# %%
import os
import sys
from pathlib import Path

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Project level imports
from ncbi_remap.notebook import Nb

# %%
pd.options.display.max_rows = 200

# %%
# Setup notebook
nbconfig = Nb.setup_notebook()

# %%
# Connect to data store
store = pd.HDFStore('../output/sra.h5', mode='r')

# %%
from pymongo import MongoClient

host = 'localhost'
mongoClient = MongoClient(host=host, port=27017)
db = mongoClient['sramongo']
ncbi = db['ncbi']

# %%
df = pd.read_parquet('../output/metadata-wf/select_library_strategy.parquet')

# %%

# %%
sel = pd.DataFrame(list(ncbi.aggregate([
    {
        '$project': {
            '_id': 0,
            'library_selection': '$library_selection',
            'srx': '$srx'
        }
    }
])))

# %%
cage_srxs = sel.query('library_selection == "CAGE"').srx.values

# %%
cage_srxs.shape

# %%
metadata = pd.read_csv('../output/geo-wf/rnaseq_metadata.tsv', sep='\t', index_col=0).drop('raw file', axis=1)

# %%
df.reindex(cage_srxs).join(metadata).to_csv('~/Downloads/fear_SRA_CAGE_samples.csv')

# %%
