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
df.to_csv('/home/fearjm/Downloads/fear_sra_ML_categories.csv')

# %%
