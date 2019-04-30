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
from pickle import loads

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# %%
df = pd.read_parquet('/home/fearjm/scratch/ncbi_remap/rnaseq_correlations_clustered.parquet')

# %%
fig, ax = plt.subplots(figsize=(20, 20))
ax.matshow(df.values, cmap='inferno', rasterized=True)

# %%

# %%
