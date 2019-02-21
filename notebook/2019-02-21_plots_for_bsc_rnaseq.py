# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 0.8.6
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
from ncbi_remap.plotting import make_figs

# %%
# Setup notebook
nbconfig = Nb.setup_notebook()

# %%
# Connect to data store
store = pd.HDFStore('../output/sra.h5', mode='r')

# %%
rnaseq_srxs = (
    pd.read_parquet('../output/metadata-wf/select_library_strategy.parquet')
    .query('Fear_et_al_library_strategy == "RNA-Seq"')
    .index.tolist()
)

# %%
df = pd.read_parquet('../output/aln-downstream-wf/gene_counts_wide_tpm.parquet', columns=rnaseq_srxs)

# %%
df.shape

# %%


# %%


# %%


# %%


# %%


# %%


# %%


# %%

