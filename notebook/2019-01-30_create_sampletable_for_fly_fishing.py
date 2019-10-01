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

# %% [markdown]
# # Sample Table for Fly Fishing.

# %% [markdown]
# I setup a new repository [fly_fishing](https://github.com/jfear/fly_fishing) to run Salmon optimizations. I want to use ovary and testis data from the SRA for this process. Here I need to generate a set of ~100 high quality samples.

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
# Setup notebook
nbconfig = Nb.setup_notebook()

# Connect to data store
store = pd.HDFStore('../output/sra.h5', mode='r')

# %%
# Pull out head data
metadata = pd.read_csv('../output/geo-wf/rnaseq_metadata.tsv', sep='\t', index_col=0)
srxs = metadata.query('tissue == ["ovary", "testis"]').index.tolist()

# %%
metadata.head()

# %%
# Grab only the stranded head data.
pct_strand = store.select('/prealn/workflow/collectrnaseqmetrics/second', where='srx == srxs')['PCT_CORRECT_STRAND_READS'].to_frame().reset_index()
stranded = pct_strand.loc[pct_strand.PCT_CORRECT_STRAND_READS > .9, ['srx', 'srr']]

# %%
# Add on library layout info for mapping
layout = store.select('layout', where='srx == srxs & layout == ["PE", "SE"]').to_frame().reset_index()
stranded_w_layout = stranded.merge(layout, on=['srx', 'srr'], how='left').dropna()
stranded_w_layout.columns = ['samplename', 'Run', 'layout']

# %%
tissue = metadata[['tissue', 'study']].reindex(srxs).reset_index()

# %%
stranded_w_layout_w_tissue = stranded_w_layout.merge(tissue, how='left', left_on='samplename', right_on='sample_name')

# %%
stranded_w_layout_w_tissue.to_csv('../output/notebook/2019-01-30_186_testis_ovary_samples_for_fly_fishing.tsv', sep='\t', index=False)
