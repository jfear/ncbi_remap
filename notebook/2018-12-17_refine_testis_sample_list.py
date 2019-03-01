# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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
sys.path.insert(0, '../lib')
from ncbi_remap.notebook import Nb
from ncbi_remap.plotting import make_figs

# Setup notebook
nbconfig = Nb.setup_notebook()

# Connect to data store
store = pd.HDFStore('../output/sra.h5', mode='r')

# %%
# Get a list of Testis Libraries
metadata = pd.read_csv('../output/geo-wf/rnaseq_metadata.tsv', sep='\t', index_col=0).query('tissue == "testis"')
testis_srxs = metadata.index.unique().tolist()
metadata.head()

# %%
# Get a list of srxs that have 75% of reads well stranded
stranded = store.select('prealn/workflow/collectrnaseqmetrics/second', where='srx == testis_srxs & PCT_CORRECT_STRAND_READS >= .75').index.get_level_values('srx').unique().tolist()
with open('../output/notebook/2018-12-17_testis_stranded_libraries.txt', 'w') as fh:
    fh.write('\n'.join(stranded))

# %%
# Get list of PE testis libraries
pe = store.select('layout', 'srx == testis_srxs & layout == "PE"').index.get_level_values('srx').unique().tolist()
with open('../output/notebook/2018-12-17_testis_pe_libraries.txt', 'w') as fh:
    fh.write('\n'.join(pe))

# %%

