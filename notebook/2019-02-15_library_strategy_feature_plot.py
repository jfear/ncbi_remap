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
import altair as alt

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
strategy = (
    pd.read_parquet('../output/metadata-wf/select_library_strategy.parquet')
    .rename(columns={'Fear_et_al_library_strategy': 'library_strategy'})
    .library_strategy
    .pipe(lambda x: x[~x.str.contains('\|')])
    .pipe(lambda x: x[x.isin(['RNA-Seq', 'WGS'])])
)

# %%
feature_names = (
    pd.read_csv('../output/metadata-wf/random_forest_library_strategy_feature_importance.tsv', sep='\t', header=None, names=['feature', 'importance'])
    .sort_values('importance', ascending=False)
    .head(10)
    .feature
    .values
    .tolist()
)

# %%
data = (
    pd.read_parquet('../output/metadata-wf/build_library_strategy_feature_set.parquet', columns=feature_names)
    .join(strategy, how='inner')
)

# %%
data.sample(10_000).to_json('bob.json', orient='records')

# %%
chart = (
    alt.Chart('bob.json')
    .mark_circle()
    .encode(
        alt.X(alt.repeat('column'), type='quantitative'),
        alt.Y(alt.repeat('row'), type='quantitative'),
        color='library_strategy:N'
    )
    .repeat(
        column=feature_names,
        row=feature_names
    )
)

# %%
chart


# %%


