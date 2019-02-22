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
# # Plots for BSC

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
nbconfig = Nb.setup_notebook(ref_dir='../references/lcdb-references')

# %%
# Connect to data store
store = pd.HDFStore('../output/sra.h5', mode='r')

# %%
from pymongo import MongoClient
try:
    with open('../output/.mongodb_host', 'r') as fh:
        host = fh.read().strip()
except FileNotFoundError:
    host = 'localhost'

mongoClient = MongoClient(host=host, port=27017)
db = mongoClient['sra']
ncbi = db['ncbi']

# %%
df = (
    pd.DataFrame(list(
        ncbi.aggregate([
            {'$unwind': {'path': '$sra.run'}},
            {
                '$project': {
                    '_id': False,
                    'srx': '$_id',
                    'srr': '$sra.run.run_id',
                    'date': '$sra.run.load_date',
                    'size_MB': '$sra.run.size_MB'
                }
            },
        ])
    ))
    .set_index(['srx', 'srr'])
    .assign(date=lambda df: pd.to_datetime(df.date))
    .sort_values('date')
    .assign(cum_sum_TB = lambda df: df.size_MB.cumsum() / 1e6)
    .assign(year=lambda df: df.date.dt.year)
    .dropna()
    .assign(year=lambda df: df.year.astype(int))
)

# %%
df.head()

# %%
fig, ax = plt.subplots(figsize=(10, 8))
df.plot('date', 'cum_sum_TB', ax=ax, legend=False)
ax.fill_between(df.date.dt.to_pydatetime(), 0, df['cum_sum_TB'])
ax.margins(0)
ax.set_xlabel('Date', fontweight='bold')
ax.set_ylabel('Date', fontweight='bold')
ax.set_title('D. melanogaster Data In the SRA', fontweight='bold', family='Cambria')
fig.autofmt_xdate(rotation=0, ha='center')


# %%



# %%



# %%



# %%



# %%



# %%



# %%
strategy = (
    pd.read_parquet('../output/metadata-wf/select_library_strategy.parquet')
    .rename(columns={'Fear_et_al_library_strategy': 'library_strategy'})
    .library_strategy
    .pipe(lambda x: x[~x.str.contains('\|')])
    #.pipe(lambda x: x[x.isin(['RNA-Seq', 'WGS'])])
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
url = '../output/notebook/2019-02-19_lib_strat_features.json'
data.sample(1_000).to_json(url, orient='records')

# %%
def test_theme():
    return {
        'config': {
            'axis': {
                'titleFont': 'Arial',
                'titleFontSize': 8,
                'labelFont': 'Arial',
                'labelFontSize': 3,
            },
            'title': {
                'font': 'Cambria',
                'fontSize': 8
            }
        }
    }

alt.themes.register('test_theme', test_theme)
alt.themes.enable('test_theme')

# %%
chart = (
    alt.Chart(url)
    .mark_circle()
    .encode(
        alt.X(alt.repeat('column'), type='quantitative'),
        alt.Y(alt.repeat('row'), type='quantitative'),
        color=alt.Color('library_strategy:N', legend=alt.Legend(title='Library Strategy'))
    )
    .properties(
        width=90,
        height=90
    )
    .repeat(
        column=feature_names,
        row=feature_names
    )
)

# %%
chart.savechart('~/Downloads/test.png', format='png')

# %%



# %%



# %%



# %%



# %%



# %%



