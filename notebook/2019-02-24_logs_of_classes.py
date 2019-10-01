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
from collections import Counter

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
db = mongoClient['sramongo']
ncbi = db['ncbi']

# %%
rnaseq = (
    pd.read_parquet('../output/metadata-wf/select_library_strategy.parquet')
    .rename(columns={'Fear_et_al_library_strategy': 'strategy'})
    .strategy
    .pipe(lambda x: x[x == 'RNA-Seq'])
    .index.unique().tolist()
)

# %%
query = ncbi.aggregate([
    {
        '$match': {
            '_id': {'$in': rnaseq}
        }
    },
    {
        '$unwind': {
            'path': '$biosample'
        }
    },
    {
        '$project': {
            'metadata': '$biosample.attributes'
        }
    },
])

names = []
for record in query:
    for dat in record['metadata']:
        names.append(dat['name'])

# %%

# %%
name_count = pd.DataFrame(Counter(names).most_common(), columns=['names', 'cnt'])
total_size = name_count.shape[0]
name_count.loc[name_count.cnt < 200, 'names'] = '...'
name_count = name_count.groupby('names').cnt.sum()
name_count['...'] = 10

order = ['source_name', 'tissue', 'sex', 'dev_stage', 'strain', 'genotype', 'cell_line', 'age', 'RNAi reagent from harvard drosophila rnai screening center', 'RNAi target gene name',
         'treatment', 'replicate', 'cell_type', 'sample_name', 'Stage', 'time point', 'genetic background', 'sample_type', 'biomaterial_provider', 'isolate', 'developemntal stage',
         'drosdel id', 'bloomington stock center id', 'carrier species', 'breed', 'cultivar', 'ecotype', 'temp', 'flag', 'drosdel deficiency', 'crosses', 'x chromosome dose', 'collection_date',
         'illumina barcode', 'plate and well id', 'flybase species id', 'ercc infomation', 'geo_loc_name', 'ENA-FIRST-PUBLIC', 'ENA-LAST-UPDATE' , '...',
]

name_count = name_count.reindex(order).to_frame().reset_index()
diff = total_size - (name_count.shape[0] - 1)

mapper = {
    'RNAi reagent from harvard drosophila rnai screening center': 'RNAi ... screening center',
    '...': f'{diff} Additional Class Names'
}

name_count.names = name_count.names.replace(mapper)


# %%
fig, ax = plt.subplots(1, 1, figsize=(8, 5))
sns.barplot('names', 'cnt', data=name_count, color='C0', ax=ax)
plt.setp(ax.get_xticklabels(), rotation=75, ha='right', va='top')

fig.savefig('../output/notebook/2019-02-24_class_counts.svg', bbox_inches='tight')

# %%
query = ncbi.aggregate([
    {
        '$match': {
            '_id': {'$in': rnaseq}
        }
    },
    {
        '$unwind': {
            'path': '$biosample'
        }
    },
    {
        '$project': {
            'metadata': '$biosample.attributes'
        }
    },
])

values = []
for record in query:
    for dat in record['metadata']:
        if dat['name'] == 'tissue':
            values.append(dat['value'])

# %%
value_count = pd.DataFrame(Counter(values).most_common(), columns=['values', 'cnt'])

total_size = value_count.shape[0]
value_count.loc[value_count.cnt < 30, 'values'] = '...'
value_count = value_count.groupby('values').cnt.sum()
value_count['...'] = 1

value_count = value_count.sort_values(ascending=False).to_frame().reset_index()
diff = total_size - (value_count.shape[0] - 1)

mapper = {
    'Entire individual except for head': 'Entire ... for head',
    'Periphery (entire body except head)': 'Periphery ...',
    'abdomen without digestive or reproductive system': 'abdomen ... system',
    '...': f'{diff} Additional Values'
}

value_count['values'] = value_count['values'].replace(mapper)

# %%
fig, ax = plt.subplots(1, 1, figsize=(8, 5))
sns.barplot('values', 'cnt', data=value_count, color='C0', ax=ax)
plt.setp(ax.get_xticklabels(), rotation=75, ha='right', va='top');
ax.text(0.5, 0.5, f'{total_size} values', transform=ax.transAxes)
fig.savefig('../output/notebook/2019-02-24_tissue_value_counts.svg', bbox_inches='tight')

# %%
