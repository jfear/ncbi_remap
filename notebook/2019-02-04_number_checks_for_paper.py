# -*- coding: utf-8 -*-
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

from joblib import memory

# Project level imports
sys.path.insert(0, '../lib')
from ncbi_remap.notebook import Nb
from ncbi_remap.plotting import make_figs
from ncbi_remap.config import memory

# Setup notebook
nbconfig = Nb.setup_notebook()

# %%
# Connect to data store
store = pd.HDFStore('../output/sra.h5', mode='r')

# %%
from pymongo import MongoClient
mongoClient = MongoClient(host='localhost', port=27017)
db = mongoClient['sra']
ncbi = db['ncbi']

# %%
store.root

# %% [markdown]
# ## Pre-alignment and Alignment Workflow

# %%
srrs = store.ids.srr.unique().shape[0]
srxs = store.ids.srx.unique().shape[0]
biosamples = len(set([x['biosample'][0]['biosample_accn'] for x in ncbi.find({'biosample.biosample_accn': {'$exists': True}}, {'biosample.biosample_accn': True, '_id': False})]))
print(f'As of June 2018, there were {srrs:,} *D. melanogaster* SRA runs (SRRs), which is equivalent to {srxs:,} unique SRA experiments (SRXs) or {biosamples:,} BioSamples.')

# %%
complete = store['prealn/complete'].srr.unique().shape[0]
print(f'A total of {complete:,} SRRs completed the pre-alignment workflow.')

# %%
download_bad = store['prealn/download_bad'].srr.unique().shape[0]
downloaded = store.ids.srr.unique().shape[0] - store['prealn/queue'].srr.unique().shape[0] - download_bad
print(f'We used fastq-dump (@REF) to download {downloaded:,} SRRs from the SRA; there were {download_bad:,} SRRs that we were unable to download from the SRA.')

# %%
se = store.select('layout', where='layout == "SE"').shape[0]
print(f'The majority of SRRs were from single end libraries ({se:,}).')

# %%
pe = store.select('layout', where='layout == "PE"').shape[0]
k1 = store.select('layout', where='layout == "keep_R1"').shape[0]
k2 = store.select('layout', where='layout == "keep_R2"').shape[0]
print(f'There were {pe + k1 + k2:,} libraries labeled as pair-ended in the SRA, but only {pe:,} had useful read pairs. ')
print(f'There were a number of SRRs that were uploaded as pair-end, but either R1 ({k2:,}) or R2 ({k1:,}) looked to be sample barcodes and not true read pairs. ')

# %%
corrupt = store['prealn/quality_scores_bad'].srr.unique().shape[0]
print(f'We found {corrupt:,} SRRs that were corrupted, with either partial reads or non-matching quality scores.')

# %%
abi = store['prealn/abi_solid'].srr.unique().shape[0]
print(f'Finally, we excluded {abi:,} SRRs which were sequenced using ABI Solid technology (@REF).')

# %%
align = store['prealn/alignment_bad'].srr.unique().shape[0]
print(f'We also removed {align:,} SRRs which had < 50% of reads mapping to the *D. melanogaster* genome.')

# %%
sense = store.select('prealn/workflow/collectrnaseqmetrics/first', where='PCT_CORRECT_STRAND_READS >= .75').index.get_level_values('srr').unique().shape[0]
antisense = store.select('prealn/workflow/collectrnaseqmetrics/second', where='PCT_CORRECT_STRAND_READS >= .75').index.get_level_values('srr').unique().shape[0]
unstranded = complete - (sense + antisense)
print(f'A library was considered stranded in the alignment workflow if ≥ 75% of reads mapped to either the sense (n = {sense:,}) or anti-sense strand (n = {antisense:,}), all other libraries were considered unstranded (n = {unstranded:,}).')

# %%
num_srr_per_srx = store['aln/complete'].groupby('srx').size().value_counts().map(lambda x: f'{x:,}').to_frame()
num_srr_per_srx.columns = ['Number of SRXs']
num_srr_per_srx.index.name = 'Number of SRRs'
num_srr_per_srx

# %%
more_than_2_srr = num_srr_per_srx.iloc[2:, 0].astype(int).sum()
print(f'The majority of SRXs had a single SRR (n = {num_srr_per_srx.iloc[0, 0]}), with {num_srr_per_srx.iloc[1, 0]} SRXs having 2 SRRs, and {more_than_2_srr:,} SRXs having more than 2 SRRs (@Table_S1).')

# %%
@memory.cache
def get_corr(srx):
    return (srx, pd.read_parquet(f'../output/prealn-wf/gene_counts/{srx}.parquet').set_index('srr', append=True)['count'].unstack().corr(method='spearman').min().min())

corrs = pd.DataFrame([
    get_corr(srx)
    for srx in store['aln/complete'].groupby('srx').size().pipe(lambda x: x[x >= 2]).index.unique()
], columns=['srx', 'min_corr']).set_index('srx')

# %%
low_corr = corrs.query('min_corr <= 0.85').shape[0]
print(f'We merged only highly correlated SRRs (Spearman\'s rho ≥ 0.85); the {low_corr:,} SRXs with poor correlation were removed (@Figure_2D). ')

# %% [markdown]
# ## Validation of Library Strategy 
#

# %%
lib_strat = pd.DataFrame([
    (record['_id'], record['sra']['experiment']['library_strategy'])
    for record in ncbi.find({'_id': {'$in': store['aln/complete'].srx.unique().tolist()}, 'sra.experiment.library_strategy': {'$exists': True}}, {'sra.experiment.library_strategy': True})
], columns=['srx', 'library_stategy']).set_index('srx').library_stategy.value_counts()

# %%
print('SRA Library Strategy')
lib_strat

# %%
print(f'For example, the data presented here consist of RNA-Seq (n = {lib_strat["RNA-Seq"]:,}), EST (n = {lib_strat["EST"]:,}), ChIP-Seq (n = {lib_strat["ChIP-Seq"]:,}), or whole genome shotgun sequencing (WGS; n = {lib_strat["WGS"]:,}) (@Table_S1).')


# %%
print(f'For example, a new technology may not exist in the curated list of library strategies, so the authors may classify their samples using the closest term or as OTHER (n = {lib_strat["OTHER"]:,}).')

# %%
lib_strat_features = pd.read_parquet('../output/metadata-wf/build_library_strategy_feature_set.parquet').columns
print(f'To validate library strategy, we selected {lib_strat_features.shape[0]:,} features from the pre-alignment and alignment workflows.')

# %%
my_lib = pd.read_parquet('../output/metadata-wf/select_library_strategy.parquet')
fear_counts = my_lib.Fear_et_al_library_strategy.value_counts()

# %%
fear_counts.to_frame()

# %%
print(f'We unambiguously identified samples as RNA-Seq (n = {fear_counts["RNA-Seq"]:,}), EST (n = {fear_counts["EST"]:,}), ChIP-Seq (n = {fear_counts["ChIP-Seq"]:,}), or WGS (n = {fear_counts["WGS"]:,}).')

# %%
ncbi_strat = [
    "AMPLICON",
    "ATAC-Seq",
    "Bisulfite-Seq",
    "ChIP-Seq",
    "CLONE",
    "CLONEEND",
    "CTS",
    "DNase-Hypersensitivity",
    "EST",
    "FAIRE-seq",
    "FINISHING",
    "FL-cDNA",
    "HI-C",
    "MBD-Seq",
    "MeDIP-Seq",
    "miRNA-Seq",
    "MNase-Seq",
    "MRE-Seq",
    "ncRNA-Seq",
    "OTHER",
    "POOLCLONE",
    "RIP-Seq",
    "RNA-Seq",
    "Synthetic-Long-Read",
    "SELEX",
    "Tn-Seq",
    "WCS",
    "WGA",
    "WGS",
    "WXS",
    # spelling differences in metadata
    'HiC-Seq',
    'ATAC-seq',
    'DNA-Seq',
    "FAIRE-Seq",
]

# %%


# %%
strategies = []
for strat in fear_counts.index:
    strategies.extend(strat.split('|'))

# %%
new_strats = []
for new in set(strategies):
    if new not in ncbi_strat:
        new_strats.append(new)

# %%
sorted(new_strats)

# %%
len(sorted(new_strats))

# %%


# %%
other_counts = fear_counts[fear_counts.index.str.contains('OTHER')]

# %%
other_counts

# %%
print(f'We also re-annotated a number of samples from OTHER to RNA-Seq (n = {other_counts["RNA-Seq|OTHER"]:,}), EST (n = {other_counts["EST|OTHER"]:,}), ChIP-Seq (n = {other_counts["ChIP-Seq|OTHER"]:,}), or WGS (n = {other_counts["WGS|OTHER"]:,}).')

# %%


# %%


# %%


# %%


# %%


# %%

