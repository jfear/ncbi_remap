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
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# %%
# !ls /home/fearjm/scratch/gffcompare/

# %%
DIR = Path('~/scratch/gffcompare').expanduser()

# %%
df = pd.read_csv(DIR / 'gffcmp_pacbio.w1118_testi1.collapsed.gff.tmap', sep='\t')

# %%
codes = ['=', 'c']
num_matched = df.query(f'class_code == {codes}').qry_gene_id.unique().shape[0]
print(f'{num_matched:,}')

# %%
codes = ['u',]
num_unique = df.query(f'class_code == {codes}').qry_gene_id.unique().shape[0]
print(f'{num_unique:,}')

# %%
print('\n'.join(df.query(f'class_code == {codes}').groupby('qry_gene_id').first().query('num_exons == "2"').index.tolist()))

# %%

# %%
codes = ['m',]
num_retained = df.query(f'class_code == {codes}').qry_gene_id.unique().shape[0]
print(f'{num_retained:,}')

# %%
print('\n'.join(df.query(f'class_code == {codes}').groupby('qry_gene_id').first().query('num_exons >= 2').index.tolist()))

# %%

# %%
codes = ['o',]
num_other = df.query(f'class_code == {codes}').qry_gene_id.unique().shape[0]
print(f'{num_other:,}')

# %%
print('\n'.join(df.query(f'class_code == {codes}').groupby('qry_gene_id').first().query('num_exons >= 2').index.tolist()))

# %%
