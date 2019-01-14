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
# # Check modENCODE in my GEO submission

# %% [markdown]
# Claire had sent me this email:
#
# > I noticed that modEncode datasets are not included in GSM3273571.  Are they submitted separated? 
#
# I want to verify that the RNA-Seq is indeed in this submission.

# %%
import pandas as pd

# %%
modencode = pd.read_csv('../output/modENCODE_sampletable.tsv', sep='\t', index_col=[0, 1])
modencode.index = modencode.index.droplevel('srr')
modencode.drop_duplicates(inplace=True)

# %%
print(modencode.shape)
modencode.head(20)

# %%
my_geo = pd.read_csv('../output/geo-wf/rnaseq_metadata.tsv', sep='\t', index_col=0)
my_geo_list = my_geo.index.tolist()
print(len(my_geo_list))

# %%
in_mine = modencode[modencode.index.isin(my_geo_list)].copy()

# %%
not_in_mine = modencode[~modencode.index.isin(my_geo_list)].copy()

# %%
in_mine.shape, not_in_mine.shape

# %%


# %%
modencode['fear_submission'] = False

# %%
modencode.loc[modencode.index.isin(my_geo_list), 'fear_submission'] = True

# %%
modencode.to_csv('../output/notebook/2019-01-14_check_modENCODE_in_submission.csv')

# %%

