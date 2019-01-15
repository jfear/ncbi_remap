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
# # Look at pybigwig

# %% [markdown]
# I am prototyping different methods for normalizing merged coverage counts. Here I am exploring what pybigwig can do and if I should be using it.

# %%
import pyBigWig
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

# %%
# open the bigwig file
plus = pyBigWig.open('../output/uber-stranded-wf/samples/merged.plus.bw')
minus = pyBigWig.open('../output/uber-stranded-wf/samples/merged.minus.bw')

# %%
# get a dictionary of chromsizes
chrom_sizes = plus.chroms()

# %%
# Get header info
#bw.header()


# grab a range from a chromosome
#vals = bw.values("chrM", 1000, 1050, numpy=True)
#vals

# calculate mean
#np.round(bw.stats("chrM", 1000, 1050), 0)[0], np.round(np.mean(vals), 0)

# calculate min
#bw.stats("chrM", 1000, 1050, type='min')[0], np.min(vals)

# calculate max
#bw.stats("chrM", 1000, 1050, type='max')[0], np.max(vals)

# calculate the fraction of reads covered
#bw.stats("chrM", 1000, 1050, type='coverage')[0], np.sum(vals > 0) / len(vals)

# calculate the std, not appears to be slightly different than numpy's std
#np.round(bw.stats("chrM", 1000, 1050, type='std')[0], 0), np.round(np.std(vals), 0)

# %%
# grab a range from a chromosome
pvals = np.nan_to_num(plus.values("chrY", 1, chrom_sizes['chrY'], numpy=True))
mvals = np.nan_to_num(minus.values("chrY", 1, chrom_sizes['chrY'], numpy=True))
_max = max(np.max(pvals), np.max(mvals))
#_max = max(plus.header()['maxVal'], minus.header()['maxVal'])

# %%
pvals_scaled = pvals / _max
mvals_scaled = mvals / _max

# %%


# %%
start = np.argmax(mvals) - 1_000
end = np.argmax(mvals) + 1_000

# %%
fig, ax = plt.subplots(1, 1, figsize=(30, 20))
plt.plot(mvals_scaled)
plt.plot(pvals_scaled)
ax.set_xlim(start, end)

# %%


# %%


# %%
plus_out = pyBigWig.open('../output/uber-stranded-wf/samples/merged.plus.scaled.bw', 'w')
minus_out = pyBigWig.open('../output/uber-stranded-wf/samples/merged.minus.scaled.bw', 'wb')

# %%
header = [(k, v) for k, v in plus.chroms().items()]
plus_out.addHeader(header)
minus_out.addHeader(header)

# %%
for chrom, size in header:
    plus_intervals = np.array(plus.intervals(chrom))
    pstarts = plus_intervals[:, 0].astype(np.int64)
    pends = plus_intervals[:, 1].astype(np.int64)
    pvals = plus_intervals[:, 2]
    pchroms = np.array([chrom] * pvals.shape[0])
    
    minus_intervals = np.array(minus.intervals(chrom))
    mstarts = minus_intervals[:, 0].astype(np.int64)
    mends = minus_intervals[:, 1].astype(np.int64)
    mvals = minus_intervals[:, 2]
    mchroms = np.array([chrom] * mvals.shape[0])
    
    _max = max(np.max(pvals), np.max(mvals))
    
    pvals = pvals / _max
    mvals = mvals / _max
    
    plus_out.addEntries(pchroms.tolist(), pstarts.tolist(), ends=pends.tolist(), values=pvals.tolist())
    minus_out.addEntries(mchroms.tolist(), mstarts.tolist(), ends=mends.tolist(), values=mvals.tolist())
    
    break
    
plus_out.close()
minus_out.close()

# %%


# %%

