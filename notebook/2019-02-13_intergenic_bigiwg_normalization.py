# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 0.8.6
#   kernelspec:
#     display_name: Python [conda env:pybigwig]
#     language: python
#     name: conda-env-pybigwig-py
# ---

# %% [markdown]
# # Intergenic BigWig Normalization

# %% [markdown]
# In the aggregated bigwigs there is a lot of noise. I want to see if we can calculate the upper quartile of the intergenic nosie and subtract that.

# %%
from pathlib import Path
from collections import namedtuple
import shutil
from hashlib import md5

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pyBigWig
from lcdblib.utils import chrom_convert

# %%
scratch = Path('/scratch/ncbi_remap')

# %%
def median_95th(bw, genic):
    results = []
    for chrom, dd in genic.groupby('chrom'): 
        if chrom not in ['X', '2L', '2R', '3L', '3R']:
            continue
        chrom_length = bw.chroms()[chrom]
        coverage = np.nan_to_num(bw.values(chrom, 0, chrom_length, numpy=True))
        for _, gene in dd.iterrows():
            coverage[gene.start:gene.end] = np.nan
        coverage = coverage[~np.isnan(coverage)]
        if coverage.shape[0] < 10:
            continue
        percentile = np.percentile(coverage, 99)
        results.append(percentile)
    return np.median(results)


Bed = namedtuple('Bed', 'chroms starts ends values')


def coverage_to_bed(chrom, length, coverage):
    chroms, starts, ends, values = [], [], [], []
    previous = None
    for loc, value in enumerate(coverage):
        if loc == 0:
            start = 0
            previous = value
            continue
            
        if value == previous and loc != length - 1:
            continue
            
        chroms.append(str(chrom))
        starts.append(int(start))
        ends.append(int(loc))
        values.append(float(previous))
        previous = value
        start = loc
        
    return  Bed(chroms, starts, ends, values)

# %%
chrom_mapper = chrom_convert.import_conversion('UCSC', 'FlyBase')

# %%
genic = (
    pd.read_csv(scratch / 'dmel_r6-11.sort.bed12', sep='\t', usecols=[0, 1, 2], header=None, names=['chrom', 'start', 'end'])
    .drop_duplicates()
    .assign(chrom = lambda df: df.chrom.map(chrom_mapper))
)

# %%
bw_plus = pyBigWig.open(str(scratch / 'fb_plus.bw'))
bw_minus = pyBigWig.open(str(scratch / 'fb_minus.bw'))

# %%
plus = median_95th(bw_plus, genic)
minus = median_95th(bw_minus, genic)

# %%
plus, minus

# %%
chrom, start, end = genic.iloc[46, :].values
coverage = np.nan_to_num(bw_plus.values(chrom, int(start), int(end), numpy=True)) - plus
coverage[coverage < 0] = 0
print(f'{chrom}:{start}..{end}')
plt.plot(coverage)
ax = plt.gca()

# %%
chrom = '3R'
start = 7924323
end = 7967408
coverage = np.nan_to_num(bw_minus.values(chrom, int(start), int(end), numpy=True)) - minus
coverage[coverage < 0] = 0
plt.plot(coverage)
ax = plt.gca()

# %%
bw_out = pyBigWig.open(str(scratch / 'fb_plus_scaled.bw'), 'w')
header = sorted([(k, v) for k, v in bw_plus.chroms().items() if k == '2L'], key=lambda x: x[0])
bw_out.addHeader(header)

for chrom, length in header:
    coverage = np.nan_to_num(bw_plus.values(chrom, 0, length, numpy=True)) - plus
    coverage[coverage < 0] = 0
    bed = coverage_to_bed(chrom, length, coverage)
    bw_out.addEntries(bed.chroms, bed.starts, ends=bed.ends, values=bed.values)
bw_out.close()

# %%
bw_plus.close()
bw_minus.close()

# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%
srxs = pd.read_csv('../output/uber-stranded-wf/uber_srxs.txt', header=None).iloc[:, 0].values

# %%



# %%
def copy_file(pth1, pth2, check=False):
    if not pth2.exists():
        shutil.copy(str(pth1), str(pth2))
    
    if check:
        with pth1.open(mode='rb') as fh1, pth2.open(mode='rb') as fh2:
            md1 = md5(fh1.read()).hexdigest()
            md2 = md5(fh2.read()).hexdigest()

        if md1 != md2:
            pth2.unlink()
            copy_file(pth1, pth2)


# %%
for srx in srxs:
    pth1 = Path(f'../output/aln-wf/samples/{srx}/{srx}.flybase.first.bw')
    pth2 = Path(f'/scratch/ncbi_remap/{srx}.flybase.first.bw')
    copy_file(pth1, pth2)
    print(srx, end=', ')

# %%
for srx in srxs:
    pth1 = Path(f'../output/aln-wf/samples/{srx}/{srx}.flybase.second.bw')
    pth2 = Path(f'/scratch/ncbi_remap/{srx}.flybase.second.bw')
    copy_file(pth1, pth2)
    print(srx, end=', ')

# %%



# %%



# %%



# %%


# %%


# %%


# %%



# %%
counts = np.zeros(length)
for i, srx in enumerate(srxs):
    bw = pyBigWig.open(f'../output/aln-wf/samples/{srx}/{srx}.flybase.first.bw')
    counts += (np.nan_to_num(bw.values(chrom, 0, length, numpy=True)) > 0).astype(int)
    bw.close()
    if i % 100 == 0:
        print(i, end=', ')

counts /= len(srxs)

# %%
sns.kdeplot(counts)

# %%



# %%
bed = coverage_to_bed(chrom, length, counts)

# %%
bw_out = pyBigWig.open(str(scratch / 'fb_plus_prop.bw'), 'w')
header = [(chrom, length)]
bw_out.addHeader(header)
bw_out.addEntries(bed.chroms, bed.starts, ends=bed.ends, values=bed.values)
bw_out.close()

# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%


