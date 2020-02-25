# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.spatial.distance import mahalanobis, pdist, cdist
from scipy.stats import probplot
from scipy.special import factorial
from sklearn.covariance import MinCovDet, ShrunkCovariance

from ncbi_remap.normalization import tpm


# %%
df = pd.read_csv("../output/agg-rnaseq-wf/gene_counts.tsv", sep="\t", index_col=0)

# %%
print(df.shape)

# %%
gene_lengths = pd.read_csv("../output/gene_ts_lengths.tsv", sep="\t", index_col=0).squeeze()
norm = np.log10(tpm(df.T, gene_lengths).dropna() + 1).T

# %%
cov = ShrunkCovariance().fit(norm)

# %%
testis = open("../output/notebook/2020-02-21_testis.txt", "r").read().strip().split()

# %%
_mean =norm.mean(axis=1).values

# %%
norm.apply(lambda x: mahalanobis(x, _mean))

# %%
probplot(md, sparams=(cov.covariance_.shape[0]), dist='chi2', plot=plt)


# %%
