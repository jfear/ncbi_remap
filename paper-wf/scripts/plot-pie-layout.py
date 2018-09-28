"""Plot a set of pie charts to show layout and strandedness."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import common

store = pd.HDFStore(snakemake.input[0], mode='r')
complete = store['aln/complete'].srx.unique().tolist()

fig, ax = plt.subplots(1, 1, figsize=(2, 2))

# Layout
summary = store.select('layout', where='srx == complete').value_counts()
total = summary.sum()

mapper = dict(SE='SE', PE='PE', keep_R1='PE-R1', keep_R2='PE-R2')
summary.index = summary.index.map(mapper)

summary.plot.pie(autopct=lambda pct: '{:,.0f}'.format(np.round(pct * total / 100, 0)),
                 textprops=dict(fontsize=10), ax=ax)

ax.set_aspect('equal')
ax.set_ylabel('')

fig.suptitle('Library Layout')
fig.subplots_adjust(left=0, right=.84, bottom=0, top=.98)
fig.savefig(snakemake.output[0])
