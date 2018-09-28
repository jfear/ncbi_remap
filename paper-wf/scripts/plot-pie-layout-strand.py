"""Plot a set of pie charts to show layout and strandedness."""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.style.use('scripts/paper_1c.mplstyle')

store = pd.HDFStore(snakemake.input[0], mode='r')
complete = store['aln/complete'].srx.unique().tolist()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(4, 2), gridspec_kw=dict(wspace=.1))

# Layout
summary = store.select('layout', where='srx == complete').value_counts()
total = summary.sum()

mapper = dict(SE='SE', PE='PE', keep_R1='PE-R1', keep_R2='PE-R2')
summary.index = summary.index.map(mapper)

summary.plot.pie(autopct=lambda pct: '{:,.0f}'.format(np.round(pct * total / 100, 0)),
                 textprops=dict(fontsize=10), ax=ax1)

ax1.set_aspect('equal')
ax1.set_ylabel('')
ax1.set_title('Library Layout')

# Strandedness
summary = store.select('strand', where='srx == complete').value_counts()
total = summary.sum()

mapper = dict(unstranded='US', same_strand='FS', opposite_strand='RS')
summary.index = summary.index.map(mapper)

summary.plot.pie(autopct=lambda pct: '{:,.0f}'.format(np.round(pct * total / 100, 0)),
                 textprops=dict(fontsize=10), ax=ax2)

ax2.set_aspect('equal')
ax2.set_ylabel('')
ax2.set_title('Strandedness')


fig.subplots_adjust(left=-0.01, right=.98, bottom=-0.01, top=.9)
fig.savefig(snakemake.output[0])
