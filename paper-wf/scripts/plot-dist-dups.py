import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns

mpl.style.use('scripts/paper_1c.mplstyle')

store = pd.HDFStore(snakemake.input[0], mode='r')

fig, ax = plt.subplots(figsize=(4, 2))
dups = store['prealn/workflow/markduplicates'].PERCENT_DUPLICATION * 100
sns.distplot(dups, ax=ax)
ax.set_xlabel('Percent Duplication (%)')
ax.set_ylabel('Density')
ax.set_title('Distribution of Duplication Rate')

fig.subplots_adjust(left=.15, right=.95, bottom=.25, top=.85)
fig.savefig(snakemake.output[0])
