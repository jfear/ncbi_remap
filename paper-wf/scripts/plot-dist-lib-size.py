import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns

mpl.style.use('scripts/paper_1c.mplstyle')

store = pd.HDFStore(snakemake.input[0], mode='r')

fig, ax = plt.subplots(figsize=(4, 2))
libsize = np.log10(store['prealn/workflow/fastq'][['libsize_R1', 'libsize_R2']].max(axis=1))
sns.distplot(libsize, ax=ax)
ax.set_xlabel('Library Size (Log10)')
ax.set_ylabel('Density')
ax.set_title('Distribution of Library Sizes')

fig.subplots_adjust(left=.15, right=.95, bottom=.25, top=.85)
fig.savefig(snakemake.output[0])
