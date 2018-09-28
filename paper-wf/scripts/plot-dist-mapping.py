import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns

mpl.style.use('scripts/paper_1c.mplstyle')

store = pd.HDFStore(snakemake.input[0], mode='r')

fig, ax = plt.subplots(figsize=(4, 2))
pct_map = store['prealn/workflow/hisat2'].per_alignment
sns.distplot(pct_map, ax=ax)
ax.set_xlabel('Percent Mapped (%)')
ax.set_ylabel('Density')
ax.set_title('Distribution of Percent Reads Mapping')

fig.subplots_adjust(left=.15, right=.95, bottom=.25, top=.85)
fig.savefig(snakemake.output[0])
