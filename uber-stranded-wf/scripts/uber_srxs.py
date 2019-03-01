"""Generate a list of high quality rnaseq stranded libraries.

"""
import numpy as np
import pandas as pd

strand_cutoff = float(snakemake.params['strand_cutoff'])
quality_cutoff = int(snakemake.params['quality_cutoff'])

# Connect to data store
store = pd.HDFStore('../output/sra.h5', mode='r')

# Get extremely well stranded libraries
well_stranded = store.select(
    'prealn/workflow/collectrnaseqmetrics/second',
    where='PCT_CORRECT_STRAND_READS >= strand_cutoff',
    columns=['PCT_CORRECT_STRAND_READS']
).index.get_level_values('srx').unique().tolist()

# Get list of high confidence RNA-Seq libraries
library_strategy = pd.read_parquet(snakemake.input['libstrat'])

rnaseq = library_strategy.reset_index()\
        .melt(id_vars='srx')\
        .groupby('srx').value\
        .value_counts()\
        .unstack()['RNA-Seq']

rnaseq = rnaseq[rnaseq == 20]

# Get list of stranded libraries that are RNA-Seq
stranded_rnaseq = []
for srx in well_stranded:
    if srx in rnaseq.index:
        stranded_rnaseq.append(srx)

# Grab only the best quality samples
quality_ranks = pd.read_csv(snakemake.input['quality'], sep='\t', index_col=0)
stranded_quality = quality_ranks.reindex(stranded_rnaseq)
stranded_quality.sort_values(by='quality_rank', inplace=True)
cutoff = np.percentile(stranded_quality.quality_rank, quality_cutoff)
uber_srx = stranded_quality[stranded_quality.quality_rank < cutoff].index.unique().tolist()

with open(snakemake.output[0], 'w') as fh:
    fh.write('\n'.join(uber_srx))
