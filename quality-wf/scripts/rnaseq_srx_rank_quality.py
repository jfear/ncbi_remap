"""Create a ranking of SRXs (RNA-Seq) by quality.

I want to rank samples by their "quality" using the various metrics created by
the prealn-wf and aln-wf. The actual rankings are not important, but give an
easy way to find the top fraction of samples. This is going to be particularly
useful for building the aggregated sets.

"""
import pandas as pd

# Connect to data store
store = pd.HDFStore('../output/sra.h5', mode='r')

# I only want to focus on the RNA-Seq samples
libstrat = pd.read_parquet(snakemake.input[0]).iloc[:, 0]
rnaseq = libstrat[libstrat == 'RNA-Seq'].index.tolist()

# Get a list of SRXs that have completed the entire workflow.
srxs = store.select('aln/complete', where='srx == rnaseq').srx.unique().tolist()

# Build Features
# Hisat2
cols = [
    'per_alignment'
]
hisat2 = store.select('aln/workflow/hisat2', where='srx == srxs', columns=cols)

# CollectRNASeqMetrics
cols = [
    'PCT_CODING_BASES',
    'PCT_INTERGENIC_BASES',
]
metrics = store.select('prealn/workflow/collectrnaseqmetrics/unstranded', where='srx == srxs', columns=cols)
metrics['INV_PCT_INTERGENIC_BASES'] = 1 - metrics.PCT_INTERGENIC_BASES
metrics.drop('PCT_INTERGENIC_BASES', axis=1, inplace=True)

# Markduplicates
cols = [
    'PERCENT_DUPLICATION',
]
dups = store.select('prealn/workflow/markduplicates', where='srx == srxs', columns=cols)
dups['INV_PERCENT_DUPLICATION'] = 1 - dups.PERCENT_DUPLICATION
dups.drop('PERCENT_DUPLICATION', axis=1, inplace=True)

# Munge features together
features = hisat2.join(metrics).join(dups)
ranks = features.rank().median(axis=1).groupby('srx').mean().rank(ascending=False).sort_values().to_frame()
ranks.columns = ['quality_rank']
ranks.index.name = 'srx'
ranks.to_csv(snakemake.output[0], sep='\t')
