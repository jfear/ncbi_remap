"""Create a feature set for machine learning.

Identify and munge a set of features from the files generated by the prealn-wf
and aln-wf.

Features include:
* CollectRNASeqMetrics
    * 'PCT_CODING_BASES',
    * 'PCT_UTR_BASES',
    * 'PCT_INTRONIC_BASES',
    * 'PCT_INTERGENIC_BASES',
    * 'PCT_MRNA_BASES',
    * 'MEDIAN_CV_COVERAGE',
    * 'MEDIAN_5PRIME_BIAS',
    * 'MEDIAN_3PRIME_BIAS'
* CollectRNASeqMetrics Gene Body Coverage
* Markduplicates
    * 'PERCENT_DUPLICATION',
* FeatureCounts Summary
    * 'Assigned',
    * 'Unassigned_Ambiguity',
    * 'Unassigned_MultiMapping',
    * 'Unassigned_NoFeatures',
    * 'Unassigned_Unmapped'
* FeatureCounts Coverage Counts
    * All Genes
    * All junctions in genes
    * All Intergenic regions
"""
import numpy as np
import pandas as pd

# Connect to data store
store = pd.HDFStore('../sra.h5', mode='r')

# Generate a list of completed SRXs
srxs = store['aln/complete'].srx.unique().tolist()

# CollectRNASeqMetrics
cols = [
    'PCT_CODING_BASES',
    'PCT_UTR_BASES',
    'PCT_INTRONIC_BASES',
    'PCT_INTERGENIC_BASES',
    'PCT_MRNA_BASES',
    'MEDIAN_CV_COVERAGE',
    'MEDIAN_5PRIME_BIAS',
    'MEDIAN_3PRIME_BIAS'
]
cm = store.select('prealn/workflow/collectrnaseqmetrics/unstranded', where='srx == srxs', columns=cols)

# CollectRNASeqMetrics Gene Body Coverage
gb = store.select('prealn/workflow/collectrnaseqmetrics/genebody', where='srx == srxs')

# Markduplicates
cols = [
    'PERCENT_DUPLICATION',
]
mark = store.select('prealn/workflow/markduplicates', where='srx == srxs', columns=cols)

# FeatureCounts Summary
cols = [
    'Assigned',
    'Unassigned_Ambiguity',
    'Unassigned_MultiMapping',
    'Unassigned_NoFeatures',
    'Unassigned_Unmapped'
]
feature_summary = store.select('prealn/workflow/feature_counts/summary', columns=cols, where='srx == srxs')

# Munge Together by SRR
dat_by_srr = cm.join(gb).join(mark).join(feature_summary)
dat_by_srr_no_na = dat_by_srr.fillna(dat_by_srr.mean().to_dict(), axis=0)
dat_by_srr_no_na.shape

# FeatureCount Coverage Counts
genic = pd.read_parquet('../output/aln-downstream-wf/aggregate_genic_counts.parquet').reindex(srxs)['count']
genic.name = 'genic'
genic = genic.astype(np.float64)

intergenic = pd.read_parquet('../output/aln-downstream-wf/aggregate_intergenic_counts.parquet').reindex(srxs)['count']
intergenic.name = 'intergenic'
intergenic = intergenic.astype(np.float64)

junctions = pd.read_parquet('../output/aln-downstream-wf/aggregate_junction_counts.parquet').reindex(srxs)['count']
junctions.name = 'junctions'
junctions = junctions.astype(np.float64)

coverage = pd.concat([genic, intergenic, junctions], axis=1)

# Make final feature set aggregated to SRX
features = dat_by_srr_no_na.reset_index().groupby('srx').mean().join(coverage)

features.to_parquet(snakemake.output[0])
