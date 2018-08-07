#!/usr/bin/env python
"""Get gene lengths.

Parse a single feature counts file to get the gene length info. I need gene
lengths for normalizing.
"""

from pathlib import Path

import pandas as pd

aln_wf = next(Path('../aln-wf/output/samples').iterdir())
fname = next(aln_wf.glob('*.bam.counts'))
df = pd.read_csv(fname, sep='\t', comment='#', usecols=['Geneid', 'Length']).set_index('Geneid')
df.index.name = 'FBgn'
df.columns = ['gene_length']

df.to_parquet(snakemake.output[0])
