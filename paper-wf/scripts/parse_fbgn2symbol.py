"""Create a dictionary mapping fbgn to gene symbol."""

from pathlib import Path
from pickle import dump

import pandas as pd

from ncbi_remap.config import REFERENCES_DIR, config

oname = snakemake.output[0]
assembly = config['assembly']
tag = config['tag']
annot_fn = Path(REFERENCES_DIR, assembly, tag, 'fb_annotation', f'{assembly}_{tag}.fb_annotation')

df = pd.read_csv(annot_fn, sep='\t', index_col=1).fillna('nan')
fbgn2symbol = df['gene_symbol'].to_dict()

with open(oname, 'wb') as fo:
    dump(fbgn2symbol, fo)
