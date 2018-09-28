"""Creates a mapping of gene symbol to FBgn."""

from pickle import load, dump

fname = snakemake.input[0]
oname = snakemake.output[0]

with open(fname, 'rb') as fh:
    fbgn2symbol = load(fh)

symbol2fbgn = {}
for k, v in fbgn2symbol.items():
    symbol2fbgn[v] = k

with open(oname, 'wb') as fo:
    dump(symbol2fbgn, fo)
