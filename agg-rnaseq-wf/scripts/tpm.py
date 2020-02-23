import os
import pandas as pd

from ncbi_remap.normalization import tpm

def main():
    gene_lengths = pd.read_csv(snakemake.input.gene_lengths, sep="\t", index_col=0).squeeze()
    header = None

    for i, df in enumerate(pd.read_csv(snakemake.input.counts, sep="\t", chunksize=2_000, index_col=0)):
        _tpm = tpm(df.T, gene_lengths).dropna().T
        if i == 0:
            header = _tpm.columns
            _tpm.to_csv(snakemake.output[0], sep="\t")
        else:
            _tpm.reindex(header).to_csv(snakemake.output[0], sep="\t", header=False, mode="a")



if __name__ == '__main__':
    if os.getenv('SNAKE_DEBUG', False):
        from ncbi_remap.debug import snakemake_debug
        snakemake = snakemake_debug(
            input=dict(
                counts="../../output/agg-rnaseq-wf/gene_counts.tsv",
                gene_lengths="../../output/gene_ts_lengths.tsv",
            )
        )

    main()
