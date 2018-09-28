#!/usr/bin/env python
"""Normalize gene level counts using various methods.

To make thing easier down the road normalizing all of the data with CPM, TPM,
and RPKM.
"""

import pandas as pd

from ncbi_remap.logging import logger
from ncbi_remap.normalization import cpm, rpkm, tpm


def normalize(func, data, oname, *args):
    norm = func(data, *args)
    norm.to_parquet(oname)


def main():
    gene_lengths = pd.read_parquet(snakemake.input.gene_lens).gene_length
    data = pd.read_parquet(snakemake.input.gene_counts)

    logger.info('Calculating CPM')
    normalize(cpm, data, snakemake.output.cpm)

    logger.info('Calculating RPKM')
    normalize(rpkm, data, snakemake.output.rpkm, gene_lengths)

    logger.info('Calculating TPM')
    normalize(tpm, data, snakemake.output.tpm, gene_lengths)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.info('Script was interrupted.')
    except Exception:
        raise
    else:
        logger.info('Script Complete')
