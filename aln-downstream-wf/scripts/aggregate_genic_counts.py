"""Aggregate gene level counts to single matrix.

Currently each gene level count has it's own parquest file. This is makes it
really slow to work with. Here I parse all of these files and combine them into
a single parquet file. A nice feature of parquet is that it is a column score,
so you can easily import a specific set of columns using:

`pd.read_parquet('../../output/aln-downstream-wf/gene_counts_wide.parquet', columns=['SRX00000', 'SRX00001'])`

This process takes a long time (~1h) and uses a lot of RAM. I process all of
the files in batches (chunks) and save the intermediate files to a cache. The
cache defaults `~/.cache`, but you can set `CACHE_DIR` as an environmental
variable to change this.

"""

from pathlib import Path

import pandas as pd

from ncbi_remap.logging import logger
from ncbi_remap.basic import grouper
from ncbi_remap.config import CACHE_DIR

CACHE = Path(CACHE_DIR, 'gene_counts_wide_cache')
CACHE.mkdir(parents=True, exist_ok=True)
logger.info(f'Saving intermediate files to: {CACHE.as_posix()}')


def clean_cache():
    if CACHE.exists():
        logger.info(f'Removing cache: {CACHE.as_posix()}')
        for fname in CACHE.iterdir():
            fname.unlink()

        CACHE.rmdir()


def run_group(i, group):
    oname = Path(CACHE, f'chunk_{i}.parquet')

    if oname.exists():
        logger.info(f'Skipping chunk {i}')
        return

    logger.info(f'Running chunk {i}')
    dfs = []
    for fname in group:
        if fname is None:
            # this is fillvalue from grouper
            continue

        _df = pd.read_parquet(fname).reset_index().pivot(index='FBgn', columns='srx', values='count')
        dfs.append(_df)

    df = pd.concat(dfs, axis=1)
    df.to_parquet(oname)


def merge_chunks():
    chunks = CACHE.glob('*.parquet')

    dfs = []
    for chunk in chunks:
        dfs.append(pd.read_parquet(chunk))

    return pd.concat(dfs, axis=1)


def main():
    files = list(Path('../output/aln-wf/gene_counts').glob('*.parquet'))

    logger.info(f'Parsing {len(files):,} files in chunks of 3,000')
    groups = grouper(files, 3000, fillvalue=None)
    for i, group in enumerate(groups):
        run_group(i + 1, group)

    logger.info('Merging chunks')
    df = merge_chunks()

    logger.info('Writing giant table to parquet.')
    df.to_parquet(snakemake.output[0])


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.info('Script was interrupted')
    except Exception:
        logger.info('Script Failed')
        raise
    else:
        logger.info('Script Complete')
        clean_cache()
