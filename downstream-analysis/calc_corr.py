#!/usr/bin/env python
"""Calculate the pairwise correlations among samples."""
import os
import sys
from pathlib import Path
from json import dumps

import pandas as pd
from scipy.stats import spearmanr

from dask import delayed
from dask.distributed import Client, LocalCluster, as_completed

sys.path.insert(0, '../lib')
from ncbi_remap.logging import logger


DEBUG = False
DEBUG_N = 1000
cache = Path('output/corr.cache')


def cleanup():
    fnames = [
        Path('output/dask_info.json'),
        Path('output/gene_counts_corr.parquet'),
    ]
    for fname in fnames:
        if fname.exists():
            fname.unlink()


def start_cluster():
    # Start up mini cluster
    cpus = int(os.environ['SLURM_CPUS_PER_TASK'])
    mem = int(os.environ['SLURM_MEM_PER_NODE']) / 1024 / cpus
    logger.info(f'Starting Mini Cluster: Cpu={cpus} Mem={mem}')

    cluster = LocalCluster(n_workers=cpus, threads_per_worker=1,
                           memory_limit=f'{mem}GB')

    client = Client(cluster)
    with open('output/dask_info.json', 'w') as fh:
        fh.write(dumps(client._scheduler_identity, indent=True))

    return client


def get_samples():
    # Connect to data store
    store = pd.HDFStore('../sra.h5', mode='r')
    samples = store['aln/complete'].srx.unique().tolist()
    store.close()
    logger.info(f'{len(samples):,} samples loaded')

    if DEBUG:
        logger.info(f'DEBUG ON: only running {DEBUG_N:,} samples')
        return samples[:DEBUG_N]

    return samples


@delayed
def build_pairwise(srx1, samples):
    work = []
    for srx2 in samples:
        work.append((srx1, srx2))
    return work


@delayed
def calc_corr(srx1, srx2):
    if srx1 == srx2:
        return srx1, srx2, 1.0

    pth1 = Path(f'../aln-wf/output/gene_counts/{srx1}.parquet')
    pth2 = Path(f'../aln-wf/output/gene_counts/{srx2}.parquet')

    df1 = pd.read_parquet(pth1).set_index('srx', append=True)
    df2 = pd.read_parquet(pth2).set_index('srx', append=True)

    df = pd.concat([df1, df2]).unstack()
    _corr, _ = spearmanr(df)
    return srx1, srx2, _corr


def save_data():
    logger.info('Munging into dataframe')
    df = pd.read_csv(cache, sep='\t', index_col=[0, 1])
    df = df.sort_index(level=['srx1', 'srx2']).unstack().copy()
    df.columns = df.columns.droplevel(0)

    logger.info('Saving results')
    df.to_parquet('output/gene_counts_corr.parquet')
    cache.unlink()


def use_corr():
    text = input("Do you want to use corr.cache? [y|n]")
    if text == 'y':
        save_data()
        sys.exit(0)
    elif text == 'n':
        return
    else:
        use_corr()


def main():
    samples = get_samples()
    client = start_cluster()

    try:
        logger.info('Building pairwise list')
        work = []
        for i, srx1 in enumerate(samples):
            work.append(build_pairwise(srx1, samples[i:]))

        futures = client.compute(work)

        logger.info('Submitting corrrelations to Mini Cluster')
        work = []
        for fut in as_completed(futures):
            for srx1, srx2 in fut.result():
                work.append(calc_corr(srx1, srx2))

        futures = client.compute(work)

        logger.info('Writing correlation list')
        with cache.open(mode='w') as fo:
            fo.write('srx1\tsrx2\tcorrelation\n')
            for fut in as_completed(futures):
                srx1, srx2, _corr = fut.result()
                out = f'{srx1}\t{srx2}\t{_corr}\n'
                if srx1 != srx2:
                    out += f'{srx2}\t{srx1}\t{_corr}\n'
                fo.write(out)

    except KeyboardInterrupt as e:
        if cache.exists():
            cache.unlink()
        raise e
    finally:
        logger.info('Shutting down Mini Cluster')
        client.close()

    save_data()


if __name__ == '__main__':
    Path('output').mkdir(exist_ok=True)
    cleanup()

    if cache.exists():
        use_corr()

    try:
        main()
    except KeyboardInterrupt:
        logger.info('Interrupted, cleaning up files')
        cleanup()
    else:
        logger.info('Script Complete')
