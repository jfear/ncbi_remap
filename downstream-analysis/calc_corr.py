#!/usr/bin/env python
"""Calculate the pairwise correlations among samples."""
import os
import sys
from pathlib import Path
from json import dumps
from itertools import product

import pandas as pd
from scipy.stats import spearmanr

from dask import delayed
from dask.distributed import Client, LocalCluster, as_completed
from dask import dataframe as dd

sys.path.insert(0, '../lib')
from ncbi_remap.logging import logger


DEBUG = True
DEBUG_N = 8


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
    return samples[::-1]


@delayed
def calc_corr(pairs):
    res = []
    for srx1, srx2 in pairs:
        if srx1 == srx2:
            res.append((srx1, srx2, 1.0))
            continue

        pth1 = Path(f'../aln-wf/output/gene_counts/{srx1}.parquet')
        pth2 = Path(f'../aln-wf/output/gene_counts/{srx2}.parquet')

        df1 = pd.read_parquet(pth1).set_index('srx', append=True)
        df2 = pd.read_parquet(pth2).set_index('srx', append=True)

        df = pd.concat([df1, df2]).unstack()
        _corr, _ = spearmanr(df)
        res.append((srx1, srx2, _corr))
        res.append((srx2, srx1, _corr))
    return res


@delayed
def save_data(srx, _corr):
    df = pd.DataFrame(_corr, columns=['srx1', 'srx2', '_corr'])
    df.set_index(['srx1', 'srx2'], inplace=True)
    df.to_parquet(f'output/samples/{srx}_corr.parquet')
    return None


def main():
    samples = get_samples()
    srxs = samples

    if DEBUG:
        logger.info(f'DEBUG ON: only running {DEBUG_N:,} samples')
        srxs = samples[:DEBUG_N]

    client = start_cluster()
    try:
        logger.info('Building Job List')
        work = []
        for i, srx in enumerate(srxs):
            if Path(f'output/samples/{srx}_corr.parquet').exists():
                continue
            pairs = product([srx,], samples[i:])
            _corr = calc_corr(pairs)
            work.append(save_data(srx, _corr))

        logger.info('Submitting jobs to Mini Cluster')
        futures = client.compute(work)
        _ = client.gather(futures)
    except KeyboardInterrupt as e:
        raise e
    finally:
        logger.info('Shutting down Mini Cluster')
        client.close()

    logger.info('Munging into dataframe')
    df = dd.read_parquet('output/samples/*_corr.parquet').compute()
    df = df.sort_index(level=['srx1', 'srx2']).unstack().copy()
    df.columns = df.columns.droplevel(0)

    logger.info('Saving results')
    df.to_parquet('output/gene_counts_corr.parquet')


if __name__ == '__main__':
    try:
        Path('output/samples').mkdir(parents=True, exist_ok=True)
        main()
    except KeyboardInterrupt:
        logger.info('Interrupted')
    else:
        logger.info('Script Complete')
