#!/usr/bin/env python
"""Initialize or update the data store.

The HDF5 data store is used as a queue system as well as a data store for
various flags and datasets.

"""
import sys
from pathlib import Path
import argparse
from yaml import load

import numpy as np
import pandas as pd
from dask import delayed
from dask.distributed import Client

sys.path.insert(0, '../lib')
# from ncbi_remap.config import DATA_STORE
from ncbi_remap.logging import logger
from ncbi_remap.snakemake import get_patterns

DATA_STORE = '../sra_new.h5'


def get_keepers():
    return [
        patterns['fastq_screen'],
        patterns['layout'],
        patterns['strand'],
        patterns['hisat2']['summary'],
        patterns['feature_counts']['summary'],
        patterns['samtools_stats'],
        patterns['samtools_idxstats'],
        patterns['bamtools_stats'],
        patterns['picard']['markduplicates']['metrics'],
    ]


def get_options():
    """Function to pull in command line arguments"""
    parser = argparse.ArgumentParser(description='Update Pre-Alignment Queue')
    parser.add_argument(
        '-j',
        dest='n_workers',
        action='store',
        required=False,
        default=2,
        type=int,
        help='Number of workers [default 2]',
    )
    parser.add_argument(
        '--host',
        dest='host',
        action='store',
        required=False,
        default='localhost',
        help='Host running mongo database [Optional].'
    )
    parser.add_argument(
        '--port',
        dest='port',
        action='store',
        required=False,
        type=int,
        default=27022,
        help='Mongo database port [Optional].'
    )
    return parser.parse_args()


def get_update_db_ids():
    """Pull out SRX, SRR from mongo."""
    from pymongo import MongoClient
    mongoClient = MongoClient(host=args.host, port=args.port)
    db = mongoClient['sra']
    ncbi = db['ncbi']

    # Dump all ids out of database
    df = pd.DataFrame(list(ncbi.aggregate([
        {
            '$unwind': '$sra.run'
        },
        {
            '$match': {
                'sra.run.run_id': {'$exists': 1}
            }
        },
        {
            '$project': {
                '_id': 0,
                'srx': '$_id',
                'srr': '$sra.run.run_id'
            }
        },
    ])))

    return df[['srx', 'srr']]


def dask_run(ids, func, *args, **kwargs):
    """Helper function to run function using dask."""
    futures = []
    for idx, (srx, srr) in ids.iterrows():
        futures.append(delayed(func)(srx, srr, *args, **kwargs))

    return daskClient.gather(daskClient.compute(futures))


def check_flag_file(srx, srr, pattern, *args, **kwargs):
    """Parse flag file and return value."""
    fname = Path(pattern.format(srx=srx, srr=srr, **kwargs))
    try:
        with fname.open() as fh:
            flag = fh.read().strip()
        return srx, srr, flag
    except FileNotFoundError:
        return srx, srr, np.nan


def check_indicator_file(srx, srr, pattern, *args, **kwargs):
    """Check if an indicator file is present."""
    fname = Path(pattern.format(srx=srx, srr=srr, **kwargs))
    if fname.exists():
        return srx, srr

    return np.nan, np.nan


def initialize_store():
    """Initialize the data store with ids and prealn/queue."""
    logger.info('Initialize data store.')
    store['ids'] = get_update_db_ids()
    store.put('prealn/queue', store['ids'], data_columns=True,
              format='table')

    # Add flag files
    for key in 'layout', 'strand':
        dat = dask_run(store['prealn/queue'], check_flag_file, patterns[key])

        df = pd.DataFrame(dat, columns=['srx', 'srr', key])\
            .set_index(['srx', 'srr'])\
            .iloc[:, 0]\
            .dropna()

        store.put(key, df, data_columns=True, format='t')

    # Add indicator files
    for key in ['download_bad', 'alignment_bad', 'quality_scores_bad',
                'abi_solid']:
        dat = dask_run(store['prealn/queue'], check_indicator_file,
                       patterns[key])

        df = pd.DataFrame(dat, columns=['srx', 'srr'])\
            .dropna()\
            .reset_index(drop=True)

        store_key = 'prealn/' + key
        store.put(store_key, df, data_columns=True, format='t')


def check_outputs(srx, srr, *args, **kwargs):
    """Check if output file is present."""
    for pattern in get_keepers():
        fname = Path(pattern.format(srx=srx, srr=srr, **kwargs))
        if not fname.exists():
            return np.nan, np.nan
    return srx, srr


def process_outputs():
    """Check if all output files are present and update queues."""
    logger.info('Checking outputs.')
    queue = store['prealn/queue'].copy()
    res = dask_run(queue, check_outputs)
    ids = pd.DataFrame(res, columns=['srx', 'srr']).dropna()

    srrs = ids.srr
    store.remove('prealn/queue', 'srr == srrs')
    store.append('prealn/complete', ids, data_columns=True, format='t')
    store.append('aln/queue', ids, data_columns=True, format='t')
    logger.info('{:,} SRRs removed from queue'.format(len(srrs)))
    logger.info('Complete')


def process_flags():
    """Check all flag and indicator files and store results."""
    # flag files
    queue = store['prealn/queue'].copy()
    for key in 'layout', 'strand':
        logger.info(f'Checking {key} files.')
        dat = dask_run(queue, check_flag_file, patterns[key])
        df = pd.DataFrame(dat, columns=['srx', 'srr', key])\
            .set_index(['srx', 'srr'])\
            .iloc[:, 0]\
            .dropna()

        idx = store[key].index.get_level_values('srr')
        srr_already_in_store = df.index.get_level_values('srr').isin(idx)
        store.append(key, df[~srr_already_in_store], data_columns=True,
                     format='t')

    # indicator files
    for key in ['download_bad', 'alignment_bad', 'quality_scores_bad',
                'abi_solid']:
        logger.info(f'Checking {key} files.')
        dat = dask_run(queue, check_indicator_file, patterns[key])
        df = pd.DataFrame(dat, columns=['srx', 'srr'])\
            .dropna()\
            .reset_index(drop=True)

        srrs = df.srr
        # remove form queue because an indicator files indicates that
        # something is something wrong.
        if len(srrs) > 0:
            store_key = 'prealn/' + key
            store.append(store_key, df, data_columns=True, format='t')
            store.remove('prealn/queue', 'srr == srrs')
            break


if __name__ == '__main__':
    args = get_options()
    daskClient = Client(n_workers=args.n_workers, threads_per_worker=1)
    store = pd.HDFStore(DATA_STORE)
    patterns = get_patterns('patterns.yaml')

    # Initialize if needed.
    if not store.__contains__('ids'):
        logger.info('Initialize data store.')
        initialize_store()

    process_flags()
    process_outputs()
