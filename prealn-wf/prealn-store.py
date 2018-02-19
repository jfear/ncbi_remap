#!/usr/bin/env python
"""Initialize or update the data store.

The HDF5 data store is used as a queue system as well as a data store for
various flags and datasets.

"""
import sys
from pathlib import Path
import argparse

import numpy as np
import pandas as pd
from dask import delayed
from dask.distributed import Client

sys.path.insert(0, '../lib')
# from ncbi_remap.config import DATA_STORE
from ncbi_remap.logging import logger
from ncbi_remap.snakemake import get_patterns

DATA_STORE = '../sra_new.h5'


def get_options():
    parser = argparse.ArgumentParser(
        description='Pre-Alignment Queue',
        usage="prealn_store.py <command> [args]")

    subparsers = parser.add_subparsers(dest='command', help='Commands')

    parser_init = subparsers.add_parser(
        'init',
        help='Initialize data store'
    )
    parser_init.add_argument(
        '--host',
        dest='host',
        action='store',
        default='localhost',
        help='Host running mongo database [Optional].'
    )
    parser_init.add_argument(
        '--port',
        dest='port',
        action='store',
        type=int,
        default=27022,
        help='Mongo database port [Optional].'
    )
    parser_init.add_argument(
        '-j',
        dest='n_workers',
        action='store',
        default=2,
        type=int,
        help='Number of workers [default 2]',
    )
    parser_init.add_argument(
        '--append',
        dest='append',
        action='store_true',
        help='Download and append new ids to the data store.',
    )

    parser_queue = subparsers.add_parser('queue',
                                         help='Interact with prealn/queue')
    parser_queue.add_argument(
        '-j',
        dest='n_workers',
        action='store',
        default=2,
        type=int,
        help='Number of workers [default 2]',
    )
    parser_queue.add_argument(
        '--update',
        dest='queue_update',
        action='store_true',
        help='Run update on prealn/queue.',
    )
    parser_queue.add_argument(
        '--print',
        dest='queue_print',
        action='store_true',
        help='Print the prealn/queue',
    )
    parser_queue.add_argument(
        '--abi',
        dest='queue_abi',
        action='store_true',
        help='Check for ABI Solid data by fastq.',
    )

    parser_queue = subparsers.add_parser('data',
                                         help='Interact with prealn/workflow')
    parser_queue.add_argument(
        '-j',
        dest='n_workers',
        action='store',
        default=2,
        type=int,
        help='Number of workers [default 2]',
    )
    parser_queue.add_argument(
        '--update',
        dest='data_update',
        action='store_true',
        help='Run update on prealn/queue.',
    )
    return parser.parse_args()


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


def check_abi_solid(srx, srr, *args, **kwargs):
    """Check Fastq is Abi Solid."""
    logger.info('Checking queue for Abi Solid')
    import gzip
    layout = Path(patterns['layout'].format(srx=srx, srr=srr))
    abi_solid = Path(patterns['abi_solid'].format(srx=srx, srr=srr))
    fname = Path(patterns['fastq']['r1'].format(srx=srx, srr=srr))
    r2 = Path(patterns['fastq']['r2'].format(srx=srx, srr=srr))

    if abi_solid.exists():
        return

    try:
        with layout.open() as fh:
            flag = fh.read().strip()

        if flag == 'keep_R2':
            fname = r2
        with gzip.open(fname, 'rt') as fh:
            _, l = fh.readline(), fh.readline()
            if l.startswith('T'):
                abi_solid.touch()
    except (KeyError, FileNotFoundError):
        pass

    logger.info('Abi Solid Complete')


def initialize_store():
    """Initialize the data store with ids and prealn/queue."""
    logger.info('Initialize data store.')
    logger.info('Querying database for ids.')
    store['ids'] = get_update_db_ids()
    logger.info('Adding ids to store.')
    store.put('prealn/queue', store['ids'], data_columns=True, format='table')
    logger.info('Initialization complete')


def append_store():
    """Download new data and append to store."""
    logger.info('Append new ids to data store.')
    curr_ids = store['ids'].copy()

    logger.info('Querying database for ids.')
    ids = get_update_db_ids()
    store['ids'] = ids

    # Sometimes an id is removed from the SRA, need to make sure to remove it
    # from the queue.
    srr_no_longer_in_store = curr_ids.srr.isin(ids.srr)
    removed_ids = curr_ids[~srr_no_longer_in_store]

    if len(removed_ids) > 0:
        logger.info(
            'There are {:,} ids no longer in the SRA.'.format(len(removed_ids))
        )
        srrs = removed_ids.srr
        store.remove('prealn/queue', 'srr == srrs')

    # Do not add ids already in the store back to the queue.
    srr_already_in_store = ids.srr.isin(curr_ids.srr)
    new_ids = ids[~srr_already_in_store]

    if len(new_ids) == 0:
        logger.info('There are no new ids.')
        return

    logger.info('There are {:,} new ids.'.format(len(new_ids)))
    logger.info('Adding ids to store.')
    store.append('prealn/queue', new_ids, data_columns=True, format='table')
    logger.info('Initialization complete')


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
    res = dask_run(store['prealn/queue'], check_outputs)
    ids = pd.DataFrame(res, columns=['srx', 'srr']).dropna()

    mask = store['prealn/queue'].srr.isin(ids.srr)
    cnt = mask.sum()
    if cnt > 0:
        store.remove('prealn/queue', mask)
        logger.info('{:,} SRRs removed from queue.'.format(cnt))
        store.append('prealn/complete', ids, data_columns=True, format='t')
        store.append('aln/queue', ids, data_columns=True, format='t')
    logger.info('Complete')


def process_flags():
    """Check all flag and indicator files and store results."""
    # flag files
    for key in 'layout', 'strand':
        logger.info(f'Checking {key} files.')
        dat = dask_run(store['prealn/queue'], check_flag_file, patterns[key])
        df = pd.DataFrame(dat, columns=['srx', 'srr', key])\
            .set_index(['srx', 'srr'])\
            .iloc[:, 0]\
            .dropna()

        # Prevent duplicate entries in the store
        if store.__contains__(key):
            idx = store[key].index.get_level_values('srr')
            srr_already_in_store = df.index.get_level_values('srr').isin(idx)
            df = df[~srr_already_in_store]

        store.append(key, df, data_columns=True, format='t')

    # indicator files
    for key in ['download_bad', 'alignment_bad', 'quality_scores_bad',
                'abi_solid']:
        logger.info(f'Checking {key} files.')
        dat = dask_run(store['prealn/queue'], check_indicator_file,
                       patterns[key])
        df = pd.DataFrame(dat, columns=['srx', 'srr'])\
            .dropna()\
            .reset_index(drop=True)

        # Prevent duplicate entries in the store
        if store.__contains__(key):
            idx = store[key].srr
            srr_already_in_store = df.srr.isin(idx)
            df = df[~srr_already_in_store]

        mask = store['prealn/queue'].srr.isin(df.srr)
        cnt = mask.sum()
        # remove form queue because an indicator files indicates that
        # something is something wrong.
        if cnt > 0:
            store_key = 'prealn/' + key
            store.append(store_key, df, data_columns=True, format='t')

            logger.info('{:,} SRRs removed from queue.'.format(cnt))
            store.remove('prealn/queue', mask)


def update_queue():
    logger.info('Updating queue.')
    process_flags()
    process_outputs()
    logger.info('Update complete')


def print_queue():
    def parse_pairs(name, key):
        try:
            return name, store[key].shape[0]
        except KeyError:
            return name, 0

    pairs = [
        ("ids in the system", 'ids'),
        ("queued", "prealn/queue"),
        ("completed", 'prealn/complete'),
        ("download bad", 'prealn/download_bad'),
        ("quality scores bad", 'prealn/quality_scores_bad'),
        ("alignment bad", 'prealn/alignment_bad'),
        ("abi solid", 'prealn/abi_solid'),
    ]
    parsed = [parse_pairs(k, v) for k, v in pairs]

    if store.__contains__('layout'):
        layout = store['layout'].value_counts()
        parsed.extend([
            ("Single End", layout.SE),
            ("Pair End", layout.PE),
            ("Really Single End R1", layout.keep_R1),
            ("Really Single End R2", layout.keep_R2),
        ])

    if store.__contains__('strand'):
        strand = store['strand'].value_counts()

        parsed.extend([
            ("first strand", strand.same_strand),
            ("second strand", strand.opposite_strand),
            ("unstranded", strand.unstranded),
        ])

    report = '\nCurrent Queue Summary\n'
    for k, v in parsed:
        report += '{:,}\t\t\t{}\n'.format(v, k)

    print(report)


def agg_small(name, pattern, func):
    logger.info(f'Checking {key} files.')
    complete = store['prealn/complete']
    dat = dask_run(complete, func, pattern)
    df = pd.concat([x for x in dat if x is not None])

    key = 'prealn/workflow/' + name
    idx = store[key].index.get_level_values('srr')
    srr_already_in_store = df.index.get_level_values('srr').isin(idx)
#     store.append(key, df[~srr_already_in_store], data_columns=True, format='t')


def store_small_data():
    """Parse smaller datasets and add them to the store.

    Smaller datasets are munged and added to their corresponding table in the
    data store.
    """
    from ncbi_remap.parser import (
        parse_fastq_summary,
        parse_fastq_screen,
        parse_hisat2,
        parse_samtools_stats,
        parse_bamtools_stats,
        parse_picard_markduplicate_metrics,
        parse_picardCollect_summary,
        parse_featureCounts_summary
    )

    agg_small('fastq', patterns['fastq']['summary'], parse_fastq_summary)
    agg_small('fastq_screen', patterns['fastq_screen'], parse_fastq_screen)
    agg_small('hisat2', patterns['hisat2']['bam'] + '.log', parse_hisat2)
    agg_small('samtools_stats', patterns['samtools_stats'],
              parse_samtools_stats)
    agg_small('bamtools_stats', patterns['bamtools_stats'],
              parse_bamtools_stats)
    agg_small('markduplicates',
              patterns['picard']['markduplicates']['metrics'],
              parse_picard_markduplicate_metrics)
    agg_small('collectrnaseqmetrics/first',
              patterns['picard']['collectrnaseqmetrics']['metrics']['first'],
              parse_picardCollect_summary)
    agg_small('collectrnaseqmetrics/second',
              patterns['picard']['collectrnaseqmetrics']['metrics']['second'],
              parse_picardCollect_summary)
    agg_small(
        'collectrnaseqmetrics/unstranded',
        patterns['picard']['collectrnaseqmetrics']['metrics']['unstranded'],
        parse_picardCollect_summary
    )
    agg_small('feature_counts/summary',
              patterns['feature_counts']['summary'],
              parse_featureCounts_summary)


def agg_large(pattern, func):
    pass


def munge_big_data():
    """Munge larger datasets for easier handling with dask.

    Larger tables are not performant when added to the datastore. Dask does
    well accessing multiple files on disk in parallel. Dask documentation
    suggests to use the apache parquet. Here we munge the data and save it out
    to the parquet format.
    """
    from ncbi_remap.parser import (
        parse_samtools_idxstats,
        parse_featureCounts_counts,
        parse_featureCounts_jcounts,
    )

    agg_large(patterns['samtools_idxstats'], parse_samtools_idxstats)
    agg_large(patterns['feature_counts']['counts'], parse_featureCounts_counts)
    agg_large(patterns['feature_counts']['jcounts'],
              parse_featureCounts_jcounts)


def agg_data():
    logger.info('Aggregating data tables.')
    store_small_data()
    munge_big_data()
    logger.info('Aggregation complete.')


if __name__ == '__main__':
    args = get_options()
    daskClient = Client(n_workers=args.n_workers, threads_per_worker=1)
    store = pd.HDFStore(DATA_STORE)
    patterns = get_patterns('patterns.yaml')

    if args.command == 'init':
        if args.append:
            append_store()
        else:
            initialize_store()
    elif args.command == 'queue':
        if args.queue_update:
            update_queue()
        elif args.queue_print:
            print_queue()
        elif args.queue_abi:
            dask_run(store['prealn/queue'], check_abi_solid)
    elif args.command == 'data':
        agg_data()

    store.close()
    daskClient.close()