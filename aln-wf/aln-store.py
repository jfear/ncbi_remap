#!/usr/bin/env python
"""Update the data store.

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
        description='Alignment Queue',
        usage="aln_store.py <command> [args]")

    subparsers = parser.add_subparsers(dest='command', help='Commands')

    parser_queue = subparsers.add_parser('queue',
                                         help='Interact with aln/queue')
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
        help='Run update on aln/queue.',
    )
    parser_queue.add_argument(
        '--print',
        dest='queue_print',
        action='store_true',
        help='Print the aln/queue',
    )

    parser_data = subparsers.add_parser('data',
                                        help='Interact with aln/workflow')
    parser_data.add_argument(
        '-j',
        dest='n_workers',
        action='store',
        default=2,
        type=int,
        help='Number of workers [default 2]',
    )
    return parser.parse_args()


def get_keepers():
    return [
        patterns['hisat2']['summary'],
        patterns['srxMerge']['feature_counts']['summary'],
        patterns['srxMerge']['feature_counts_intergenic']['summary'],
        patterns['srxMerge']['samtools_stats'],
        patterns['srxMerge']['samtools_idxstats'],
        patterns['srxMerge']['bamtools_stats'],
        patterns['srxMerge']['bigWig'],
        patterns['srxMerge']['bigWigFlyBase'],
    ]


def dask_run(ids, parser, *args, **kwargs):
    """Helper function to run function using dask."""
    futures = []
    for idx, (srx, srr) in ids.iterrows():
        futures.append(delayed(parser)(srx, srr, *args, **kwargs))

    return daskClient.gather(daskClient.compute(futures))


def dask_run_large_data(name, pattern, munger, parser):
    ids = store['prealn/complete']
    futures = []
    for idx, dat in ids.groupby('srx'):
        srx = idx
        srrs = dat.srr.tolist()
        futures.append(delayed(munger)(srx, srrs, name, pattern, parser))

    return daskClient.gather(daskClient.compute(futures))


def check_indicator_file(srx, srr, pattern, *args, **kwargs):
    """Check if an indicator file is present."""
    fname = Path(pattern.format(srx=srx, srr=srr, **kwargs))
    if fname.exists():
        return srx, srr

    return np.nan, np.nan


def check_outputs(srx, srr, *args, **kwargs):
    """Check if output file is present."""
    for pattern in get_keepers():
        if 'strand' in pattern:
            for strand in ['first', 'second']:
                fname = Path(pattern.format(srx=srx, srr=srr, strand=strand,
                                            **kwargs))
                if not fname.exists():
                    return np.nan, np.nan
        else:
            fname = Path(pattern.format(srx=srx, srr=srr, **kwargs))
            if not fname.exists():
                return np.nan, np.nan
    return srx, srr


def process_outputs():
    """Check if all output files are present and update queues."""
    logger.info('Checking outputs.')
    res = dask_run(store['aln/queue'], check_outputs)
    ids = pd.DataFrame(res, columns=['srx', 'srr']).dropna()

    mask = store['aln/queue'].srr.isin(ids.srr)
    cnt = mask.sum()
    if cnt > 0:
        store.remove('aln/queue', mask)
        logger.info('{:,} SRRs removed from queue.'.format(cnt))
        store.append('aln/complete', ids, data_columns=True, format='t')
    logger.info('Complete')


def process_flags():
    """Check all flag and indicator files and store results."""
    # indicator files
    for key in ['alignment_bad', ]:

        queue = store['aln/queue']

        logger.info(f'Checking {key} files.')
        if store.__contains__(key):
            ids = queue[~queue.srr.isin(store[key].srr)]
        else:
            ids = queue
        dat = dask_run(ids, check_indicator_file, patterns[key])
        df = pd.DataFrame(dat, columns=['srx', 'srr'])\
            .dropna()\
            .reset_index(drop=True)

        # Prevent duplicate entries in the store
#         if store.__contains__(key):
#             idx = store[key].srr
#             srr_already_in_store = df.srr.isin(idx)
#             df = df[~srr_already_in_store]

        mask = queue.srr.isin(df.srr)
        cnt = mask.sum()
        # remove form queue because an indicator files indicates that
        # something is something wrong.
        if cnt > 0:
            store_key = 'aln/' + key
            store.append(store_key, df, data_columns=True, format='t')

            logger.info('{:,} SRRs removed from queue.'.format(cnt))
            store.remove('aln/queue', mask)


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
        ("queued", "aln/queue"),
        ("completed", 'aln/complete'),
        ("alignment bad", 'aln/alignment_bad'),
    ]
    parsed = [parse_pairs(k, v) for k, v in pairs]

    report = '\nCurrent Queue Summary\n'
    for k, v in parsed:
        report += '{:,}\t\t\t{}\n'.format(v, k)

    print(report)


def agg_small(name, pattern, func):
    complete = store['aln/complete']
    key = 'aln/workflow/' + name
    logger.info(f'Checking {key} files.')

    if store.__contains__(key):
        srrs_already_processed = complete.srr.isin(store[key].index.get_level_values('srr'))
        toRun = complete[~srrs_already_processed]
    else:
        toRun = complete

    dat = dask_run(toRun, func, pattern)
    dfs = [x for x in dat if x is not None]
    if len(dfs) > 0:
        df = pd.concat(dfs)
        try:
            store.append(key, df, data_columns=True, format='t')
        except ValueError:
            # Sometimes data types don't match, use pandas to coerce.
            tmp = store[key]
            store.append(key, pd.concat([tmp, df]), data_columns=True,
                         format='t', append=False)


def munge_as_series(srx, srrs, name, pattern, parser):
    srx_dirname = Path(pattern).parents[1].as_posix().format(srx=srx)
    oname = Path(srx_dirname, name)

    if oname.exists():
        return

    srs = []
    for srr in srrs:
        sr = parser(srx, srr, pattern)
        sr.index = sr.index.droplevel([0, 1])
        sr.name = srr
        srs.append(sr)

    if len(srs) > 1:
        df = pd.concat(srs, axis=1)
    else:
        df = sr.to_frame()

    df.to_parquet(oname.as_posix(), engine='pyarrow')


def munge_as_dataframe(srx, srrs, name, pattern, parser):
    srx_dirname = Path(pattern).parents[1].as_posix().format(srx=srx)
    oname = Path(srx_dirname, name)

    if oname.exists():
        return

    dfs = []
    for srr in srrs:
        df = parser(srx, srr, pattern)
        dfs.append(df)

    if len(dfs) > 1:
        df = pd.concat(dfs)

    df.to_parquet(oname.as_posix(), engine='pyarrow')


def add_small_data_to_store():
    """Parse smaller datasets and add them to the store.

    Smaller datasets are munged and added to their corresponding table in the
    data store.
    """
    from ncbi_remap.parser import (
        parse_hisat2,
        parse_samtools_stats,
        parse_bamtools_stats,
        parse_featureCounts_summary
    )

    agg_small('hisat2', patterns['hisat2']['bam'] + '.log', parse_hisat2)

    agg_small('samtools_stats', patterns['samtools_stats'],
              parse_samtools_stats)

    agg_small('bamtools_stats', patterns['bamtools_stats'],
              parse_bamtools_stats)

    agg_small('feature_counts/summary',
              patterns['feature_counts']['summary'],
              parse_featureCounts_summary)


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

    logger.info(f'Building samtools idxstats files')
    dask_run_large_data('srr_samtools_idxstats.parquet',
                        patterns['samtools_idxstats'],
                        munge_as_dataframe,
                        parse_samtools_idxstats)

    logger.info(f'Building counts files')
    dask_run_large_data('srr_counts.parquet',
                        patterns['feature_counts']['counts'],
                        munge_as_series,
                        parse_featureCounts_counts)

    logger.info(f'Building junction counts files')
    dask_run_large_data('srr_jcounts.parquet',
                        patterns['feature_counts']['jcounts'],
                        munge_as_dataframe,
                        parse_featureCounts_jcounts)


def updated_data_tables():
    logger.info('Aggregating data tables.')
    add_small_data_to_store()
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
        updated_data_tables()

    store.close()
    daskClient.close()
