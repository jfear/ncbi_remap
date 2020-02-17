#!/usr/bin/env python
"""Interact with the alignment workflow data store.

The HDF5 data store is used as a queue system as well as a data store for
various flags and datasets generated by the alignment workflow. This script
provides easy access to the alignment workflow queue as well as easy updating
of the queue and data. It also parses large files and saves them as a
data frame. Specifically this script can:

* Print a summary of the queue `./aln-store.py queue --print`
* Scan the `samples` directory and remove completed items from the queue
  `./aln-store.py queue --update -j 8`
* Scan the `samples` directory and parse data and add small datasets and flags
  to the data sore and save large datasets as a data frame `./aln-store.py data
  -j 8`

Flags checked:

* alignment_bad
* bigwig_bad

Small datasets added to the data store:

* `samtools stats`
* `bamtools stats`
* Genic `feature counts` summary
* Intergenic `feature counts` summary

Large datasets parsed to data frame and saved:

* `samtools idxstats`
* Genic `feature counts` counts
* Genic `feature counts` junction counts
* Intergenic `feature counts` counts
* Intergenic `feature counts` junction counts
"""

from pathlib import Path
import argparse

import numpy as np
import pandas as pd
from dask import delayed
from dask.distributed import Client

from ncbi_remap.logging import logger
from ncbi_remap.snakemake import get_patterns
from ncbi_remap.queue import (
    dask_run_srr_checker, dask_run_srx_checker, dask_run_srx_parser,
    dask_run_srr_parser, check_indicator_file,
)

WORKDIR = Path(__file__).absolute().parent
DATA_STORE = Path(WORKDIR, '../output/sra.h5').as_posix()


def get_options():
    parser = argparse.ArgumentParser(
        description=__doc__,
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
        patterns['srxMerge']['feature_counts']['summary'],
        patterns['srxMerge']['feature_counts_intergenic']['summary'],
        patterns['srxMerge']['samtools_stats'],
        patterns['srxMerge']['samtools_idxstats'],
        patterns['srxMerge']['bamtools_stats'],
        patterns['srxMerge']['bigWig'],
        patterns['srxMerge']['bigWigFlyBase'],
    ]


def dask_run_large_data(name, pattern, munger, parser):
    srxs = store['aln/complete'].srx.unique()
    futures = []
    for srx in srxs:
        futures.append(delayed(munger)(srx, name, pattern, parser))

    return daskClient.gather(daskClient.compute(futures))


def check_outputs(srx):
    """Check if output file is present."""
    for pattern in get_keepers():
        if 'strand' in pattern:
            for strand in ['first', 'second']:
                fname = Path(pattern.format(srx=srx, strand=strand))
                if not fname.exists():
                    return np.nan
        else:
            fname = Path(pattern.format(srx=srx))
            if not fname.exists():
                return np.nan

    return srx


def process_outputs():
    """Check if all output files are present and update queues."""
    logger.info('Checking outputs.')
    dat = dask_run_srx_checker(store['aln/queue'].srx.unique(), check_outputs,
                               daskClient)
    df = pd.DataFrame(dat, columns=['srx']).dropna()

    mask = store['aln/queue'].srx.isin(df.srx)
    cnt = mask.sum()
    if cnt > 0:
        ids = store['aln/queue'][mask]
        store.remove('aln/queue', mask)
        logger.info('{:,} SRRs removed from queue.'.format(cnt))
        store.append('aln/complete', ids, data_columns=True, format='t')
    logger.info('Complete')


def process_flags():
    """Check all flag and indicator files and store results."""
    # indicator files
    for key in ['atropos_bad', 'alignment_bad', 'bigwig_bad']:
        queue = store['aln/queue']
        logger.info(f'Checking {key} files.')
        if store.__contains__(key):
            ids = queue[~queue.srr.isin(store[key].srr)]
        else:
            ids = queue
        dat = dask_run_srr_checker(ids, check_indicator_file, patterns[key],
                                   daskClient)

        df = pd.DataFrame(dat, columns=['srx', 'srr'])\
            .dropna()\
            .reset_index(drop=True)

        # remove form queue because an indicator files indicates that
        # something is something wrong. Using mask because querying a list of
        # ids was not working.
        mask = queue.srr.isin(df.srr)
        cnt = mask.sum()
        if cnt > 0:
            store_key = 'aln/' + key
            store.append(store_key, df, data_columns=True, format='t',
                         min_itemsize={'srx': 10, 'srr': 10})

            logger.info('{:,} SRRs removed from queue.'.format(cnt))
            store.remove('aln/queue', mask)


def munge_as_dataframe(srx, name, pattern, parser):
    try:
        oname = Path(name.format(srx=srx))

        if oname.exists():
            return

        fname = pattern.format(srx=srx)
        df = parser(fname)
        df = df.assign(srx=srx)

        oname.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(oname.as_posix(), engine='pyarrow')
    except Exception as e:
        logger.error(f'There was a problem parsing {fname}')
        raise e


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
    dask_run_large_data(Path(WORKDIR, '../output/aln-wf/samtools_idxstats/{srx}.parquet').as_posix(),
                        patterns['srxMerge']['samtools_idxstats'],
                        munge_as_dataframe,
                        parse_samtools_idxstats)

    logger.info(f'Building counts files')
    dask_run_large_data(Path(WORKDIR, '../output/aln-wf/gene_counts/{srx}.parquet').as_posix(),
                        patterns['srxMerge']['feature_counts']['counts'],
                        munge_as_dataframe,
                        parse_featureCounts_counts)

    logger.info(f'Building junction counts files')
    dask_run_large_data(Path(WORKDIR, '../output/aln-wf/junction_counts/{srx}.parquet').as_posix(),
                        patterns['srxMerge']['feature_counts']['jcounts'],
                        munge_as_dataframe,
                        parse_featureCounts_jcounts)

    logger.info(f'Building intergenic counts files')
    dask_run_large_data(
        Path(WORKDIR, '../output/aln-wf/intergenic_counts/{srx}.parquet').as_posix(),
        patterns['srxMerge']['feature_counts_intergenic']['counts'],
        munge_as_dataframe,
        parse_featureCounts_counts
    )


def agg_small(name, pattern, parser):
    complete = store['aln/complete']
    key = 'aln/workflow/' + name
    logger.info(f'Checking {key} files.')

    if store.__contains__(key):
        srrs_already_processed = complete.srr.isin(
            store[key].index.get_level_values('srr')
        )
        toRun = complete[~srrs_already_processed]
    else:
        toRun = complete

    dat = dask_run_srr_parser(toRun, parser, pattern, daskClient)
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


def agg_srx_small(name, pattern, parser):
    complete = store['aln/complete'].srx.unique()
    key = 'aln/workflow/' + name
    logger.info(f'Checking {key} files.')

    if store.__contains__(key):
        srx_already_processed = store[key].index.get_level_values('srx').tolist()
        toRun = [x for x in complete if x not in srx_already_processed]
    else:
        toRun = complete

    dat = dask_run_srx_parser(toRun, parser, pattern, daskClient)
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


def add_small_data_to_store():
    """Parse smaller datasets and add them to the store.

    Smaller datasets are munged and added to their corresponding table in the
    data store.
    """
    from ncbi_remap.parser import (
        parse_hisat2,
        parse_samtools_stats,
        parse_bamtools_stats,
        parse_featureCounts_summary,
    )

    agg_small('hisat2', patterns['hisat2']['bam'] + '.log', parse_hisat2)

    agg_srx_small('samtools_stats', patterns['srxMerge']['samtools_stats'],
                  parse_samtools_stats)

    agg_srx_small('bamtools_stats', patterns['srxMerge']['bamtools_stats'],
                  parse_bamtools_stats)

    agg_srx_small('feature_counts/summary',
                  patterns['srxMerge']['feature_counts']['summary'],
                  parse_featureCounts_summary)

    agg_srx_small('feature_counts_intergenic/summary',
                  patterns['srxMerge']['feature_counts_intergenic']['summary'],
                  parse_featureCounts_summary)


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
        ("atropos bad", 'aln/atropos_bad'),
        ("alignment bad", 'aln/alignment_bad'),
        ("bigwig bad", 'aln/bigwig_bad'),
    ]
    parsed = [parse_pairs(k, v) for k, v in pairs]

    report = '\nCurrent Queue Summary\n'
    for k, v in parsed:
        report += '{:,}\t\t\t{}\n'.format(v, k)

    print(report)


def update_queue():
    logger.info('Updating queue.')
    process_flags()
    process_outputs()
    logger.info('Update complete')


def updated_data_tables():
    logger.info('Aggregating data tables.')
    add_small_data_to_store()
    munge_big_data()
    logger.info('Aggregation complete.')


if __name__ == '__main__':
    args = get_options()
    daskClient = Client(n_workers=args.n_workers, threads_per_worker=1)
    store = pd.HDFStore(DATA_STORE)
    patterns = get_patterns(Path(WORKDIR, 'patterns.yaml').as_posix())

    if args.command == 'queue':
        if args.queue_update:
            update_queue()
        elif args.queue_print:
            print_queue()
    elif args.command == 'data':
        updated_data_tables()

    store.close()
    daskClient.close()