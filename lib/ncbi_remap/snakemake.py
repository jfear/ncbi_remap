#!/usr/bin/env python
"""Set of helpers for use with snakemake."""
import os
from pathlib import Path
from itertools import zip_longest

import pandas as pd
from .io import add_id, remove_id
from .logging import logger

def wrapper_for(path):
    return 'file:' + path


def put_flag(fname, flag):
    """Export flag from file.

    This is a little helper script to write a flag to a file.

    """
    with open(fname, 'w') as fh:
        fh.write(flag)


def get_flag(fname):
    """Import flag from file.

    This is a little helper script to import a flag stored in a file.

    """
    with open(fname) as fh:
        return fh.read().strip()


def check_download(store, pattern, **kwargs):
    """Checks FASTQ logs determine if bad DB entry.

    There are some SRRs that are no longer accessible in the SRA database. I
    get an error in the fastq-dump log that I want to translate into the
    DOWNLOAD_BAD flag.

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    patter : str
        File naming pattern for the ALIGNEMNT_BAD file.
    **kwargs
        Keywords needed to fill the pattern.

    """
    logName = Path(pattern.format(**kwargs))
    flagBad = Path(logName.parent, 'DOWNLOAD_BAD')
    if logName.exists():
        dat = logName.read_text()
        if 'failed to resolve accession' in dat:
            flagBad.touch()

    if flagBad.exists():
        remove_id(store, 'prealn/queue', **kwargs)
        add_id(store, 'prealn/download_bad', **kwargs)
        return True


def check_alignment(store, pattern, **kwargs):
    """Checks for ALIGNMENT_BAD file.

    If there is an ALIGNMENT_BAD file then remove from the queue, add to
    complete, and add to 'prealn/alignment_bad'.

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    patter : str
        File naming pattern for the ALIGNEMNT_BAD file.
    **kwargs
        Keywords needed to fill the pattern.

    """
    ab = pattern.format(**kwargs)
    if os.path.exists(ab):
        remove_id(store, 'prealn/queue', **kwargs)
        add_id(store, 'prealn/alignment_bad', **kwargs)
        return True


def check_layout(store, pattern, **kwargs):
    """Checks LAYOUT file.

    Parses layout file and adds to the corresponding hdf5 lists.

        * 'layout/SE'
        * 'layout/PE'
        * 'layout/keep_R1'
        * 'layout/keep_R2'

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    patter : str
        File naming pattern for the ALIGNEMNT_BAD file.
    **kwargs
        Keywords needed to fill the pattern.

    """
    with open(pattern.format(**kwargs)) as fh:
        strand = fh.read().strip()
        if strand == 'SE':
            key = 'layout/SE'
        elif strand == 'PE':
            key = 'layout/PE'
        elif strand == 'keep_R1':
            key = 'layout/keep_R1'
        elif strand == 'keep_R2':
            key = 'layout/keep_R2'
        add_id(store, key, **kwargs)


def check_strand(store, pattern, **kwargs):
    """Checks STRAND file.

    Parses strand file and adds to the corresponding hdf5 lists.

        * 'strand/first'
        * 'strand/second'
        * 'strand/unstranded'

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    patter : str
        File naming pattern for the ALIGNEMNT_BAD file.
    **kwargs
        Keywords needed to fill the pattern.

    """
    with open(pattern.format(**kwargs)) as fh:
        strand = fh.read().strip()
        if (strand == 'first_strand') | (strand == 'same_strand'):
            key = 'strand/first'
        elif (strand == 'second_strand') | (strand == 'opposite_strand'):
            key = 'strand/second'
        elif strand == 'unstranded':
            key = 'strand/unstranded'
        add_id(store, key, **kwargs)


def combine(func, pattern, row):
    """Parse input file using parser function.

    Uses a parser function to return a data frame.

    Parameters
    ----------
    func : .parser.parser_*
        A parser function that returns a dataframe.
    pattern : str
        A file name pattern that can be filled with row.
    row : pd.Series
        A row from a sample table.

    Returns
    -------
    pd.DataFrame
        parse dataframe.

    """
    return func(row.srx, row.srr, pattern.format(**row.to_dict()))


def agg(store, key, func, pattern, df, large=False):
    """Aggregator to import tables and dump into a hdf5 store.

    Parameters
    ----------
    store : pd.HDFStore
        Data store to save results.
    key : str
        Node in the data store to save results.
    func : .parser.parser_*
        A parser function that returns a dataframe.
    pattern : str
        A file name pattern that can be filled with row.
    df : pd.DataFrame
        A sample table containing samples to parse.
    large : bool
        If True import datasets one at time because they are large.

    Returns
    -------
    None

    """
    done = []
    if store.get_node(key):
        done = set()        # change to set to keep memory down.
        # Iterate over chunks and grab srrs that are already there.
        for chunk in store.select(key, chunksize=1e5):
            idx = chunk.index.names.index('srr')
            done |= set(chunk.index.levels[idx])
        done = list(done)

    dfs = []
    for _, row in df.iterrows():
        if row.srr in done:
            continue
        try:
            if large:
                dd = combine(func, pattern, row)
                store.append(key, dd)
            else:
                dat = combine(func, pattern, row)
                if dat is None:
                    continue
                dfs.append(dat)
        except ValueError:
            logger.error('Error parsing {}->{}'.format(row.srx, row.srr))

    if dfs:
        ddf = pd.concat(dfs)
        store.append(key, ddf)


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks.

    This function is from the itertools recipes section.
    https://docs.python.org/3/library/itertools.html#itertools-recipes

    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

