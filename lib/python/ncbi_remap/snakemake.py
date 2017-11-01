#!/usr/bin/env python
"""Set of helpers for use with snakemake."""
import sys
from itertools import zip_longest

import pandas as pd
from dask import delayed, compute

from .io import add_table

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
    df = func(row.srr, pattern.format(**row.to_dict()))
    df.index.name = 'srr'
    df['srx'] = row.srx
    return df.reset_index().set_index(['srx', 'srr']).reset_index()


def agg(store, key, func, pattern, df):
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

    Returns
    -------
    None

    """
    if store.get_node(key):
        done = store.select_column(key, 'srr').unique().tolist()
    else:
        done = []

    dfs = []
    for i, row in df.iterrows():
        if row.srr in done:
            continue
        dfs.append(delayed(combine)(func, pattern, row))

    if dfs:
        ddf = pd.concat(compute(*dfs), ignore_index=True)
        add_table(store, key, data=ddf, columns=['srx', 'srr'])


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks.

    This function is from the itertools recipes section.
    https://docs.python.org/3/library/itertools.html#itertools-recipes

    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def agg_big(store, key, func, pattern, df, n=100):
    """Aggregator for larger tables that won't all fit in memory.

    Some of the tables that I am trying to aggregate have thousands of rows. If
    I try to import too many at one time I run out of memory. This adds a
    grouper step that only imports 100 tables at 1 time.

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
    n : int
        Size of group to import at one time.

    Returns
    -------
    None

    """

    if store.get_node(key):
        done = store.select_column(key, 'srr').unique().tolist()
    else:
        done = []

    dfs = []
    for i, row in df.iterrows():
        if row.srr in done:
            continue

        dfs.append(delayed(combine)(func, pattern, row))

    if dfs:
        for dd in grouper(dfs, n):
            ddf = pd.concat(compute(*dd), ignore_index=True)
            add_table(store, key, data=ddf, columns=['srx', 'srr'])
