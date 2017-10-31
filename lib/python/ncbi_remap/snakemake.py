#!/usr/bin/env python
"""Set of helpers for use with snakemake."""
import sys

import pandas as pd
from dask import delayed, compute

from .io import add_table

def wrapper_for(path):
    return 'file:' + path


def put_flag(fname, flag):
    with open(fname, 'w') as fh:
        fh.write(flag)


def get_flag(fname):
    with open(fname) as fh:
        return fh.read().strip()


def combine(func, pattern, row):
    df = func(row.srr, pattern.format(**row.to_dict()))
    df.index.name = 'srr'
    df['srx'] = row.srx
    return df.reset_index().set_index(['srx', 'srr']).reset_index()


def agg(store, key, func, pattern, df):
    if store.get_node(key):
        done = store[key].srr.tolist()
    else:
        done = []

    dfs = []
    for i, row in df.iterrows():
        if row.srr in done:
            continue
        dfs.append(delayed(combine)(func, pattern, row))

    if dfs:
        ddf = pd.concat(compute(*dfs), ignore_index=True)
        try:
            add_table(store, key, data=ddf, columns=['srx', 'srr'])
        except ValueError as e:
            print(ddf.head(), ddf.tail())
            raise e

