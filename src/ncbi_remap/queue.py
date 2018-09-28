from pathlib import Path

import numpy as np
from dask import delayed


def dask_run_srr_checker(ids, checker, pattern, client):
    """Helper function to run function using dask."""
    lazy = []
    for idx, (srx, srr) in ids.iterrows():
        lazy.append(delayed(checker)(srx, srr, pattern))
    return client.gather(client.compute(lazy))


def dask_run_srx_checker(srxs, checker, client):
    """Helper function to run function using dask."""
    lazy = []
    for srx in srxs:
        lazy.append(delayed(checker)(srx))
    return client.gather(client.compute(lazy))


@delayed
def build_parsed(fname, parser, **kwargs):
    df = parser(fname)
    if df is None:
        return None
    return df.assign(**kwargs).set_index(sorted(kwargs.keys(), reverse=True))


def dask_run_srr_parser(ids, parser, pattern, client):
    """Helper function to run function using dask."""
    lazy_dataframes = []
    for idx, (srx, srr) in ids.iterrows():
        fname = pattern.format(srx=srx, srr=srr)
        df = build_parsed(fname, parser, srx=srx, srr=srr)
        lazy_dataframes.append(df)

    return client.gather(client.compute(lazy_dataframes))


def dask_run_srx_parser(srxs, parser, pattern, client):
    """Helper function to run function using dask."""
    lazy_dataframes = []
    for srx in srxs:
        fname = pattern.format(srx=srx)
        df = build_parsed(fname, parser, srx=srx)
        lazy_dataframes.append(df)

    return client.gather(client.compute(lazy_dataframes))


def check_indicator_file(srx, srr, pattern):
    """Check if an indicator file is present."""
    fname = Path(pattern.format(srx=srx, srr=srr))
    if fname.exists():
        return srx, srr

    return np.nan, np.nan


