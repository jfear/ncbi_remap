from pathlib import Path
from typing import Callable, List, Tuple

import numpy as np
from pandas import DataFrame
from dask.distributed import Client


def dask_run_srr_checker(
    ids: DataFrame, checker: Callable, pattern: str, client: Client
) -> List[DataFrame]:
    """Helper function to run function using dask."""
    futures = client.map(lambda x: checker(x["srx"], x["srr"], pattern), ids.to_dict("records"))
    return client.gather(futures)


def dask_run_srx_checker(srxs: List, checker: Callable, client: Client) -> List[DataFrame]:
    """Helper function to run function using dask."""
    futures = client.map(checker, srxs)
    return client.gather(futures)


def build_parsed(fname: str, parser: Callable, **kwargs) -> DataFrame:
    df = parser(fname)
    if df is None:
        return None
    return df.assign(**kwargs).set_index(sorted(kwargs.keys(), reverse=True))


def dask_run_srr_parser(
    ids: DataFrame, parser: Callable, pattern: str, client: Client
) -> List[DataFrame]:
    """Helper function to run function using dask."""
    futures = client.map(lambda x: build_parsed(parser, pattern, **x), ids.to_dict("records"))
    return client.gather(futures)


def dask_run_srx_parser(
    srxs: List[str], parser: Callable, pattern: str, client: Client
) -> List[DataFrame]:
    """Helper function to run function using dask."""
    futures = client.map(lambda srx: build_parsed(parser, pattern, srx=srx), srxs)
    return client.gather(futures)


def check_indicator_file(srx: str, srr: str, pattern: str) -> Tuple:
    """Check if an indicator file is present."""
    fname = Path(pattern.format(srx=srx, srr=srr))
    if fname.exists():
        return srx, srr

    return np.nan, np.nan
