from pathlib import Path
from typing import Callable, List, Tuple, Optional
import pickle
import re

import numpy as np
from pandas import DataFrame
from dask.distributed import Client


def get_samples(
    sample_file: str, done_folder: str, subset_file: Optional[str] = None, size: int = 100
) -> List[str]:
    """Create of list of samples that need to be run.

    Parameters
    ----------
    sample_file : str
        A pickle file with a set of samples.
    done_folder : str
        A folder where all filenames are indicators if a sample has completed
        a workflow.
    subset_file : str, optional
        An optional pickle file with a subset of samples to focus on.
    size : int, optional
        The number of samples to return, by default 100

    Returns
    -------
    A list of samples.

    """
    queued = pickle.load(open(sample_file, "rb"))  # type: set

    pth = Path(done_folder)
    if pth.exists():
        completed = {x.stem for x in pth.iterdir()} # type: set
        queued = queued - completed

    if subset_file:
        focus = pickle.load(open(subset_file, "rb"))  # type: set
        queued = queued.intersection(focus)

    return sorted(list(queued), key=sort_accession)[:size]


def sort_accession(x):
    match = re.findall(r"\d+", x)[0]
    return int(match)


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
