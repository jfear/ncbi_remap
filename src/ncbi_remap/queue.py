from pathlib import Path
from typing import Optional, Union, List, Tuple, Set, Callable, TYPE_CHECKING
import pickle
import re
from collections.abc import Sequence

import numpy as np
import pandas as pd
from dask.distributed import Client


class Queue:
    """A file based queue system.

    Attributes
    ----------
    srxs : list
        List of SRXs in the queue.
    srrs : list
        List of SRRs in the queue.
    sample_table : pd.DataFrame
        A pandas DataFrame of queued SRXs and SRRs
    n_srx: int
        The total number of SRXs that could be queue.
    n_srr: int
        The total number of SRRs that could be queue.

    """

    _module_path = Path(__file__).absolute().parent
    _srx_pattern = re.compile(r"^[SDE]RX\d+$")
    _srr_pattern = re.compile(r"^[SDE]RR\d+$")

    def __init__(
        self,
        targets: Optional[Union[str, List[str]]] = None,
        completed: Optional[Union[str, List[str]]] = None,
        problems: Optional[Union[str, List[str]]] = None,
        subset: Optional[Union[str, List[str]]] = None,
        srx2srr: Optional[str] = None,
        size: int = 100,
    ):
        """File based queue system.

        For the parameters [targets, completed, problems, subset] you can pass:
            * A pickle file containing list of SRXs or SRRs
            * A folder containing files named as SRXs or SRRs
            * A list of SRXs or SRRs

        Parameters
        ----------
        targets : str or list, optional
            A list of SRXs|SRRs that to queue.
        completed : str or list, optional
            A list of SRXs|SRRs that have already completed and can be skipped.
        problems : str or list, optional
            A list of SRXs|SRRs that are problematic and can be skipped.
        subset : str or list, optional
            A list of SRXs|SRRs that should be run if in targets.
        srx2srr : str, optional
            A file containing all SRX|SRR pairs.
        size : int, optional
            Number of samples to queue, by default 100

        """
        self._srxs = set()
        self._srrs = set()
        self.srx2srr = None
        self.size = size

        self._load_srx2srr(srx2srr)
        self._load_targets(targets)
        self._load_completed(completed)
        self._load_problems(problems)
        self._load_subset(subset)

    @property
    def srxs(self):
        return sorted(self._srxs, key=self.sort_accession)[: self.size]

    @property
    def srrs(self):
        return self.sample_table.srr.unique().tolist()

    @property
    def sample_table(self):
        return self.srx2srr.query(f"srx == {self.srxs} & srr == {list(self._srrs)}").sort_values(
            ["srx", "srr"]
        )

    @property
    def n_srxs(self):
        return len(self._srxs)

    @property
    def n_srrs(self):
        return len(self._srrs)

    def get_srrs(self, srx: str) -> List[str]:
        return self.sample_table.query(f"srx == '{srx}'").srr.unique().tolist()

    def get_srx(self, srr: str) -> str:
        return self.sample_table.query(f"srr == '{srr}'").srx.values[0]

    def _load_srx2srr(self, file_name):
        file_name = file_name or self._module_path / "../../output/srx2srr.csv"
        self.srx2srr = pd.read_csv(file_name)

    def _load_targets(self, targets):
        if targets is None:
            self._srxs |= set(self.srx2srr.srx.values)
            self._srrs |= set(self.srx2srr.srr.values)
            return

        values = self._get_values(targets)
        srxs, srrs = self._translate_sample(values)
        self._srxs |= srxs
        self._srrs |= srrs

    def _load_completed(self, completed):
        if completed is None:
            return

        values = self._get_values(completed)
        srxs, srrs = self._translate_sample(values)
        self._srxs -= srxs
        self._srrs -= srrs

    def _load_problems(self, problems):
        if problems is None:
            return

        values = self._get_values(problems)
        srxs, srrs = self._translate_sample(values)
        self._srxs -= srxs
        self._srrs -= srrs

    def _load_subset(self, subset):
        if subset is None:
            return

        values = self._get_values(subset)
        srxs, srrs = self._translate_sample(values)
        self._srxs &= srxs
        self._srrs &= srrs

    def _get_values(self, name) -> Set[str]:
        if isinstance(name, str):
            if re.match(self._srx_pattern, name) or re.match(self._srr_pattern, name):
                return set([name])

            if name.endswith(".pkl"):
                return set(pickle.load(open(name, "rb")))

            if Path(name).is_dir():
                return {x.stem for x in Path(name).iterdir()}

            return set([name])

        if isinstance(name, Sequence):
            return set().union(*(self._get_values(value) for value in name))

    def _translate_sample(self, values: Set[str]) -> Tuple[Set[str], Set[str]]:
        value = next(iter(values))
        if re.match(self._srx_pattern, value):
            srxs = values
            srrs = set(self.srx2srr.query(f"srx == {list(srxs)}").srr.values)
        elif re.match(self._srr_pattern, value):
            srrs = values
            srxs = set(self.srx2srr.query(f"srr == {list(srrs)}").srx.values)
        else:
            return set(), set()

        return srxs, srrs

    @staticmethod
    def sort_accession(x):
        match = re.findall(r"\d+", x)[0]
        return int(match)


def dask_run_srr_checker(
    ids: pd.DataFrame, checker: Callable, pattern: str, client: Client
) -> List[pd.DataFrame]:
    """Helper function to run function using dask."""
    futures = client.map(lambda x: checker(x["srx"], x["srr"], pattern), ids.to_dict("records"))
    return client.gather(futures)


def dask_run_srx_checker(srxs: List, checker: Callable, client: Client) -> List[pd.DataFrame]:
    """Helper function to run function using dask."""
    futures = client.map(checker, srxs)
    return client.gather(futures)


def build_parsed(fname: str, parser: Callable, **kwargs) -> pd.DataFrame:
    df = parser(fname)
    if df is None:
        return None
    return df.assign(**kwargs).set_index(sorted(kwargs.keys(), reverse=True))


def dask_run_srr_parser(
    ids: pd.DataFrame, parser: Callable, pattern: str, client: Client
) -> List[pd.DataFrame]:
    """Helper function to run function using dask."""
    futures = client.map(lambda x: build_parsed(parser, pattern, **x), ids.to_dict("records"))
    return client.gather(futures)


def dask_run_srx_parser(
    srxs: List[str], parser: Callable, pattern: str, client: Client
) -> List[pd.DataFrame]:
    """Helper function to run function using dask."""
    futures = client.map(lambda srx: build_parsed(parser, pattern, srx=srx), srxs)
    return client.gather(futures)


def check_indicator_file(srx: str, srr: str, pattern: str) -> Tuple:
    """Check if an indicator file is present."""
    fname = Path(pattern.format(srx=srx, srr=srr))
    if fname.exists():
        return srx, srr

    return np.nan, np.nan
