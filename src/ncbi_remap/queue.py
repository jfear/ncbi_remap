import pickle
import re
from collections.abc import Sequence
from pathlib import Path
from textwrap import dedent
from typing import TYPE_CHECKING, Callable, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd


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
        size: Union[int, str] = 100,
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
        size : int or 'none', optional
            Number of samples to queue. If 'none' is provided then all
            samples will be used. [default 100]

        """
        self._srxs = set()
        self._srrs = set()
        self._srxs_short = None
        self._srrs_short = None
        self._sample_table = None
        self.srx2srr = None
        self.size = size

        self._load_srx2srr(srx2srr)
        self._load_targets(targets)
        self._load_completed(completed)
        self._load_problems(problems)
        self._load_subset(subset)
        self.update_filter_set()

    @property
    def srxs(self):
        return self._srxs_short

    @property
    def srrs(self):
        return self._srrs_short

    @property
    def sample_table(self):
        return self._sample_table

    @property
    def n_srxs(self):
        return len(self._srxs)

    @property
    def n_srrs(self):
        return len(self._srrs)

    def get_srrs(self, srx: str) -> List[str]:
        """Get a list of SRRs from the provided SRX."""
        return self._sample_table.query(f"srx == '{srx}'").srr.unique().tolist()

    def get_srx(self, srr: str) -> str:
        """Get an SRX from the given SRR."""
        return self._sample_table.query(f"srr == '{srr}'").srx.values[0]

    def update_filter_set(self, size=None):
        """Updates filter set using new size.

        Parameters
        ----------
        size : int or 'none', optional
            If size is provided then it updates the filtered feature set to
            include this many SRXs. If 'none is provided then all SRXs will
            be included.
        """
        if isinstance(size, int) | (size == "none"):
            self.size = size

        if self.size == "none":
            self._srxs_short = sorted(self._srxs, key=self.sort_accession)
        else:
            self._srxs_short = sorted(self._srxs, key=self.sort_accession)[: self.size]

        self._srrs_short = (
            self.srx2srr.query(f"srx == {self._srxs_short} & srr == {list(self._srrs)}")
            .srr.unique()
            .tolist()
        )
        self._sample_table = self.srx2srr.query(f"srr == {self._srrs_short}").sort_values(
            ["srx", "srr"]
        )

    @staticmethod
    def sort_accession(x):
        match = re.findall(r"\d+", x)[0]
        return int(match)

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

            if name.endswith(".txt"):
                return set(open(name, "r").read().strip().split())

            if Path(name).is_dir():
                return {x.stem for x in Path(name).iterdir()}

            return set()

        if isinstance(name, Sequence):
            return set().union(*(self._get_values(value) for value in name))

    def _translate_sample(self, values: Set[str]) -> Tuple[Set[str], Set[str]]:
        if len(values) > 0:
            value = next(iter(values))
            if re.match(self._srx_pattern, value):
                # values are SRXs
                return values, set(self.srx2srr.query(f"srx == {list(values)}").srr.values)

            if re.match(self._srr_pattern, value):
                # values are SRRs
                return set(self.srx2srr.query(f"srr == {list(values)}").srx.values), values

        return set(), set()

    def __str__(self):
        return dedent(
            f"""
            Total SRXs:\t{self.n_srxs:,}
            Total SRRs:\t{self.n_srrs:,}
            Queued SRXs:\t{len(self._srxs_short):,}
            Queued SRRs:\t{len(self._srrs_short):,}
            """
        )

    def __repr__(self):
        return (
            f"SRXs: ({self.n_srxs:,}, {len(self._srxs_short):,}), "
            f"SRRs: ({self.n_srrs:,}, {len(self._srrs_short):,}), "
        )

