import csv
import re
from multiprocessing.pool import ThreadPool as Pool
from pathlib import Path
from typing import Callable, Generator, Iterable, Tuple, Optional

from more_itertools import grouper
import pandas as pd


def file_name_to_sample_name(
    file_names: Iterable, pattern: re.Pattern
) -> Generator[Tuple[str, str], None, None]:
    """Uses regex to pull out sample name form file name"""
    for file_name in file_names:
        sample_name = re.findall(pattern, file_name)[0]
        yield file_name, sample_name


class RnaSeqAggregator:
    """Aggregate feature counts and appends to a data store"""

    sample_pattern = re.compile(r"[SDE]R[XR]\d+")

    def __init__(
        self,
        count_files: Iterable,
        data_store: str,
        parser: Callable,
        sample_type: str = "srx",
        threads: int = 1,
    ):
        self.count_files = count_files
        self.data_store = data_store
        self.parser = parser
        self.sample_type = sample_type
        self.pool = Pool(threads)
        self.header = None

        self.pattern = self._get_sample_pattern(self.sample_type)
        self.complete = self._load_already_processed(self.data_store)

        self.aggregate()

    @staticmethod
    def _get_sample_pattern(sample_type) -> re.Pattern:
        if sample_type == "srx":
            return re.compile(r"[SDE]RX\d+")

        if sample_type == "srr":
            return re.compile(r"[SDE]RR\d+")

        raise ValueError(f"sample_type: {sample_type} is not known.")

    @staticmethod
    def _load_already_processed(file_name) -> set:
        """Gets a set of samples already in the data store"""
        if Path(file_name).exists():
            with open(file_name, "r") as handler:
                reader = csv.reader(handler, delimiter="\t")
                return {row[0] for row in reader}

        return set()

    def _filter(self) -> Generator[Tuple[str, str], None, None]:
        """Remove samples already in the data store"""
        files_w_samples = file_name_to_sample_name(self.count_files, self.pattern)
        if len(self.complete) == 0:
            yield from files_w_samples
        else:
            for file_name, sample_name in files_w_samples:
                if not sample_name in self.complete:
                    yield file_name, sample_name

    def _load_data(self, args) -> pd.DataFrame:
        """Reads CSV using multiple threads."""
        file_name, sample_name = args
        if file_name is None:
            return
        return self.parser(file_name, sample_name, label=self.sample_type).reindex(
            columns=self.header
        )

    def aggregate(self):
        """Aggregate samples and appends to data store"""
        files_w_samples = self._filter()

        # Pull out header from data store
        if not Path(self.data_store).exists():
            # Create the data store first if it does not exist.
            file_name, sample_name = next(files_w_samples)
            df = self.parser(file_name, sample_name)
            df.to_csv(self.data_store, sep="\t")
            self.header = df.columns
        else:
            self.header = pd.read_csv(self.data_store, sep="\t", index_col=0, nrows=1).columns

        # Split into groups of 1,000 samples to keep memory use down.
        for group in grouper(files_w_samples, n=1_000, fillvalue=(None, None)):
            df = pd.concat(self.pool.map(self._load_data, group))
            df.to_csv(self.data_store, header=False, mode="a", sep="\t")
