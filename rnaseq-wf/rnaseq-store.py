#!/usr/bin/env python
"""Aggregate rnaseq workflow data"""
from pathlib import Path
from typing import Optional, Set, Tuple

import pandas as pd

RNASEQ_OUTPUT = Path(__file__).absolute().parents[1] / "output/rnaseq-wf"
OUTPUTS = [
    "aln_stats",
    "atropos",
    "hisat2",
    "merge_summary",
]


def main():
    workflow_samples = load_complete_samples()
    for output in OUTPUTS:
        print(f"Aggregating: {output:>20}", end="\t")
        aggregate_data_store(
            workflow_samples, RNASEQ_OUTPUT / output, RNASEQ_OUTPUT / f"{output}.parquet"
        )


def load_complete_samples() -> Set[str]:
    with (RNASEQ_OUTPUT / "done.txt").open() as fh:
        srxs = [srx.strip() for srx in fh]
        srrs = pd.read_csv(RNASEQ_OUTPUT / "../srx2srr.csv").query("srx == @srxs").srr.to_list()
    return set(srxs).union(set(srrs))


def aggregate_data_store(workflow_samples: set, data_folder_pth: Path, data_store_pth: Path):
    old_samples = load_data_store(data_store_pth)
    new_samples = find_new_samples(workflow_samples, old_samples, data_folder_pth)

    print(f"({len(new_samples):,})")
    new_data = load_data_folder(new_samples, data_folder_pth)

    if new_data is not None:
        updated_data = update_data_store(data_store_pth, new_data)
        updated_data.to_parquet(data_store_pth)


def load_data_store(data_store_pth: Path) -> Set[str]:
    if data_store_pth.exists():
        data_store = pd.read_parquet(data_store_pth)
        return set(data_store.index.unique())
    return set()


def find_new_samples(workflow_samples: set, old_samples: set, data_folder_pth: Path) -> Set[str]:
    dir_content = set([file_name.stem for file_name in data_folder_pth.iterdir()])
    return workflow_samples.intersection(dir_content - old_samples)


def load_data_folder(samples: set, data_pth: Path) -> Optional[pd.DataFrame]:
    if len(samples) == 0:
        return None

    return pd.concat(
        [pd.read_parquet(data_pth / f"{sample}.parquet") for sample in samples],
        sort=False,
    )


def update_data_store(data_store_pth: Path, new_data: pd.DataFrame) -> pd.DataFrame:
    try:
        return pd.concat([pd.read_parquet(data_store_pth), new_data], sort=False)
    except OSError:
        return new_data


if __name__ == "__main__":
    main()
