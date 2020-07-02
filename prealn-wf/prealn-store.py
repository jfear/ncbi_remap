#!/usr/bin/env python
"""Aggregate pre-alignment workflow data"""
from pathlib import Path
from typing import Optional, Set, Tuple

import pandas as pd

PREALN_OUTPUT = Path(__file__).absolute().parents[1] / "output/prealn-wf"
OUTPUTS = [
    "aln_stats",
    "atropos",
    "count_summary",
    "fastq_screen",
    "genebody_coverage",
    "hisat2",
    "markduplicates",
    "rnaseqmetrics",
    "strand",
]


def main():
    workflow_samples = load_complete_samples()
    for output in OUTPUTS:
        print(f"Aggregating: {output:>20}", end="\t")
        aggregate_data_store(
            workflow_samples, PREALN_OUTPUT / output, PREALN_OUTPUT / f"{output}.parquet"
        )


def load_complete_samples() -> Set[str]:
    with (PREALN_OUTPUT / "done.txt").open() as fh:
        return {srr.strip() for srr in fh}


def aggregate_data_store(workflow_samples: set, data_folder_pth: Path, data_store_pth: Path):
    data_store, old_samples = load_data_store(data_store_pth)
    new_samples = find_new_samples(workflow_samples, old_samples, data_folder_pth)

    print(f"({len(new_samples):,})")
    new_data = load_data_folder(new_samples, data_folder_pth)

    if new_data is not None:
        updated_data = update_data_store(data_store, new_data)
        updated_data.to_parquet(data_store_pth)


def load_data_store(data_store_pth: Path) -> Tuple[Optional[pd.DataFrame], Set[str]]:
    if data_store_pth.exists():
        data_store = pd.read_parquet(data_store_pth)
        old_samples = set(data_store.index.unique())
    else:
        data_store, old_samples = None, set()
    return data_store, old_samples


def find_new_samples(workflow_samples: set, old_samples: set, data_folder_pth: Path) -> Set[str]:
    dir_content = {file_name.stem for file_name in data_folder_pth.iterdir()}
    return workflow_samples.intersection(dir_content - old_samples)


def load_data_folder(samples: set, data_pth: Path) -> Optional[pd.DataFrame]:
    if not samples:
        return None

    return pd.concat(
        [pd.read_parquet(data_pth / f"{sample}.parquet") for sample in samples], sort=False
    )


def update_data_store(data_store: Optional[pd.DataFrame], new_data: pd.DataFrame) -> pd.DataFrame:
    if data_store:
        return pd.concat([data_store, new_data], sort=False)
    return new_data


if __name__ == "__main__":
    main()
