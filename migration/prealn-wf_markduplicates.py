import sys
from pathlib import Path
from collections import namedtuple

import numpy as np
import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_picard_markduplicate_metrics


PREALN_PATH = Path("../output/prealn-wf/samples")
DTYPES = {
    "UNPAIRED_READS_EXAMINED": np.int64,
    "READ_PAIRS_EXAMINED": np.int64,
    "UNPAIRED_READ_DUPLICATES": np.int64,
    "READ_PAIR_DUPLICATES": np.int64,
    "PERCENT_DUPLICATION": np.float64,
    "ESTIMATED_LIBRARY_SIZE": np.int64,
}

Files = namedtuple("Files", "srr idx bam metrics log output")


def main():
    for srr_pth in PREALN_PATH.glob("**/SRR*"):
        if not srr_pth.is_dir():
            continue
        srr = srr_pth.name
        files = Files(
            srr,
            pd.Index([srr], name="srr"),
            (srr_pth / f"{srr}.hisat2.bam.picard.markduplicates.bam"),
            (srr_pth / f"{srr}.hisat2.bam.picard.markduplicates.metrics"),
            (srr_pth / f"{srr}.hisat2.bam.picard.markduplicates.metrics.log"),
            f"../output/prealn-wf/markduplicates/{srr}.parquet",
        )

        if files.metrics.exists():
            parse_table(files)
            remove_file(files.bam)
            remove_file(files.metrics)
            remove_file(files.log)


def parse_table(files):
    df = parse_picard_markduplicate_metrics(files.metrics)[DTYPES.keys()].fillna(0).astype(DTYPES)
    df.PERCENT_DUPLICATION = df.PERCENT_DUPLICATION * 100

    df.columns = [col.lower() for col in df.columns]
    df.index = files.idx

    df.to_parquet(files.output)


def remove_file(file_name):
    if Path(file_name).exists():
        Path(file_name).unlink()


if __name__ == "__main__":
    main()
