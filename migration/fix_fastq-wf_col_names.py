"""Change 'SRR' to 'srr'.

I accidentally changed the colname srr to SRR. I need to go through a fix these.
"""
from pathlib import Path

import pandas as pd


def main():
    for pth in Path("../output/fastq-wf/layout").iterdir():
        fix_column(pth)

    for pth in Path("../output/fastq-wf/libsize").iterdir():
        fix_column(pth)


def fix_column(file_name: Path):
    df = pd.read_parquet(file_name)
    if df.index.name == "SRR":
        df.index.name = "srr"
        df.to_parquet(file_name)


if __name__ == "__main__":
    main()
