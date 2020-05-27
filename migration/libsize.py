"""Convert libsize tables

From:
srr:index
libsize_R1:int
avgLen_R1:float
libsize_R2:int
avgLen_R2:float

To:
srr:index
libsize:int
avgLen_R1:float
avgLen_R2:float
"""
from pathlib import Path

import pandas as pd

PROJECT_DIR = Path(__file__).absolute().parents[1]
LIBSIZE_DIR = PROJECT_DIR / "output/fastq-wf/libsize"


def main():
    for file_name in LIBSIZE_DIR.iterdir():
        df = pd.read_parquet(file_name)
        if "libsize_R1" in df.columns:
            df["libsize"] = df[["libsize_R1", "libsize_R2"]].max(axis=1)
            df[["libsize", "avgLen_R1", "avgLen_R2"]].astype(
                {"libsize": int, "avgLen_R1": float, "avgLen_R2": float,}
            ).to_parquet(file_name)


if __name__ == "__main__":
    main()
