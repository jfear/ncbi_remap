"""Migrates to fastq-wf/fastq_info/{srr}/[LAYOUT, summary.tsv]"""
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_srr
from ncbi_remap.snakemake import get_flag

FASTQ_INFO_PATH = Path("../output/fastq-wf/fastq_info")


def main():
    convert_layout()
    convert_summary()
    rmdir()

def convert_layout():
    outdir = Path("../output/fastq-wf/layout").mkdir(exist_ok=True)
    for pth in FASTQ_INFO_PATH.glob("**/LAYOUT"):
        srr = parse_srr(pth.as_posix())
        print(f"Moving Layout:\t{srr}")
        idx = pd.Index([srr], name="srr")
        layout = get_flag(pth)
        df = pd.DataFrame([[layout]], index=[idx], columns=["layout"])
        df.to_parquet(f"../output/fastq-wf/layout/{srr}.parquet")
        pth.unlink()


def convert_summary():
    outdir = Path("../output/fastq-wf/libsize").mkdir(exist_ok=True)
    dtypes = {
        "libsize_R1": np.int64,
        "avgLen_R1": np.float64,
        "libsize_R2": np.int64,
        "avgLen_R2": np.float64,
    }

    for pth in FASTQ_INFO_PATH.glob("**/summary.tsv"):
        srr = parse_srr(pth.as_posix())
        print(f"Moving Libsize:\t{srr}")
        idx = pd.Index([srr,], name="srr")
        df = pd.read_table(pth)[dtypes.keys()].fillna(0).astype(dtypes)
        df.index = idx
        df.to_parquet(f"../output/fastq-wf/libsize/{srr}.parquet")
        pth.unlink()

def rmdir():
    for pth in FASTQ_INFO_PATH.iterdir():
        if not pth.is_dir():
            continue
        pth.rmdir()


if __name__ == "__main__":
    main()
