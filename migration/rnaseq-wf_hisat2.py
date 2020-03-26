import sys
from collections import namedtuple
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_hisat2

CLEAN_UP = False
RNASEQ_PATH = Path("../output/rnaseq-wf/samples")

Files = namedtuple("Files", "srr idx log bam bai output")


class Hisat2Exception(Exception):
    """Basic exception when there are problems running Atropos"""


def main():
    for srr_pth in RNASEQ_PATH.glob("**/SRR*"):

        if not srr_pth.is_dir():
            continue

        srr = srr_pth.name
        files = Files(
            srr,
            pd.Index([srr], name="srr"),
            (srr_pth / f"{srr}.fq.bam.log"),
            (srr_pth / f"{srr}.fq.bam"),
            (srr_pth / f"{srr}.fq.bam.bai"),
            f"../output/rnaseq-wf/hisat2/{srr}.parquet",
        )
        if files.log.exists():
            convert_log(files)
            if CLEAN_UP:
                Path(files.log).unlink()

    remove_alignment_bad()


def convert_log(files):
    try:
        try:
            df = pd.read_table(files.log, index_col=[0, 1]).reset_index(drop=True).fillna(0)
        except IndexError:
            df = parse_hisat2(files.log).fillna(0)
        except pd.errors.EmptyDataError:
            raise Hisat2Exception

        if (df == 0).all().all():
            raise Hisat2Exception

        df.index = pd.Index([files.srr], name="srr")
        df.to_parquet(files.output)

        uniquely_aligned = (
            df.iloc[0, :]["num_concordant_reads_uniquely_aligned"]
            + df.iloc[0, :]["num_uniquely_aligned"]
        )

        if (df.iloc[0, :]["per_alignment"] < 1) | (uniquely_aligned < 1000):
            raise Hisat2Exception

    except Hisat2Exception:
        alignment_bad_path = Path("../output/rnaseq-wf/alignment_bad")
        alignment_bad_path.mkdir(exist_ok=True)
        Path(alignment_bad_path, files.srr).touch()

        if Path(files.bam).exists() and CLEAN_UP:
            Path(files.bam).unlink()

        if Path(files.bai).exists() and CLEAN_UP:
            Path(files.bai).unlink()


def remove_alignment_bad():
    srrs = set([pth.name for pth in Path("../output/rnaseq-wf/alignment_bad").iterdir()])
    for pth in Path("../output/rnaseq-wf/samples/").glob("**/SRR*"):
        if pth.name in srrs:
            for file_pth in pth.iterdir():
                if not "fastq_screen" in file_pth.name and CLEAN_UP:
                    file_pth.unlink()


if __name__ == "__main__":
    main()
