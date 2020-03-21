from pathlib import Path
from collections import namedtuple

import numpy as np
import pandas as pd


PREALN_PATH = Path("../output/prealn-wf/samples")

Files = namedtuple("Files", "srr idx counts jcounts summary log output")


class LogException(Exception):
    """File was not complete"""


def main():
    for srr_pth in PREALN_PATH.glob("**/SRR*"):
        if not srr_pth.is_dir():
            continue
        srr = srr_pth.name
        files = Files(
            srr,
            pd.Index([srr], name="srr"),
            (srr_pth / f"{srr}.hisat2.bam.feature_counts.counts"),
            (srr_pth / f"{srr}.hisat2.bam.feature_counts.counts.jcounts"),
            (srr_pth / f"{srr}.hisat2.bam.feature_counts.counts.summary"),
            (srr_pth / f"{srr}.hisat2.bam.feature_counts.counts.log"),
            f"../output/prealn-wf/count_summary/{srr}.parquet",
        )

        if not files.log.exists() & files.counts.exists() & files.jcounts.exists():
            continue

        try:
            check_log(files.log)
            parse_table(files)
        except (FileNotFoundError, LogException):
            print(files.log)
            pass

        remove_file(files.counts)
        remove_file(files.jcounts)
        remove_file(files.summary)
        remove_file(files.log)


def check_log(file_name):
    # Check log for completeness
    with open(file_name, "r") as fh:
        log = fh.read()
        if not ("Read assignment finished" in log or "Summary of counting results" in log):
            raise LogException


def get_counts(file_name):
    col_name = pd.read_table(file_name, comment="#", nrows=1).columns[-1]
    return pd.read_table(file_name, comment="#", usecols=[col_name], dtype=np.int64).values


def parse_table(files):
    gene_counts = get_counts(files.counts)
    genic_reads = gene_counts.sum()
    percent_genes_on = (gene_counts > 0).mean() * 100

    junction_counts = get_counts(files.jcounts)
    junction_reads = junction_counts.sum()
    number_junctions_on = junction_counts.shape[0]

    df = pd.DataFrame(
        [[genic_reads, percent_genes_on, junction_reads, number_junctions_on]],
        columns=[
            "number_genic_reads",
            "percent_genes_on",
            "number_junction_reads",
            "number_junctions_on",
        ],
        index=files.idx,
    )

    df.to_parquet(files.output)


def remove_file(file_name):
    if Path(file_name).exists():
        Path(file_name).unlink()


if __name__ == "__main__":
    main()
