"""Looks at R1 and R2 to determine library layout"""
import os
import sys
from pathlib import Path
from collections import namedtuple

import numpy as np
import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.fastq import fastq_is_empty, fastq_read_stats, fastq_is_abi_solid

ReadSummary = namedtuple("ReadSummary", "library_size read_length")


class DownloadException(Exception):
    """Basic exception for problems downloading from SRA"""


class AbiException(Exception):
    """Basic exception when ABI file was downloaded from SRA"""


def main():
    try:
        verify_sra_download(snakemake.input.download_log)
        r1 = summarize(snakemake.input.r1)
        r2 = summarize(snakemake.input.r2)
        create_flags(r1, r2)
        save_summary(r1, r2)
    except AbiException:
        Path(snakemake.params.abi_solid, snakemake.wildcards.srr).touch()
        remove_sample_due_to_problems()
    except DownloadException:
        Path(snakemake.params.download_bad, snakemake.wildcards.srr).touch()
        remove_sample_due_to_problems()


def verify_sra_download(download_log: str):
    with open(download_log) as fh:
        text = fh.read()
        if "downloaded successfully" not in text:
            raise DownloadException


def summarize(fastq: str) -> ReadSummary:
    """Summarize FASTQ file."""
    with open(fastq, "r") as fastq_handle:
        if fastq_is_empty(fastq_handle):
            return ReadSummary(0, 0.0)

        fastq_handle.seek(0, 0)  # Go back to the begining of the file
        if fastq_is_abi_solid(fastq_handle):
            raise AbiException

        fastq_handle.seek(0, 0)  # Go back to the begining of the file
        library_size, avg_read_length = fastq_read_stats(fastq_handle)
        return ReadSummary(library_size, avg_read_length)


def create_flags(r1: ReadSummary, r2: ReadSummary):
    if (r1.library_size > 1000) & (r1.read_length > 10) & (r2.library_size == 0):
        # SE because R1 looks good and R2 is empty
        layout = "SE"
    elif (
        (r1.library_size > 1000)
        & (r1.read_length > 10)
        & (r2.library_size > 1000)
        & (r2.read_length > 10)
    ):
        if r1.library_size == r2.library_size:
            # PE because R1 and R2 look good and have the same number of reads
            layout = "PE"
        else:
            # SE/keep_R1 because R1 and R2 look good, but uneven number of reads
            layout = "keep_R1"
    elif (r1.library_size > 1000) & (r1.read_length > 10):
        # SE/keep_R1 because only R1 looks good
        layout = "keep_R1"
    elif (r2.library_size > 1000) & (r2.read_length > 10):
        # SE/keep_R2 because only R2 looks good
        layout = "keep_R2"
    else:
        # We have a problem because R1 and R2 look bad
        raise DownloadException

    idx = pd.Index([snakemake.wildcards.srr], name="srr")
    df = pd.DataFrame([[layout]], index=[idx], columns=["layout"])
    df.to_parquet(snakemake.output.layout)


def save_summary(r1: ReadSummary, r2: ReadSummary):
    idx = pd.Index([snakemake.wildcards.srr], name="srr")
    df = pd.DataFrame(
        [[r1.library_size, r1.read_length, r2.library_size, r2.read_length,]],
        index=idx,
        columns=["libsize_R1", "avgLen_R1", "libsize_R2", "avgLen_R2"],
    )
    df.to_parquet(snakemake.output.summary)


def remove_sample_due_to_problems():
    try:
        Path(snakemake.input.sra).unlink()
        Path(snakemake.input.r1).unlink()
        Path(snakemake.input.r2).unlink()
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                sra="../../output/fastq-wf/sra_cache/SRR031744.sra",
                r1="../../output/fastq-wf/fastqs/SRR031744_1.fastq",
                r2="../../output/fastq-wf/fastqs/SRR031744_2.fastq",
                download_log="../../output/fastq-wf/sra_download_logs/SRR031744.log",
            ),
            params=dict(
                download_bad="../../output/fastq-wf/download_bad",
                abi_solid="../../output/fastq-wf/abi_solid",
            ),
            output=dict(
                layout="../output/fastq-wf/layout/SRR031744.parquet",
                summary="../output/fastq-wf/libsize/SRR031744.parquet"
            ),
        )

    Path(snakemake.params.download_bad).mkdir(exist_ok=True)
    Path(snakemake.params.abi_solid).mkdir(exist_ok=True)
    main()
