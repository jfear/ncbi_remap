"""Looks at R1 and R2 to determine library layout"""
import os
import sys
from pathlib import Path
from collections import namedtuple

from ncbi_remap.fastq import check_fastq, fastq_stats, fastq_abi_solid
from ncbi_remap.snakemake import put_flag

DIRNAME = Path(snakemake.input.r1).parent
OUTPUT_HEADER = "\t".join(
    ["md5_R1", "libsize_R1", "avgLen_R1", "md5_R2", "libsize_R2", "avgLen_R2"]
)

ReadSummary = namedtuple("ReadSummary", "md5 size len abi")


class DownloadException(Exception):
    """Basic exception for problems downloading from SRA"""


class AbiException(Exception):
    """Basic exception when ABI file was downloaded from SRA"""


def main():
    try:
        verify_download(snakemake.input.download_log)
        r1 = summarize(snakemake.input.r1)
        r2 = summarize(snakemake.input.r2)
        create_flags(r1, r2)
        save_summary(r1, r2)
    except (DownloadException, AbiException):
        Path(snakemake.input.sra).unlink()
        Path(snakemake.input.r1).unlink()
        Path(snakemake.input.r2).unlink()


def verify_download(file_name):
    with open(file_name) as fh:
        text = fh.read()
        if "download successfully" not in text:
            Path(DIRNAME / "DOWNLOAD_BAD").touch()
            raise DownloadException


def summarize(file_name) -> ReadSummary:
    # TODO: handle GZ files.
    if check_fastq(file_name):
        return ReadSummary(None, None, None)

    lib_size, avg_len = fastq_stats(file_name)
    abi = fastq_abi_solid(file_name)
    return ReadSummary("Deprecated", lib_size, avg_len, abi)


def create_flags(r1: ReadSummary, r2: ReadSummary):
    if (r1.size > 1000) & (r2.size is None) & (r1.len > 10) & (r2.len is None):
        # R2 does not exists so SE
        put_flag(snakemake.output.layout, "SE")

    elif (r2.size is None) & (r2.len is None):
        # SE but R1 looks bad, flag as download bad
        Path(DIRNAME, "DOWNLOAD_BAD").touch()
        raise DownloadException

    elif (r1.size > 1000) & (r2.size > 1000) & (r1.len > 10) & (r2.len > 10):
        if r1.size == r2.size:
            # Both reads look good so PE
            put_flag(snakemake.output.layout, "PE")
        else:
            # There are an uneven number of reads between R1 and R2.
            # Instead of messing with this, just consider SE and use R1.
            put_flag(snakemake.output.layout, "keep_R1")

    elif (r1.size > 1000) & (r1.len > 10):
        # Only R1 looks ok, so SE
        put_flag(snakemake.output.layout, "keep_R1")

    elif (r2.size > 1000) & (r2.len > 10):
        # Only R2 looks ok, consider single-end
        put_flag(snakemake.output.layout, "keep_R2")

    else:
        # R1 and R2 look bad, flag as download bad
        Path(DIRNAME, "DOWNLOAD_BAD").touch()
        raise DownloadException

    if r1.abi | r2.abi:
        Path(DIRNAME, "ABI_SOLID").touch()
        raise AbiException


def save_summary(r1: ReadSummary, r2: ReadSummary):
    with open(snakemake.output.summary, "w") as file_out:
        file_out.write(OUTPUT_HEADER + "\n")
        file_out.write("\t".join([*r1[:-1], *r2[:-1]] + "\n"))


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                sra="../../output/fastq-wf/sra_cache/SRR031744.sra",
                r1="../../output/fastq-wf/fastqs/SRR031744_1.fastq.gz",
                r2="../../output/fastq-wf/fastqs/SRR031744_2.fastq.gz",
                download_log="../../output/fastq-wf/sra_download_logs/SRR031744.log",
            ),
            output=dict(
                layout="../../output/fastq-wf/fastq_info/SRR031744/LAYOUT",
                summary="../../output/fastq-wf/fastq_info/SRR031744/summary.tsv",
            ),
        )

    main()
