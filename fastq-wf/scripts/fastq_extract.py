"""Extract SRA file into Gziped FASTQs"""
import os
import sys
import shutil
from gzip import open as gopen
from pathlib import Path
from subprocess import CalledProcessError

import pandas as pd
from snakemake.shell import shell
from snakemake import logger

sys.path.insert(0, "../src")
from ncbi_remap.fastq import Fastq, UnequalNumberReadsException, MixedUpReadsException

TMPDIR = os.getenv("TMPDIR", "/tmp")


class DownloadException(Exception):
    """Basic exception for problems downloading from SRA"""


class AbiException(Exception):
    """Basic exception when ABI file was downloaded from SRA"""


def main():
    try:
        verify_sra_download(snakemake.input.log)
        staged_data = stage_data(snakemake.input.sra)
        extract_fastq(staged_data)

        if Path(TMPDIR, f"{snakemake.wildcards.srr}.sra.fastq").exists():
            logger.info("Processing FASTQ as Single-End")
            R1 = Path(TMPDIR, f"{snakemake.wildcards.srr}.sra.fastq").absolute().as_posix()
            R2 = None
            fq = Fastq(R1)
            check_se(fq, snakemake.output.r1, snakemake.output.r2)
        else:
            logger.info("Processing FASTQ as Pair-End")
            R1 = Path(TMPDIR, f"{snakemake.wildcards.srr}.sra_1.fastq").absolute().as_posix()
            R2 = Path(TMPDIR, f"{snakemake.wildcards.srr}.sra_2.fastq").absolute().as_posix()
            fq = Fastq(R1, R2)
            check_pe(fq, snakemake.output.r1, snakemake.output.r2)

        logger.info(fq)
        save_layout(fq, snakemake.output.layout)
        save_summary(fq, snakemake.output.summary)

    except AbiException:
        logger.warning(f"Flagging {snakemake.wildcards.srr} as ABI Solid")
        Path(snakemake.params.abi_solid, snakemake.wildcards.srr).touch()
        remove_file(snakemake.output.r1)
        remove_file(snakemake.output.r2)
        remove_file(snakemake.input.sra)
    except DownloadException:
        logger.warning(f"Flagging {snakemake.wildcards.srr} as Download Bad")
        Path(snakemake.params.download_bad, snakemake.wildcards.srr).touch()
        remove_file(snakemake.output.r1)
        remove_file(snakemake.output.r2)
        remove_file(snakemake.input.sra)
    finally:
        logger.info(f"Cleaning up {TMPDIR}")
        remove_file(R1)
        remove_file(R2)
        remove_file(staged_data)


def verify_sra_download(download_log: str):
    with open(download_log) as fh:
        text = fh.read()
        if "downloaded successfully" not in text:
            logger.error(f"SRA file for {snakemake.wildcards.srr} was not downloaded successfully")
            raise DownloadException


def stage_data(file_name) -> str:
    """Copy data over to local scratch space."""
    data_pth = Path(file_name)
    tmp_pth = Path(TMPDIR, data_pth.name)

    logger.info(f"Staging SRA data to: {tmp_pth.as_posix()}")
    shutil.copy2(data_pth, tmp_pth)

    return tmp_pth.as_posix()


def extract_fastq(sra_file):
    try:
        shell(
            f"fasterq-dump {sra_file} "
            "--split-3 "
            "--skip-technical "
            "--min-read-len 20 "
            f"-O {TMPDIR} "
            f"-t {TMPDIR} "
            f"-e {snakemake.threads} "
        )
    except CalledProcessError:
        shell(
            f"fastq-dump {sra_file} "
            "--split-3 "
            "--skip-technical "
            "--minReadLen 20 "
            f"-O {TMPDIR} "
        )


def check_se(fq: Fastq, R1_out: str, R2_out: str):
    Path(R2_out).touch()  # always output R2
    with gopen(R1_out, "wb") as file_out1:
        for read in fq.process():
            file_out1.write(read)

    if "download_bad" in fq.flags:
        logger.error(f"{snakemake.wildcards.srr} had no data.")
        raise DownloadException

    if fq.libsize < 1000:
        logger.warning(f"{snakemake.wildcards.srr} had less than 1,000 reads.")
        raise DownloadException

    if "abi_solid" in fq.flags:
        raise AbiException


def check_pe(fq: Fastq, R1_out: str, R2_out: str):
    try:
        with gopen(R1_out, "wb") as file_out1, gopen(R2_out, "wb") as file_out2:
            for read1, read2 in fq.process():
                file_out1.write(read1)
                file_out2.write(read2)

        if "download_bad" in fq.flags:
            logger.error(f"{snakemake.wildcards.srr} had no data.")
            raise DownloadException

        if fq.libsize < 1000:
            logger.warning(f"{snakemake.wildcards.srr} had less than 1,000 reads.")
            raise DownloadException

        if "abi_solid" in fq.flags:
            raise AbiException

    except (UnequalNumberReadsException, MixedUpReadsException):
        logger.warning(f"{snakemake.wildcards.srr} switching from PE to SE.")
        remove_file(R1_out)
        remove_file(R2_out)
        check_se(fq, R1_out, R2_out)


def save_layout(fq: Fastq, layout_file):
    layout = fq.flags.intersection(set(["SE", "PE", "keep_R1", "keep_R2"])).pop()
    idx = pd.Index([snakemake.wildcards.srr], name="srr")
    df = pd.DataFrame([[layout]], index=[idx], columns=["layout"])
    df.to_parquet(layout_file)


def save_summary(fq: Fastq, summary_file):
    if isinstance(fq.avgReadLen, list):
        r1, r2 = fq.avgReadLen
    else:
        r1, r2 = fq.avgReadLen, 0.0

    idx = pd.Index([snakemake.wildcards.srr], name="srr")
    df = pd.DataFrame(
        [[fq.libsize, r1, r2]],
        index=idx,
        columns=["libsize", "avgLen_R1", "avgLen_R2"],
    )
    df.to_parquet(summary_file)


def remove_file(file_name: str):
    if file_name is None:
        return 

    if Path(file_name).exists():
        Path(file_name).unlink()


if __name__ == "__main__":
    main()
