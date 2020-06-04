"""Extract SRA file into Gziped FASTQs"""
import os
import shutil
import sys
from gzip import open as gopen
from pathlib import Path
from subprocess import CalledProcessError
from typing import Optional, Tuple

import pandas as pd
from snakemake import logger
from snakemake.shell import shell

sys.path.insert(0, "../src")
from ncbi_remap.fastq import Fastq, MixedUpReadsException, UnequalNumberReadsException
from ncbi_remap.io import remove_file
from ncbi_remap.snakemake import StepLogger

LOG = StepLogger(str(snakemake.log))
SRR = snakemake.wildcards.srr
TMPDIR = Path(os.getenv("TMPDIR", "/tmp"), SRR)
TMPDIR.mkdir(exist_ok=True)


def main():
    sra_file = get_sra(SRR, snakemake.params.sra)
    r1, r2 = fastq_dump(SRR, sra_file)
    r1_gz, r2_gz, fq = check_and_compress_fastq(r1, r2)
    save_output(
        r1_gz,
        r2_gz,
        fq,
        snakemake.output.r1,
        snakemake.output.r2,
        snakemake.output.layout,
        snakemake.output.summary,
    )


class DownloadException(Exception):
    """Basic exception for problems downloading from SRA"""


class AbiException(Exception):
    """Basic exception when ABI file was downloaded from SRA"""


def get_sra(SRR: str, sra_file: str) -> Path:
    if Path(sra_file).exists():
        return stage_data(sra_file)
    return sra_prefetch(SRR, sra_file)


def stage_data(file_name: str) -> Path:
    scratch = TMPDIR / Path(file_name).name
    shutil.copy2(file_name, scratch)
    return scratch


def sra_prefetch(SRR: str, sra_file: str) -> Path:
    """Runs SRA Prefetch"""
    scratch = TMPDIR / f"{SRR}.sra"
    log = TMPDIR / "sra.log"
    shell(f"prefetch --output-file {scratch} " f"--max-size 40G {SRR} " f"> {log} 2>&1")
    LOG.append("prefetch", log)
    verify_sra_download(log)
    log.unlink()

    # copy from scratch
    Path(sra_file).parent.mkdir(exist_ok=True)
    shutil.copy2(scratch, sra_file)
    return scratch


def verify_sra_download(download_log: Path) -> None:
    with download_log.open() as fh:
        text = fh.read()
        if "downloaded successfully" not in text:
            raise DownloadException("Prefetch Failed")


def fastq_dump(SRR: str, sra_file: Path) -> Tuple[str, Optional[str]]:
    log = TMPDIR / f"fastq_dump.log"
    try:
        shell(
            f"fastq-dump {sra_file} "
            "--split-3 "
            "--skip-technical "
            "--minReadLen 20 "
            f"-O {TMPDIR} "
            f"> {log} 2>&1"
        )
    finally:
        LOG.append("fastq-dump", log)
        log.unlink()
        sra_file.unlink()

    # Determine naming scheme
    se_R1 = TMPDIR / f"{SRR}.fastq"
    if se_R1.exists():
        return se_R1.as_posix(), None

    pe_R1 = TMPDIR / f"{SRR}_1.fastq"
    pe_R2 = TMPDIR / f"{SRR}_2.fastq"
    if pe_R1.exists() and pe_R2.exists():
        return pe_R1.as_posix(), pe_R2.as_posix()

    raise DownloadException("No FASTQs extracted")


def check_and_compress_fastq(r1: str, r2: Optional[str]) -> Tuple[Path, Path, Fastq]:
    fq = Fastq(r1, r2)
    r1_gz = TMPDIR / f"{SRR}_1.fastq.gz"
    r2_gz = TMPDIR / f"{SRR}_2.fastq.gz"

    if r2 is None:
        logger.info("Processing FASTQ as Single-End")
        run_as_se(fq, r1_gz, r2_gz)
    else:
        logger.info("Processing FASTQ as Pair-End")
        run_as_pe(fq, r1_gz, r2_gz)
    logger.info(fq)
    LOG.append("FASTQ Check", text=str(fq))
    return r1_gz, r2_gz, fq


def run_as_se(fq: Fastq, R1_out: Path, R2_out: Path) -> None:
    R2_out.touch()  # always output R2
    with gopen(R1_out, "wb") as file_out1:
        for read in fq.process():
            file_out1.write(read)

    if "abi_solid" in fq.flags:
        raise AbiException

    if "download_bad" in fq.flags:
        raise DownloadException("Empty FASTQ")

    if fq.libsize < 100_000:
        raise DownloadException("<100,000 reads")


def run_as_pe(fq: Fastq, R1_out: Path, R2_out: Path) -> None:
    try:
        with gopen(R1_out, "wb") as file_out1, gopen(R2_out, "wb") as file_out2:
            for read1, read2 in fq.process():
                file_out1.write(read1)
                file_out2.write(read2)

        if "abi_solid" in fq.flags:
            raise AbiException

        if "download_bad" in fq.flags:
            raise DownloadException("Empty FASTQ")

        if fq.libsize < 100_000:
            raise DownloadException("<100,000 reads")

    except (UnequalNumberReadsException, MixedUpReadsException):
        remove_file(R1_out)
        remove_file(R2_out)
        run_as_se(fq, R1_out, R2_out)


def save_output(
    r1_gz: Path,
    r2_gz: Path,
    fq: Fastq,
    r1_out: Path,
    r2_out: Path,
    layout_out: str,
    summary_out: str,
) -> None:
    shutil.copy2(r1_gz, r1_out)
    shutil.copy2(r2_gz, r2_out)
    save_layout(fq, layout_out)
    save_summary(fq, summary_out)


def save_layout(fq: Fastq, layout_file) -> None:
    layout = fq.flags.intersection(set(["SE", "PE", "keep_R1", "keep_R2"])).pop()
    idx = pd.Index([SRR], name="SRR")
    df = pd.DataFrame([[layout]], index=[idx], columns=["layout"])
    df.to_parquet(layout_file)


def save_summary(fq: Fastq, summary_file) -> None:
    if isinstance(fq.avgReadLen, list):
        r1, r2 = fq.avgReadLen
    else:
        r1, r2 = fq.avgReadLen, 0.0

    idx = pd.Index([SRR], name="SRR")
    df = pd.DataFrame(
        [[fq.libsize, r1, r2]], index=idx, columns=["libsize", "avgLen_R1", "avgLen_R2"],
    )
    df.to_parquet(summary_file)


def remove_outputs(outputs) -> None:
    for output in outputs:
        remove_file(output)


if __name__ == "__main__":
    try:
        main()
    except AbiException:
        logger.warning(f"Flagging {SRR} as ABI Solid")
        LOG.append("Exception", text="Abi Solid")
        pth = Path(snakemake.params.abi_solid)
        pth.mkdir(exist_ok=True)
        (pth / SRR).touch()
        remove_outputs(snakemake.output)
    except DownloadException as error:
        logger.warning(f"Flagging {SRR} as Download Bad")
        logger.warning(str(error))
        LOG.append("Exception", text="Download Bad")
        pth = Path(snakemake.params.download_bad)
        pth.mkdir(exist_ok=True)
        (pth / SRR).touch()
        remove_outputs(snakemake.output)
    finally:
        shutil.rmtree(TMPDIR)
