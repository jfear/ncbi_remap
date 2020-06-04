import os
import shutil
import sys
import stat
from time import sleep
from pathlib import Path
from subprocess import SubprocessError
from typing import Optional, Tuple

from snakemake.logging import logger
from snakemake.shell import shell
import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.io import remove_folder, remove_file
from ncbi_remap.snakemake import StepLogger
from ncbi_remap.parser import parse_picardCollect_summary, parse_picardCollect_hist


LOG = StepLogger(str(snakemake.log))
SRR = snakemake.wildcards.srr
THREADS = snakemake.threads
MEM = int(snakemake.resources.get("mem_gb", 4))

TMPDIR = Path(os.getenv("TMPDIR", "/tmp"), f"{SRR}/collectrnaseqmetrics")
TMPDIR.mkdir(parents=True, exist_ok=True)

TMPREF = Path(os.getenv("TMPDIR", "/tmp"), f"references")
TMPREF.mkdir(exist_ok=True)


def main():
    bam, ref_flat = stage_data(snakemake.input.bam, snakemake.input.refflat)
    unstranded, first, second = picard(bam, ref_flat)
    remove_file(bam)
    summarize(
        unstranded,
        first,
        second,
        snakemake.output.strand,
        snakemake.output.table,
        snakemake.output.genebody_coverage,
    )


def stage_data(bam: str, ref_flat: str) -> Tuple[Path, Path]:
    bam_local = TMPDIR / f"{SRR}.bam"
    shutil.copy2(bam, bam_local)

    ref_flat_local = TMPREF / "dmel.refflat"
    # Stage references if not present
    if ref_flat_local.exists():
        # Already on scratch wait to make sure copied
        sleep(5)
    else:
        shutil.copy(ref_flat, ref_flat_local)
        os.chmod(
            ref_flat_local,
            stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH | stat.S_IWOTH,
        )

    return bam_local, ref_flat_local


def picard(bam: Path, ref_flat: Path) -> Tuple[Path, Path, Path]:
    cmd = (
        f"picard -Xmx{MEM}g CollectRnaSeqMetrics REF_FLAT={ref_flat} "
        f"INPUT={bam} OUTPUT={{out_file}} STRAND={{strand}} > {{log}} 2>&1"
    )
    LOG.append("Picard General Command", text=cmd)

    # Unstranded
    un_out = TMPDIR / "unstranded.txt"
    log = TMPDIR / "unstranded.log"
    shell(cmd.format(out_file=un_out, strand="NONE", log=log))
    LOG.append("Picard Unstranded", log)
    _check_log(log)
    remove_file(log)

    # First stranded
    first_out = TMPDIR / "first_stranded.txt"
    log = TMPDIR / "first_stranded.log"
    shell(cmd.format(out_file=first_out, strand="FIRST_READ_TRANSCRIPTION_STRAND", log=log))
    LOG.append("Picard Unstranded", log)
    _check_log(log)
    remove_file(log)

    # Second stranded
    second_out = TMPDIR / "second_stranded.txt"
    log = TMPDIR / "second_stranded.log"
    shell(cmd.format(out_file=second_out, strand="SECOND_READ_TRANSCRIPTION_STRAND", log=log))
    LOG.append("Picard Unstranded", log)
    _check_log(log)
    remove_file(log)

    return un_out, first_out, second_out


def _check_log(log: Path) -> None:
    with log.open() as fh:
        if "CollectRnaSeqMetrics done" not in fh.read():
            raise PicardException(f"{log.stem} not complete")


def summarize(
    unstranded: Path, first: Path, second: Path, strand: str, table: str, genebody: str
) -> None:
    # Parse the flags
    if _parse_stranded(first):
        strand_ = "same_strand"
    elif _parse_stranded(second):
        strand_ = "opposite_strand"
    else:
        strand_ = "unstranded"

    # Write strand flag
    LOG.append("Strand", text=strand_)
    idx = pd.Index([SRR], name="srr")
    df = pd.DataFrame([[strand_]], index=idx, columns=["strand"])
    df.to_parquet(strand)

    # Parse main table
    df = _parse_table(unstranded)
    df.index = idx
    df.to_parquet(table)

    # Parse genome coverage histogram
    df = parse_picardCollect_hist(unstranded)
    df.index = idx
    df.to_parquet(genebody)


def _parse_stranded(file_name: Path) -> bool:
    return (parse_picardCollect_summary(file_name).PCT_CORRECT_STRAND_READS >= 0.75)[0]


def _parse_table(file_name: Path) -> pd.DataFrame:
    df = parse_picardCollect_summary(file_name)[
        [
            "PCT_CODING_BASES",
            "PCT_UTR_BASES",
            "PCT_INTRONIC_BASES",
            "PCT_INTERGENIC_BASES",
            "PCT_MRNA_BASES",
            "MEDIAN_CV_COVERAGE",
            "MEDIAN_5PRIME_BIAS",
            "MEDIAN_3PRIME_BIAS",
            "MEDIAN_5PRIME_TO_3PRIME_BIAS",
        ]
    ].fillna(0.0)
    df.PCT_CODING_BASES = df.PCT_CODING_BASES * 100
    df.PCT_UTR_BASES = df.PCT_UTR_BASES * 100
    df.PCT_INTRONIC_BASES = df.PCT_INTRONIC_BASES * 100
    df.PCT_INTERGENIC_BASES = df.PCT_INTERGENIC_BASES * 100
    df.PCT_MRNA_BASES = df.PCT_MRNA_BASES * 100

    df.columns = [x.lower() for x in df.columns]
    df.columns = [x.replace("pct_", "percent_") for x in df.columns]
    return df


class PicardException(Exception):
    """Picard Processing Exception"""


if __name__ == "__main__":
    try:
        main()
    except PicardException as error:
        logger.warning(f"{SRR}: {error}")
        LOG.append("Exception", text=str(error))
    finally:
        remove_folder(TMPDIR)
