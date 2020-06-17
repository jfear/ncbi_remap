import os
import shutil
import sys
from time import sleep
from pathlib import Path
from typing import Optional, Tuple

from snakemake.logging import logger
from snakemake.shell import shell
import pandas as pd
import numpy as np

sys.path.insert(0, "../src")
from ncbi_remap.io import remove_folder, remove_file
from ncbi_remap.snakemake import StepLogger
from ncbi_remap.parser import parse_featureCounts_counts, parse_featureCounts_jcounts


LOG = StepLogger(str(snakemake.log))
SRX = snakemake.wildcards.srx
THREADS = snakemake.threads

TMPDIR = Path(os.getenv("TMPDIR", "/tmp"), f"{SRX}/{snakemake.rule}")
TMPDIR.mkdir(parents=True, exist_ok=True)

TMPREF = Path(os.getenv("TMPDIR", "/tmp"), f"references")
TMPREF.mkdir(exist_ok=True)


def main():
    bam, gtf = stage_data(snakemake.input.bam, snakemake.input.gtf)
    counts, jcounts = featurecount(bam, gtf, snakemake.input.layout, snakemake.input.strand)

    process_counts(counts, snakemake.output.counts)
    process_jcounts(jcounts, snakemake.output.get("jcounts"))


def stage_data(bam: str, gtf: str) -> Tuple[Path, Path]:
    # Copy Bam
    bam_local = TMPDIR / f"{SRX}.bam"
    shutil.copy2(bam, bam_local)

    # Stage GTF if not present
    gtf_local = TMPREF / Path(gtf).name
    if gtf_local.exists():
        # Already on scratch wait to make sure copied
        sleep(5)
    else:
        shutil.copy(gtf, gtf_local)

    return bam_local, gtf_local


def featurecount(bam: Path, gtf: Path, layout: str, strand: str) -> Tuple[Path, Path]:
    # Look up Layout
    layout_ = pd.read_parquet(layout).layout[0]
    if layout_ == "PE":
        params = snakemake.params.extra_pe
    else:
        params = snakemake.params.extra_se

    # Look up strand
    strand_ = pd.read_parquet(strand).strand[0]
    if strand_ == "same_strand":
        params += "-s 1"
    elif strand_ == "opposite_strand":
        params += "-s 2"
    else:
        params += "-s 0"

    counts = TMPDIR / f"{SRX}.counts"
    jcounts = TMPDIR / f"{SRX}.counts.jcounts"
    summary = TMPDIR / f"{SRX}.counts.summary"
    log = TMPDIR / f"featurecounts.log"

    try:
        shell(f"featureCounts -T {THREADS} {params} -a {gtf} -o {counts} {bam} > {log} 2>&1")
    finally:
        LOG.append("Feature Counts", log)

    _check_log(log)
    remove_file(log)
    remove_file(summary)

    return counts, jcounts


def _check_log(log: Path):
    with log.open() as fh:
        log_text = fh.read()
        if not (
            "Read assignment finished" in log_text or "Summary of counting results" in log_text
        ):
            raise FeatureCountsException(f"{log.stem} not complete")


def process_counts(counts: Path, output_file: str) -> None:
    df = parse_featureCounts_counts(counts, SRX)
    df.to_parquet(output_file)


def process_jcounts(jcounts: Path, output_file: Optional[str]) -> None:
    if not output_file:
        return

    df = parse_featureCounts_jcounts(jcounts, SRX)
    df.to_parquet(output_file)


class FeatureCountsException(Exception):
    """Fastq Screen Processing Exception"""


if __name__ == "__main__":
    try:
        main()
    except FeatureCountsException as error:
        logger.warning(f"{SRX}: {error}")
        LOG.append("Exception", text=error)

        raise SystemExit
    finally:
        remove_folder(TMPDIR)
