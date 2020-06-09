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
import numpy as np

sys.path.insert(0, "../src")
from ncbi_remap.io import remove_folder, remove_file
from ncbi_remap.snakemake import StepLogger
from ncbi_remap.parser import parse_picard_markduplicate_metrics


LOG = StepLogger(str(snakemake.log))
SRR = snakemake.wildcards.srr
THREADS = snakemake.threads
MEM = int(snakemake.resources.get("mem_gb", 4))

TMPDIR = Path(os.getenv("TMPDIR", "/tmp"), f"{SRR}/markduplicates")
TMPDIR.mkdir(parents=True, exist_ok=True)


def main():
    bam = stage_data(snakemake.input[0])
    metrics = markduplicates(bam)
    summarize(metrics, snakemake.output[0])


def stage_data(bam: str) -> Path:
    bam_local = TMPDIR / f"{SRR}.bam"
    shutil.copy2(bam, bam_local)
    return bam_local


def markduplicates(bam: Path) -> Path:
    dedup_bam = TMPDIR / f"{SRR}.dedup.bam"
    metrics = TMPDIR / f"{SRR}.metrics"
    log = TMPDIR / f"markduplicates.log"

    try:
        shell(
            f"picard -Xmx{MEM}g MarkDuplicates INPUT={bam} "
            f"OUTPUT={dedup_bam} METRICS_FILE={metrics} >{log} 2>&1"
        )
    finally:
        LOG.append("Markduplicates", log)

    _check_log(log)
    remove_file(log)
    remove_file(dedup_bam)
    return metrics


def _check_log(log: Path) -> None:
    with log.open() as fh:
        if not "MarkDuplicates done" in fh.read():
            raise PicardException(f"{log.stem} not complete")


def summarize(metrics: Path, output_file: str):
    dtypes = {
        "UNPAIRED_READS_EXAMINED": np.int64,
        "READ_PAIRS_EXAMINED": np.int64,
        "UNPAIRED_READ_DUPLICATES": np.int64,
        "READ_PAIR_DUPLICATES": np.int64,
        "PERCENT_DUPLICATION": np.float64,
        "ESTIMATED_LIBRARY_SIZE": np.int64,
    }

    df = parse_picard_markduplicate_metrics(metrics)[dtypes.keys()].fillna(0).astype(dtypes)
    df.PERCENT_DUPLICATION = df.PERCENT_DUPLICATION * 100
    df.columns = [col.lower() for col in df.columns]
    df.index = pd.Index([SRR], name="srr")
    df.to_parquet(output_file)


class PicardException(Exception):
    """Picard Processing Exception"""


if __name__ == "__main__":
    try:
        main()
    except PicardException as error:
        logger.warning(f"{SRR}: {error}")
        LOG.append("Exception", text=str(error))

        raise SystemExit
    finally:
        remove_folder(TMPDIR)
