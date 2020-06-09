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


LOG = StepLogger(str(snakemake.log))
SRR = snakemake.wildcards.srr
THREADS = snakemake.threads

TMPDIR = Path(os.getenv("TMPDIR", "/tmp"), f"{SRR}/featurecounts")
TMPDIR.mkdir(parents=True, exist_ok=True)

TMPREF = Path(os.getenv("TMPDIR", "/tmp"), f"references")
TMPREF.mkdir(exist_ok=True)


def main():
    bam, gtf = stage_data(snakemake.input.bam, snakemake.input.gtf)
    counts, jcounts = featurecount(bam, gtf, snakemake.input.layout, snakemake.input.strand)
    summarize(counts, jcounts, snakemake.output[0])


def stage_data(bam: str, gtf: str) -> Tuple[Path, Path]:
    # Copy Bam
    bam_local = TMPDIR / f"{SRR}.bam"
    shutil.copy2(bam, bam_local)

    # Stage GTF if not present
    gtf_local = TMPREF / "dmel.gtf"
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
        params = "-p -P -C -J "
    else:
        params = "-J "

    # Look up strand
    strand_ = pd.read_parquet(strand).strand[0]
    if strand_ == "same_strand":
        params += "-s 1"
    elif strand_ == "opposite_strand":
        params += "-s 2"
    else:
        params += "-s 0"

    counts = TMPDIR / f"{SRR}.counts"
    jcounts = TMPDIR / f"{SRR}.counts.jcounts"
    summary = TMPDIR / f"{SRR}.counts.summary"
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


def summarize(counts: Path, jcounts: Path, output_file: str) -> None:
    gene_counts = _get_counts(counts)
    genic_reads = gene_counts.sum()
    percent_genes_on = (gene_counts > 0).mean() * 100

    junction_counts = _get_counts(jcounts)
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
        index=pd.Index([SRR], name="srr"),
    )

    df.to_parquet(output_file)


def _get_counts(file_name: Path) -> pd.DataFrame:
    col_name = pd.read_table(file_name, comment="#", nrows=1).columns[-1]
    return pd.read_table(file_name, comment="#", usecols=[col_name], dtype=np.int64).values


class FeatureCountsException(Exception):
    """Fastq Screen Processing Exception"""


if __name__ == "__main__":
    try:
        main()
    except FeatureCountsException as error:
        logger.warning(f"{SRR}: {error}")
        LOG.append("Exception", text=error)

        raise SystemExit
    finally:
        remove_folder(TMPDIR)
