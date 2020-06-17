import os
import shutil
import sys
from pathlib import Path
from typing import Tuple, Iterable, List

from snakemake.logging import logger
from snakemake.shell import shell
import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.io import remove_folder, remove_file
from ncbi_remap.snakemake import StepLogger
from ncbi_remap.parser import parse_bamtools_stats, parse_samtools_stats


LOG = StepLogger(str(snakemake.log))
SRX = snakemake.wildcards.srx
THREADS = snakemake.threads
MEM = int(snakemake.resources.get("mem_gb", 4)) - 1

TMPDIR = Path(os.getenv("TMPDIR", "/tmp"), f"{SRX}/merge")
TMPDIR.mkdir(parents=True, exist_ok=True)


def main():
    bams = stage_data(snakemake.input)

    bam, bai = merge(bams)
    save_output(bam, bai, snakemake.output.bam, snakemake.output.bai)

    alignment_stats(bam, snakemake.output.aln_stats)


def stage_data(inputs: Iterable) -> List[Path]:
    bams = []
    for bam_file in inputs:
        bam_local = TMPDIR / Path(bam_file).name
        shutil.copy2(bam_file, bam_local)
        bams.append(bam_local)

    return bams


def merge(bams: List[Path]) -> Tuple[Path, Path]:
    bams_ = " ".join(bams)
    bam = TMPDIR / f"{SRX}.bam"
    log = TMPDIR / "merge.log"
    try:
        shell(f"samtools merge -f --threads {THREADS} {bam} {bams_} > {log} 2>&1")
    finally:
        LOG.append("Merge", log)
        [remove_file(bam_) for bam_ in bams]
        remove_file(log)

    # Sort BAM
    sorted_bam = TMPDIR / f"{SRX}.sorted.bam"
    log = TMPDIR / "bam_sort.log"
    tmp = TMPDIR / "samtools_sort"
    try:
        shell(
            f"samtools sort -l 9 -m {MEM}G --output-fmt BAM "
            f"-T {tmp} --threads {THREADS} -o {sorted_bam} {bam} 2> {log}"
        )
    finally:
        LOG.append("Sort Bam", log)
        remove_file(bam)
        remove_file(log)

    # Index BAM
    sorted_bai = TMPDIR / f"{SRA}.sorted.bam.bai"
    log = TMPDIR / "bam_index.log"
    try:
        shell(f"samtools index {sorted_bam} 2> {log}")
    finally:
        LOG.append("Index Bam", log)
        remove_file(log)

    return sorted_bam, sorted_bai


def save_output(bam: Path, bai: Path, bam_out: str, bai_out: str) -> None:
    shutil.copy2(bam, bam_out)
    shutil.copy2(bai, bai_out)


def alignment_stats(bam: Path, output_file: str) -> None:
    samtools_stats = TMPDIR / "samtools.stats"
    log = TMPDIR / "samtools.log"
    try:
        shell(f"samtools stats {bam} > {samtools_stats} 2>{log}")
    finally:
        LOG.append("Samtools Stats", log)
        remove_file(log)

    log = TMPDIR / "bamtools.log"
    bamtools_stats = TMPDIR / "bamtools.stats"
    try:
        shell(f"bamtools stats -in {bam} > {bamtools_stats} 2>{log}")
    finally:
        LOG.append("Bamtools Stats", log)
        remove_file(log)

    # Summarize
    df = pd.concat([_samtools(samtools_stats), _bamtools(bamtools_stats)], axis=1, sort=False)
    df.index = pd.Index([SRX], name="srx")
    df.to_parquet(output_file)


def _samtools(stats_file: Path) -> pd.DataFrame:
    return parse_samtools_stats(stats_file)[
        [
            "reads_MQ0",
            "average_quality",
            "insert_size_average",
            "insert_size_standard_deviation",
            "inward_oriented_pairs",
            "outward_oriented_pairs",
            "pairs_with_other_orientation",
            "pairs_on_different_chromosomes",
        ]
    ]


def _bamtools(stats_file: Path) -> pd.DataFrame:
    return parse_bamtools_stats(stats_file)[["Percent Forward", "Percent Reverse"]]


if __name__ == "__main__":
    try:
        main()
    finally:
        remove_folder(TMPDIR)
