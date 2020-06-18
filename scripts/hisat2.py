import os
import shutil
import sys
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
from ncbi_remap.parser import parse_hisat2, parse_bamtools_stats, parse_samtools_stats


LOG = StepLogger(str(snakemake.log))
SRA = snakemake.wildcards.get("srr", snakemake.wildcards.get("srx"))
SRA_TYPE = "srr" if "RR" in SRA else "srx"
THREADS = snakemake.threads
MEM = int(snakemake.resources.get("mem_gb", 4)) - 1

TMPDIR = Path(os.getenv("TMPDIR", "/tmp"), f"{SRA}/hisat2")
TMPDIR.mkdir(parents=True, exist_ok=True)

TMPREF = Path(os.getenv("TMPDIR", "/tmp"), f"references")
TMPREF.mkdir(exist_ok=True)


def main():
    r1, r2, reference = stage_data(
        snakemake.input.r1, snakemake.input.r2, snakemake.params.reference
    )

    log, sam = hisat2(
        snakemake.input.layout,
        r1,
        r2,
        reference,
        snakemake.input.get("strand", None),
        snakemake.params.get("splice_sites", None),
    )
    check_hisat(log, snakemake.output.hisat_summary)

    bam, bai = compress_sort_and_index(sam)
    save_output(bam, bai, snakemake.output.bam, snakemake.output.bai)

    alignment_stats(bam, snakemake.output.aln_stats)


def stage_data(r1: str, r2: str, reference: str) -> Tuple[Path, Path, str]:
    # Copy FASTQ
    r1_local = TMPDIR / f"{SRA}_1.fastq.gz"
    shutil.copy2(r1, r1_local)

    r2_local = TMPDIR / f"{SRA}_2.fastq.gz"
    shutil.copy2(r2, r2_local)

    # Copy Reference
    ref_dir_local = TMPREF / "hisat2/dmel"
    ref_prefix_local = (ref_dir_local / Path(reference).name).as_posix()

    # Stage references if not present
    if ref_dir_local.exists():
        # Already on scratch wait to make sure copied
        sleep(15)
    else:
        shutil.copytree(
            Path(reference).parent, ref_dir_local, ignore=shutil.ignore_patterns("*.fasta")
        )

    return r1_local, r2_local, ref_prefix_local


def hisat2(
    layout: str,
    r1: Path,
    r2: Path,
    reference: str,
    strand: Optional[str] = None,
    splice: Optional[str] = None,
) -> Tuple[Path, Path]:
    # Look up Layout
    layout_ = pd.read_parquet(layout).layout[0]
    if layout == "PE":
        fastqs = f"-1 {r1} -2 {r2}"
    elif layout == "keep_R2":
        fastqs = f"-U {r2}"
    else:
        fastqs = f"-U {r1}"

    # Look up strand if it is there
    strand_param = ""
    if strand:
        strand_ = pd.read_parquet(strand).strand[0]
        if (layout == "PE") & (strand_ == "first_strand"):
            strand_param = "--rna-strandness FR"
        elif (layout == "PE") & (strand_ == "second_strand"):
            strand_param = "--rna-strandness RF"
        elif strand_ == "first_strand":
            strand_param = "--rna-strandness F"
        elif strand_ == "second_strand":
            strand_param = "--rna-strandness R"

    # Use known splice sites
    splice_param = ""
    if splice:
        splice_param = f"--known-splicesite-infile {splice}"

    # Run Hisat2
    log = TMPDIR / "hisat2.log"
    sam = TMPDIR / f"{SRA}.sam"
    cmd = (
        "hisat2 "
        f"-x {reference} "
        f"{fastqs} "
        f"--threads {THREADS} "
        "--max-intronlen 300000 "
        f"{strand_param} "
        f"{splice_param} "
        f"-S {sam} "
        f">{log} 2>&1"
    )
    LOG.append("Hisat2 Command", text=cmd)

    try:
        shell(cmd)
    finally:
        LOG.append("Hisat2", log)

    remove_file(r1)
    remove_file(r2)

    return log, sam


def check_hisat(log: Path, output_file: str) -> None:
    """Example Hisat2 output
        # Example SE
        83440 reads; of these:
        83440 (100.00%) were unpaired; of these:
            3605 (4.32%) aligned 0 times
            76302 (91.45%) aligned exactly 1 time
            3533 (4.23%) aligned >1 times
        95.68% overall alignment rate

        # Example PE
        1338 reads; of these:
        1338 (100.00%) were paired; of these:
            468 (34.98%) aligned concordantly 0 times
            839 (62.71%) aligned concordantly exactly 1 time
            31 (2.32%) aligned concordantly >1 times
            ----
            468 pairs aligned concordantly 0 times; of these:
            22 (4.70%) aligned discordantly 1 time
            ----
            446 pairs aligned 0 times concordantly or discordantly; of these:
            892 mates make up the pairs; of these:
                807 (90.47%) aligned 0 times
                74 (8.30%) aligned exactly 1 time
                11 (1.23%) aligned >1 times
        69.84% overall alignment rate
    """
    df = parse_hisat2(log).fillna(0)
    df.index = pd.Index([SRA], name=SRA_TYPE)
    df.to_parquet(output_file)
    remove_file(log)

    uniquely_aligned = (
        df.iloc[0, :]["num_concordant_reads_uniquely_aligned"]
        + df.iloc[0, :]["num_uniquely_aligned"]
    )

    per_aligned = df.iloc[0, :]["per_alignment"]

    if (per_aligned < 1) | (uniquely_aligned < 1000):
        raise HisatException(f"Poor alignment: {uniquely_aligned:,} ({per_aligned}%)")


def compress_sort_and_index(sam: Path) -> Tuple[Path, Path]:
    # Convert to BAM
    bam = TMPDIR / f"{SRA}.bam"
    log = TMPDIR / "sam_to_bam.log"
    try:
        shell(f"samtools view -Sb -q20 --threads {THREADS} {sam} > {bam} 2>{log}")
    finally:
        LOG.append("Sam to Bam", log)
    remove_file(sam)
    remove_file(log)

    # Sort BAM
    sorted_bam = TMPDIR / f"{SRA}.sorted.bam"
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
    df.index = pd.Index([SRA], name=SRA_TYPE)
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


class HisatException(Exception):
    """Hisat2 Processing Exception"""


if __name__ == "__main__":
    try:
        main()
    except HisatException as error:
        logger.warning(f"Flagging {SRA} as Hisat Bad: {error}")
        LOG.append("Exception", text=str(error))

        # Add flag
        pth = Path(snakemake.params.alignment_bad)
        pth.mkdir(exist_ok=True)
        (pth / SRR).touch()

        # Remove outputs
        remove_folder(Path(snakemake.output.bam).parent)
        remove_file(snakemake.output.hisat_summary)
        remove_file(snakemake.output.aln_stats)
        
        raise SystemExit
    finally:
        remove_folder(TMPDIR)
