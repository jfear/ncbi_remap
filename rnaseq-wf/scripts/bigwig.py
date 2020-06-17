import os
import shutil
import ssl
import sys
from pathlib import Path
from time import sleep
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import pybedtools
from snakemake.logging import logger
from snakemake.shell import shell

sys.path.insert(0, "../src")
from ncbi_remap.io import remove_file, remove_folder
from ncbi_remap.snakemake import StepLogger


LOG = StepLogger(str(snakemake.log))
SRX = snakemake.wildcards.srx
CHROM_SOURCE = snakemake.params.get("chrom_source", "ucsc")
THREADS = snakemake.threads

TMPDIR = Path(os.getenv("TMPDIR", "/tmp"), f"{SRX}/bigwig")
TMPDIR.mkdir(parents=True, exist_ok=True)

TMPREF = Path(os.getenv("TMPDIR", "/tmp"), f"references")
TMPREF.mkdir(exist_ok=True)


def main():
    bam, bai, chrom_sizes = stage_data(
        snakemake.input.bam, snakemake.input.bai, snakemake.input.chromSizes
    )

    # BedGraph
    first_bg, second_bg = run_bam_coverage(bam, bai)
    if CHROM_SOURCE == "flybase":
        convert_ucsc_to_flybase(first_bg, second_bg)

    # BigWig
    first_bw, second_bw = convert_bedgraph_to_bigwig(first_bg, second_bg, chrom_sizes)
    save_output(first_bw, second_bw, snakemake.output.first, snakemake.output.second)


def stage_data(bam: str, bai: str, chrom_sizes: str) -> Tuple[Path, Path, Path]:
    # Copy Bam
    bam_local = TMPDIR / f"{SRX}.bam"
    shutil.copy2(bam, bam_local)

    # Copy Bai
    bai_local = TMPDIR / f"{SRX}.bam.bai"
    shutil.copy2(bai, bai_local)

    # Stage Chrom Sizes if not present
    sizes_local = TMPREF / Path(chrom_sizes).name
    if sizes_local.exists():
        # Already on scratch wait to make sure copied
        sleep(5)
    else:
        shutil.copy(chrom_sizes, sizes_local)

    return bam_local, bai_local, sizes_local


def run_bam_coverage(bam: Path, bai: Path) -> Tuple[Path, Path]:
    try:
        first = _bam_coverage(bam, "first")
        second = _bam_coverage(bam, "second")
    finally:
        remove_file(bam)
        remove_file(bai)
    return first, second


def convert_ucsc_to_flybase(first: Path, second: Path) -> None:
    # I was getting an SSL error.
    # From: https://stackoverflow.com/a/56230607/4605992
    ssl._create_default_https_context = ssl._create_unverified_context

    # Temp file for conversion
    first_tmp = first.with_suffix(".tmp.bg")
    second_tmp = second.with_suffix(".tmp.bg")

    # Run Coversion
    conversion_table = _import_conversion()
    _pybedtools_convert(first, first_tmp, conversion_table)
    _pybedtools_convert(second, second_tmp, conversion_table)

    # Overwrite original file with converted temp
    shutil.move(first_tmp, first)
    shutil.move(second_tmp, second)


def convert_bedgraph_to_bigwig(first: Path, second: Path, chrom_sizes: Path) -> Tuple[Path, Path]:
    first_bw = _bedgraph_to_bigwig(first, chrom_sizes)
    second_bw = _bedgraph_to_bigwig(second, chrom_sizes)
    return first_bw, second_bw


def save_output(first: Path, second: Path, first_out: str, second_out: str) -> None:
    shutil.copy2(first, first_out)
    shutil.copy2(second, second_out)


def _bam_coverage(bam: Path, strand: str) -> Path:
    mapper = {"first": "forward", "second": "reverse"}
    bedgraph = TMPDIR / f"{SRX}.{strand}.bg"
    log = TMPDIR / f"bam_coverage_{strand}.log"
    try:
        shell(
            "bamCoverage "
            "--outFileFormat bedgraph "
            "--binSize 1 "
            "--effectiveGenomeSize 129000000 "
            "--normalizeUsing RPGC "
            "--ignoreForNormalization chrX "
            f"-p {THREADS} --filterRNAstrand {mapper[strand]} "
            f"--bam {bam} -o {bedgraph} > {log} 2>&1"
        )
    finally:
        LOG.append(f"BedGraph {strand.title()} Strand", log)
        remove_file(log)
    return bedgraph


def _import_conversion() -> dict:
    """Download conversion table from NCBI.

        col#   Header
        0      Sequence-Name
        1      Sequence-Role
        2      Assigned-Molecule
        3      Assigned-Molecule-Location/Type
        4      GenBank-Accn
        5      Relationship
        6      RefSeq-Accn
        7      Assembly-Unit
        8      Sequence-Length
        9      UCSC-style-name
    """
    file_name = TMPREF / "chrom_conversion.tsv"
    if file_name.exists():
        df = pd.read_table(file_name, header=None)
    else:
        url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_assembly_report.txt"
        # This table has a small error were FlyBase mitochondrion_genome is MT
        df = pd.read_table(url, comment="#", header=None).replace("MT", "mitochondrion_genome")
        df.to_csv(file_name, sep="\t", header=False, index=False)

    mapping = {"FlyBase": 0, "UCSC": 9, "GenBank": 4, "RefSeq": 6}
    return {k: v for k, v in df[[mapping["UCSC"], mapping["FlyBase"]]].values}


def _pybedtools_convert(bedgraph: Path, output_file: Path, mapper: dict) -> None:
    """ Use pybedtools to convert chromosomes from UCSC (chrX) to FlyBase format (X)"""
    bt = pybedtools.BedTool(bedgraph)
    bt.each(_convertFeature, mapper).saveas(output_file)


def _convertFeature(f: pybedtools.Interval, mapper: dict) -> pybedtools.Interval:
    f.chrom = mapper[f.chrom]
    return f


def _bedgraph_to_bigwig(bedgraph: Path, chrom_sizes: Path) -> Tuple[Path, Path]:
    try:
        sorted_bedgraph = bedgraph.with_suffix(".sorted.bg")
        bigwig = bedgraph.with_suffix(".bw")
        log = bedgraph.with_suffix(".bw.log")
        shell(
            f"LC_COLLATE=C bedSort {bedgraph} {sorted_bedgraph} > {log} 2>&1 && "
            f"bedGraphToBigWig {sorted_bedgraph} {chrom_sizes} {bigwig} >> {log} 2>&1"
        )
    finally:
        LOG.append(f"Convert {bedgraph.name} to BigWig", log)
        remove_file(log)
        remove_file(bedgraph)
        remove_file(sorted_bedgraph)
    return bigwig


if __name__ == "__main__":
    try:
        main()
    finally:
        remove_folder(TMPDIR)
