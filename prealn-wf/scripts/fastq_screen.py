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
from ncbi_remap.io import remove_folder
from ncbi_remap.snakemake import StepLogger
from ncbi_remap.parser import parse_fastq_screen


LOG = StepLogger(str(snakemake.log))
SRR = snakemake.wildcards.srr
THREADS = snakemake.threads
TMPDIR = Path(os.getenv("TMPDIR", "/tmp"), f"{SRR}/fastq_screen")
TMPDIR.mkdir(parents=True, exist_ok=True)
TMPREF = Path(os.getenv("TMPDIR", "/tmp"), f"references")
TMPREF.mkdir(exist_ok=True)
ALIGNER = "bowtie2"


def main():
    fastq, references = stage_data(snakemake.input[0], snakemake.params.references)
    config = create_config(references)
    results = fastq_screen(fastq, config)
    summarize(results, snakemake.output[0])


def stage_data(fastq_file: str, references: dict) -> Tuple[Path, dict]:
    # Copy FASTQ
    fastq_local = TMPDIR / f"{SRR}.fastq.gz"
    shutil.copy2(fastq_file, fastq_local)

    # Copy References
    references_local = {}
    for ref_name, ref_prefix in references.items():
        ref_dir_local = TMPREF / f"{ALIGNER}/{ref_name}"
        references_local[ref_name] = (ref_dir_local / Path(ref_prefix).name).as_posix()

        # Stage references if not present
        if ref_dir_local.exists():
            # Already on scratch wait to make sure copied
            sleep(15)
        else:
            shutil.copytree(Path(ref_prefix).parent, ref_dir_local, ignore=shutil.ignore_patterns("*.fasta"))

    return fastq_local, references_local


def create_config(references: dict) -> Path:
    config_file = TMPDIR / "fastq_screen_config.txt"

    with config_file.open("w") as fout:
        for ref_name, ref_dir in references.items():
            fout.write("\t".join(["DATABASE", ref_name, ref_dir, ALIGNER.upper()]) + "\n")
    LOG.append("Fastq Screen Config", config_file)
    return config_file


def fastq_screen(fastq: Path, config: Path) -> Path:
    log = TMPDIR / "fastq_screen.log"
    shell(
        f"fastq_screen --outdir {TMPDIR} "
        "--force "
        f"--aligner {ALIGNER} "
        f"--conf {config} "
        "--subset 100000 "
        f"--threads {THREADS} "
        f"{fastq} "
        f"> {log} 2>&1 "
    )
    LOG.append("Fastq Screen", log)

    # Make sure processing completed
    with log.open() as fh:
        if "Processing complete" not in fh.read():
            raise FastqScreenException

    return TMPDIR / f"{SRR}_screen.txt"


def summarize(results: Path, output_file: str) -> None:
    """Summarizes fastq screen results
        Cacluates the number of reads mapping to each specific reference. Ignores reads the map to multiple references.

        Goes from this:
        | reference   |   multiple_hits_multiple_libraries_count |   multiple_hits_multiple_libraries_percent |   multiple_hits_one_library_count |   multiple_hits_one_library_percent |   one_hit_multiple_libraries_count |   one_hit_multiple_libraries_percent |   one_hit_one_library_count |   one_hit_one_library_percent |   reads_processed_count |   unmapped_count |   unmapped_percent |
        |:------------|-----------------------------------------:|-------------------------------------------:|----------------------------------:|------------------------------------:|-----------------------------------:|-------------------------------------:|----------------------------:|------------------------------:|------------------------:|-----------------:|-------------------:|
        | adapters    |                                       48 |                                       0.05 |                                 0 |                                0    |                                  0 |                                 0    |                           0 |                          0    |                   99973 |            99925 |              99.95 |
        | dm6         |                                     1713 |                                       1.71 |                              6278 |                                6.28 |                                224 |                                 0.22 |                       88393 |                         88.42 |                   99973 |             3365 |               3.37 |
        | ecoli       |                                        1 |                                       0    |                                 0 |                                0    |                                  0 |                                 0    |                           2 |                          0    |                   99973 |            99970 |             100    |
        ...

        To this:
        |            |   adapters_pct_reads_mapped |   dm6_pct_reads_mapped |   ecoli_pct_reads_mapped |   ercc_pct_reads_mapped |   hg19_pct_reads_mapped |   phix_pct_reads_mapped |   rRNA_pct_reads_mapped |   wolbachia_pct_reads_mapped |   yeast_pct_reads_mapped |
        |:-----------|----------------------------:|-----------------------:|-------------------------:|------------------------:|------------------------:|------------------------:|------------------------:|-----------------------------:|-------------------------:|
        | SRR0000001 |                           0 |                94.6966 |               0.00200054 |                       0 |               0.0160043 |                       0 |              0.00100027 |                            0 |               0.00500135 |

    """
    df = parse_fastq_screen(results).set_index("reference").fillna(0)
    summarized = (
        (
            (df.one_hit_one_library_count + df.multiple_hits_one_library_count)
            / df.reads_processed_count
            * 100
        )
        .rename(SRR)
        .rename_axis("")
        .to_frame()
        .T.rename_axis("srr")
    )

    summarized.columns = [f"{col}_pct_reads_mapped" for col in summarized.columns]
    summarized.to_parquet(output_file)


class FastqScreenException(Exception):
    """Fastq Screen Processing Exception"""


if __name__ == "__main__":
    try:
        main()
    except FastqScreenException:
        logger.warning(f"{SRR}: fastq screen did not complete")
        LOG.append("Exception", text="fastq screen did not complete")
    finally:
        remove_folder(TMPDIR)
