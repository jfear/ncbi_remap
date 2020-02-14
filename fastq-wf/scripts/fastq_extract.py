"""Extract SRA file into Gziped FASTQs"""
import os
import gzip
import shutil
from pathlib import Path
from multiprocessing import Pool

from snakemake.shell import shell

TMPDIR = os.getenv("TMPDIR", "/tmp")


def main():
    extract_fastq()

    if Path(TMPDIR, f"{snakemake.wildcards.srr}.sra.fastq").exists():
        # Single ended
        shutil.move(f"{TMPDIR}/{snakemake.wildcards.srr}.sra.fastq", snakemake.output.r1)
        Path(snakemake.output.r2).touch()  # always out R2 even for SE.
    elif Path(TMPDIR, f"{snakemake.wildcards.srr}.sra_1.fastq").exists():
        # Paired ended
        shutil.move(f"{TMPDIR}/{snakemake.wildcards.srr}.sra_1.fastq", snakemake.output.r1)
        shutil.move(f"{TMPDIR}/{snakemake.wildcards.srr}.sra_2.fastq", snakemake.output.r2)


def extract_fastq():
    shell(
        f"fasterq-dump {snakemake.input[0]} "
        "--split-3 "
        "--skip-technical "
        "--min-read-len 20 "
        f"-O {TMPDIR} "
        f"-t {TMPDIR} "
        f"-e {snakemake.threads} "
    )



if __name__ == "__main__":
    main()
