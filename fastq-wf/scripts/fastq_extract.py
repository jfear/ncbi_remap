"""Extract SRA file into Gziped FASTQs"""
import os
import gzip
import shutil
from pathlib import Path
from multiprocessing import Pool

from snakemake.shell import shell

TMPDIR = os.getenv("TMPDIR", "/tmp")


def main():
    staged_data = stage_data()
    extract_fastq(staged_data)

    if Path(TMPDIR, f"{snakemake.wildcards.srr}.sra.fastq").exists():
        # Single ended
        shutil.move(f"{TMPDIR}/{snakemake.wildcards.srr}.sra.fastq", snakemake.output.r1)
        Path(snakemake.output.r2).touch()  # always out R2 even for SE.
    elif Path(TMPDIR, f"{snakemake.wildcards.srr}.sra_1.fastq").exists():
        # Paired ended
        shutil.move(f"{TMPDIR}/{snakemake.wildcards.srr}.sra_1.fastq", snakemake.output.r1)
        shutil.move(f"{TMPDIR}/{snakemake.wildcards.srr}.sra_2.fastq", snakemake.output.r2)


def stage_data() -> str:
    """Copy data over to local scratch space."""
    data_pth = Path(snakemake.input[0])
    tmp_pth = Path(TMPDIR, data_pth.name)
    shutil.copy2(data_pth, tmp_pth)
    return tmp_pth.as_posix()


def extract_fastq(sra_file):
    shell(
        f"fasterq-dump {sra_file} "
        "--split-3 "
        "--skip-technical "
        "--min-read-len 20 "
        f"-O {TMPDIR} "
        f"-t {TMPDIR} "
        f"-e {snakemake.threads} "
    )



if __name__ == "__main__":
    main()
