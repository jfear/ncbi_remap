"""Extract SRA file into Gziped FASTQs"""
import os
import gzip
import shutil
from pathlib import Path
from multiprocessing import Pool
from subprocess import CalledProcessError

from snakemake.shell import shell

TMPDIR = os.getenv("TMPDIR", "/tmp")


def main():
    staged_data = stage_data()
    extract_fastq(staged_data)

    for file_name in Path(TMPDIR).glob(f"{snakemake.wildcards.srr}*.fastq"):
        if "_1" in file_name.as_posix():
            shutil.move(file_name.absolute().as_posix(), snakemake.output.r1)
        elif "_2" in file_name.as_posix():
            shutil.move(file_name.absolute().as_posix(), snakemake.output.r2)
        else:
            # Single-ended read
            shutil.move(file_name.absolute().as_posix(), snakemake.output.r1)
            Path(snakemake.output.r2).touch()  # always out R2 even for SE.


def stage_data() -> str:
    """Copy data over to local scratch space."""
    data_pth = Path(snakemake.input[0])
    tmp_pth = Path(TMPDIR, data_pth.name)
    shutil.copy2(data_pth, tmp_pth)
    return tmp_pth.as_posix()


def extract_fastq(sra_file):
    try:
        shell(
            f"fasterq-dump {sra_file} "
            "--split-3 "
            "--skip-technical "
            "--min-read-len 20 "
            f"-O {TMPDIR} "
            f"-t {TMPDIR} "
            f"-e {snakemake.threads} "
        )
    except CalledProcessError:
        shell(
            f"fastq-dump {sra_file} "
            "--split-3 "
            "--skip-technical "
            "--minReadLen 20 "
            f"-O {TMPDIR} "
        )


if __name__ == "__main__":
    main()
