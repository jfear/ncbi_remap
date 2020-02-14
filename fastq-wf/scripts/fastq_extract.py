"""Extract SRA file into Gziped FASTQs"""
import os
import gzip
import shutil
from pathlib import Path
from multiprocessing import Pool

from snakemake.shell import shell

TMPDIR = os.getenv("TMPDIR", "/tmp")


def main():
    r1_tmp = f"{TMPDIR}/{snakemake.wildcards.srr}_1.fastq"
    r2_tmp = f"{TMPDIR}/{snakemake.wildcards.srr}_2.fastq"
    extract_fastq()
    rename_files(r1_tmp, r2_tmp)

    with Pool(2) as pool:
        pool.map(gzip_file, [(r1_tmp, snakemake.output.r1), (r2_tmp, snakemake.output.r2)])

    clean_up(r1_tmp, r2_tmp)


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


def rename_files(r1, r2):
    if Path(TMPDIR, f"{snakemake.wildcards.srr}.sra.fastq").exists():
        # Single ended
        shutil.move(f"{TMPDIR}/{snakemake.wildcards.srr}.sra.fastq", r1)
        Path(r2).touch()  # always out R2 even for SE.
    elif Path(TMPDIR, f"{snakemake.wildcards.srr}.sra_1.fastq").exists():
        # Single ended
        shutil.move(f"{TMPDIR}/{snakemake.wildcards.srr}.sra_1.fastq", r1)
        shutil.move(f"{TMPDIR}/{snakemake.wildcards.srr}.sra_2.fastq", r2)


def gzip_file(files):
    file_in, file_out = files
    with open(file_in, "rb") as f_in:
        with gzip.open(file_out, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def clean_up(r1, r2):
    Path(r1).unlink()
    Path(r2).unlink()


if __name__ == "__main__":
    main()
