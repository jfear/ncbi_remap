__author__ = "Justin Fear"
__copyright__ = "Copyright 2016, Justin Fear"
__email__ = "justin.fear@nih.gov"
__license__ = "MIT"

"""Example usage.

rule fastq_dump:
    output:
        fq1=patterns['fastq']['r1'],
        fq2=patterns['fastq']['r2'],
        flag=patterns['layout'],
        summary=patterns['fastq']['summary'],
    log: patterns['fastq']['r1'] + '.log'
    wrapper:
        wrapper_for('wrappers/fastq_dump')

"""

import os
import sys
import gzip
import shutil as sh
from pathlib import Path
from collections import namedtuple

import pandas as pd

from snakemake.shell import shell

sys.path.insert(0, "../src")
from ncbi_remap.fastq import check_fastq, md5sum, fastq_stats, fastq_abi_solid
from ncbi_remap.snakemake import put_flag

TMPDIR = os.environ["TMPDIR"]

ReadSummary = namedtuple("ReadSummary", "md5 size len abi")


def main():
    # Dump FASTQ to a TMPDIR
    shell(
        "fastq-dump -O {tmpdir} -M 0 --split-files {srr} {log}".format(
            tmpdir=TMPDIR, srr=snakemake.wildcards.srr, log=snakemake.log_fmt_shell()
        )
    )

    r1 = summarize_r1(os.path.join(TMPDIR, snakemake.wildcards.srr + "_1.fastq"))
    r2 = summarize_r2(os.path.join(TMPDIR, snakemake.wildcards.srr + "_2.fastq"))

    # Save summaries
    pd.DataFrame(
        [[r1.md5, r1.size, r1.len, r2.md5, r2.size, r2.len]],
        columns=["md5_R1", "libsize_R1", "avgLen_R1", "md5_R2", "libsize_R2", "avgLen_R2"],
    ).to_csv(snakemake.output.summary, sep="\t", index=False)

    create_flags(r1, r2)


def summarize_r1(file_name):
    md5 = md5sum(file_name)
    lib_size, avg_len = fastq_stats(file_name)
    abi = fastq_abi_solid(file_name)

    with open(file_name, "rb") as f_in:
        with gzip.open(snakemake.output.fq1, "wb") as f_out:
            sh.copyfileobj(f_in, f_out)

    return ReadSummary(md5, lib_size, avg_len, abi)


def summarize_r2(file_name):
    if check_fastq(file_name):
        # Get md5sum
        md5 = md5sum(file_name)
        lib_size, avg_len = fastq_stats(file_name)
        abi = fastq_abi_solid(file_name)

        with open(file_name, "rb") as f_in:
            with gzip.open(snakemake.output.fq2, "wb") as f_out:
                sh.copyfileobj(f_in, f_out)
    else:
        md5 = lib_size = avg_len = None
        abi = False
        Path(snakemake.output.fq2).touch()

    return ReadSummary(md5, lib_size, avg_len, abi)


def create_flags(r1, r2):
    if (r1.size > 1000) & (r2.size is None) & (r1.len > 10) & (r2.len is None):
        # R2 does not exists so SE
        put_flag(snakemake.output.flag, "SE")

    elif (r2.size is None) & (r2.len is None):
        # SE but R1 looks bad, flag as download bad
        fname = os.path.join(os.path.dirname(snakemake.output.fq1), "DOWNLOAD_BAD")
        Path(fname).touch()

    elif (r1.size > 1000) & (r2.size > 1000) & (r1.len > 10) & (r2.len > 10):
        if r1.size == r2.size:
            # Both reads look good so PE
            put_flag(snakemake.output.flag, "PE")
        else:
            # There are an uneven number of reads between R1 and R2.
            # Instead of messing with this, just consider SE and use R1.
            put_flag(snakemake.output.flag, "keep_R1")

    elif (r1.size > 1000) & (r1.len > 10):
        # Only R1 looks ok, so SE
        put_flag(snakemake.output.flag, "keep_R1")

    elif (r2.size > 1000) & (r2.len > 10):
        # Only R2 looks ok, consider single-end
        put_flag(snakemake.output.flag, "keep_R2")

    else:
        # R1 and R2 look bad, flag as download bad
        fname = os.path.join(os.path.dirname(snakemake.output.fq1), "DOWNLOAD_BAD")
        Path(fname).touch()

    if r1.abi | r2.abi:
        fname = os.path.join(os.path.dirname(snakemake.output.fq1), "ABI_SOLID")
        Path(fname).touch()


if __name__ == "__main__":

    main()
