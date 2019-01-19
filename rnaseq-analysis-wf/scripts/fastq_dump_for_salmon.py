"""Merge together fastqs from multiple SRRs.

Given an SRX accession, this script concatenates FASTQs from all of the SRRs. If the SRR was determined to be PE
in the prealn-wf than R1 and R2 FASTQs will be created. Otherwise only a single read FASTQ will be output.

"""
from pathlib import Path
from tempfile import TemporaryDirectory, mkdtemp
from collections import namedtuple

import pandas as pd
from snakemake import shell


SRX = snakemake.wildcards.srx
LOG = snakemake.log_fmt_shell()
OUTPUT = snakemake.output[0]
TMPDIR = TemporaryDirectory().name

FastqPair = namedtuple('FastqPair', 'R1 R2')


def main():
    srrs = get_list_of_srrs()

    fnames = [
        dump_fastq(srr)
        for srr in srrs
    ]

    R1 = [fname.R1 for fname in fnames]
    R2 = [fname.R2 for fname in fnames if fname.R2 is not None]

    # Concatenate Read 1
    concat_fastqs(R1, read_number=1)

    # Concatenate Read 2 if same number of files as Read 1
    if len(R1) == len(R2):
        concat_fastqs(R2, read_number=2)


def get_list_of_srrs():
    store = pd.HDFStore('../output/sra.h5', mode='r')
    srrs = store['ids'][store['ids'].srx == SRX].srr.values.tolist()
    store.close()
    return srrs


def dump_fastq(srr: str):
    shell("fastq-dump -O {TMPDIR} -M 0 --split-files {srr} --gzip {LOG}")
    R1 = str(Path(TMPDIR, f'{srr}_1.fastq.gz'))
    R2 = str(Path(TMPDIR, f'{srr}_2.fastq.gz'))

    layout: str = get_layout(srr)
    if layout == 'SE' or layout == 'keep_R1':
        R2 = None

    if layout == 'keep_R2':
        R1 = R2
        R2 = None

    return FastqPair(R1, R2)


def get_layout(srr: str):
    with open(f'../output/prealn-wf/samples/{SRX}/{srr}/LAYOUT') as fh:
        return fh.read().strip()


def concat_fastqs(file_list: list, read_number: int = 1):
    if not file_list:
        return

    oname = OUTPUT
    if read_number == 2:
        oname = oname.replace('_1.fastq.gz', '_2.fastq.gz')

    file_names = ' '.join(file_list)
    shell(f'cat {file_names} > {oname}')


if __name__ == '__main__':
    main()
