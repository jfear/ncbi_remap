import os
from pathlib import Path
import snakemake
import pandas as pd

#SRX = snakemake.wildcards.srx
SRX = 'SRX008017'
os.chdir('../')


def main():
    srrs = get_list_of_srrs()
    fnames = []
    for srr in srrs:
        dump_fastq(srr)


def get_list_of_srrs():
    store = pd.HDFStore('../output/sra.h5', mode='r')
    srrs = store['ids'][store['ids'].srx == "SRX008017"].srr.values.tolist()
    store.close()
    return srrs


def dump_fastq(srr):
    snakemake.shell("""
    fastq-dump -
    """)


if __name__ == '__main__':
    main()
