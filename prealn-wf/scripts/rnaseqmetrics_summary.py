import sys

import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_picardCollect_summary, parse_picardCollect_hist
from ncbi_remap.snakemake import put_flag


def main():
    # Parse the flags
    if parse_stranded(snakemake.input.first):
        strand = "same_strand"
    elif parse_stranded(snakemake.input.second):
        strand = "same_strand"
    else:
        strand = "unstranded"

    put_flag(snakemake.output.flag, strand)

    # Parse main table
    df = parse_unstranded()
    df["strand"] = strand
    df.index = [snakemake.wildcards.srr]
    df.to_parquet(snakemake.output.table)

    # Parse genome coverage histogram
    df = parse_picardCollect_hist(snakemake.input.unstranded)
    df.index = [snakemake.wildcards.srr]
    df.to_parquet(snakemake.output.gene_coverage)


def parse_stranded(file_name):
    parse_picardCollect_summary(file_name).PCT_CORRECT_STRAND_READS >= 0.75


def parse_unstranded():
    df = parse_picardCollect_summary(snakemake.input.unstranded)[
        [
            "PCT_CODING_BASES",
            "PCT_UTR_BASES",
            "PCT_INTRONIC_BASES",
            "PCT_INTERGENIC_BASES",
            "PCT_MRNA_BASES",
            "MEDIAN_CV_COVERAGE",
            "MEDIAN_5PRIME_BIAS",
            "MEDIAN_3PRIME_BIAS",
            "MEDIAN_5PRIME_TO_3PRIME_BIAS",
        ]
    ]
    df.columns = [x.lower() for x in df.columns]
    df.columns = [x.replace("pct_", "percent_") for x in df.columns]
    return df


if __name__ == "__main__":
    main()
