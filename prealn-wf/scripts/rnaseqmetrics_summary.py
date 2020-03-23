import sys

import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_picardCollect_summary, parse_picardCollect_hist


def main():
    idx = pd.Index([snakemake.wildcards.srr], name="srr")

    # Parse the flags
    if parse_stranded(snakemake.input.first):
        strand = "same_strand"
    elif parse_stranded(snakemake.input.second):
        strand = "opposite_strand"
    else:
        strand = "unstranded"

    # Write strand flag
    df = pd.DataFrame([[strand]], index=idx, columns=["strand"])
    df.to_parquet(snakemake.output.strand)

    # Parse main table
    df = parse_table(snakemake.input.unstranded)
    df.index = idx
    df.to_parquet(snakemake.output.table)

    # Parse genome coverage histogram
    df = parse_picardCollect_hist(snakemake.input.unstranded)
    df.index = idx
    df.to_parquet(snakemake.output.genebody_coverage)


def parse_stranded(file_name):
    return (parse_picardCollect_summary(file_name).PCT_CORRECT_STRAND_READS >= 0.75)[0]


def parse_table(file_name):
    df = parse_picardCollect_summary(file_name)[
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
    ].fillna(0.0)
    df.PCT_CODING_BASES = df.PCT_CODING_BASES * 100
    df.PCT_UTR_BASES = df.PCT_UTR_BASES * 100
    df.PCT_INTRONIC_BASES = df.PCT_INTRONIC_BASES * 100
    df.PCT_INTERGENIC_BASES = df.PCT_INTERGENIC_BASES * 100
    df.PCT_MRNA_BASES = df.PCT_MRNA_BASES * 100

    df.columns = [x.lower() for x in df.columns]
    df.columns = [x.replace("pct_", "percent_") for x in df.columns]
    return df


if __name__ == "__main__":
    main()
