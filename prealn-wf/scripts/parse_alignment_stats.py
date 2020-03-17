import sys

import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_bamtools_stats, parse_samtools_stats


def main():
    df = pd.concat([_samtools(), _bamtools()], axis=1, sort=False)
    df.index = [snakemake.wildcards.srr]
    df.to_parquet(snakemake.output[0])


def _samtools():
    return parse_samtools_stats(snakemake.input.samtools)[
        [
            "reads_MQ0",
            "average_quality",
            "insert_size_average",
            "insert_size_standard_deviation",
            "pairs_on_different_chromosomes",
            "percentage_of_properly_paired_reads_(%)",
        ]
    ]


def _bamtools():
    return parse_bamtools_stats(snakemake.input.bamtools)[["Percent Forward", "Percent Reverse"]]


if __name__ == "__main__":
    main()
