import sys

import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_bamtools_stats, parse_samtools_stats


def main():
    df = pd.concat([_samtools(), _bamtools()], axis=1, sort=False)
    df.index = pd.Index([snakemake.wildcards.srr], name="srr")
    df.to_parquet(snakemake.output[0])


def _samtools():
    return parse_samtools_stats(snakemake.input.samtools)[
        [
            "reads_MQ0",
            "average_quality",
            "insert_size_average",
            "insert_size_standard_deviation",
            "inward_oriented_pairs",
            "outward_oriented_pairs",
            "pairs_with_other_orientation",
            "pairs_on_different_chromosomes",
        ]
    ]


def _bamtools():
    return parse_bamtools_stats(snakemake.input.bamtools)[["Percent Forward", "Percent Reverse"]]


if __name__ == "__main__":
    main()
