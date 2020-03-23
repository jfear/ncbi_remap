import sys
from pathlib import Path
from collections import namedtuple

import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_bamtools_stats, parse_samtools_stats


PREALN_PATH = Path("../output/prealn-wf/samples")

Files = namedtuple("Files", "srr idx idxstats samtools bamtools output")


def main():
    for srr_pth in PREALN_PATH.glob("**/SRR*"):
        if not srr_pth.is_dir():
            continue
        srr = srr_pth.name
        files = Files(
            srr,
            pd.Index([srr], name="srr"),
            (srr_pth / f"{srr}.hisat2.bam.samtools.idxstats"),
            (srr_pth / f"{srr}.hisat2.bam.samtools.stats"),
            (srr_pth / f"{srr}.hisat2.bam.bamtools.stats"),
            f"../output/prealn-wf/aln_stats/{srr}.parquet",
        )

        if files.samtools.exists() and files.bamtools.exists():
            convert_table(files)
            files.samtools.unlink()
            files.bamtools.unlink()

        if files.idxstats.exists():
            files.idxstats.unlink()


def convert_table(files):
    try:
        df = pd.concat([_samtools(files.samtools), _bamtools(files.bamtools)], axis=1, sort=False)
        df.index = pd.Index([files.srr], name="srr")
        df.to_parquet(files.output)
    except Exception as e:
        print(files.samtools)
        raise e


def _samtools(file_name):
    return parse_samtools_stats(file_name)[
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


def _bamtools(file_name):
    return parse_bamtools_stats(file_name)[["Percent Forward", "Percent Reverse"]]


if __name__ == "__main__":
    main()
