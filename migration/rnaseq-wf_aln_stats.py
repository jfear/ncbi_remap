import sys
from pathlib import Path
from collections import namedtuple

import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_bamtools_stats, parse_samtools_stats


RNASEQ_PATH = Path("../output/rnaseq-wf/samples")

Files = namedtuple("Files", "srx idx idxstats samtools bamtools output")


def main():
    for srx_pth in RNASEQ_PATH.iterdir():
        if not srx_pth.is_dir():
            continue
        srx = srx_pth.name
        files = Files(
            srx,
            pd.Index([srx], name="srx"),
            (srx_pth / f"{srx}.bam.samtools.idxstats"),
            (srx_pth / f"{srx}.bam.samtools.stats"),
            (srx_pth / f"{srx}.bam.bamtools.stats"),
            f"../output/rnaseq-wf/aln_stats/{srx}.parquet",
        )

        if files.samtools.exists() and files.bamtools.exists():
            if files.samtools.stat().st_size > 0 and files.bamtools.stat().st_size > 0:
                convert_table(files)
            else:
                print(files.samtools, files.bamtools)

            files.samtools.unlink()
            files.bamtools.unlink()

        if files.idxstats.exists():
            files.idxstats.unlink()


def convert_table(files):
    try:
        df = pd.concat([_samtools(files.samtools), _bamtools(files.bamtools)], axis=1, sort=False)
        df.index = pd.Index([files.srx], name="srx")
        df.to_parquet(files.output)
    except Exception as e:
        print(files.samtools)
        print(files.bamtools)
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
