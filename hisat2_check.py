"""Example Hisat2 output
# Example SE
83440 reads; of these:
  83440 (100.00%) were unpaired; of these:
    3605 (4.32%) aligned 0 times
    76302 (91.45%) aligned exactly 1 time
    3533 (4.23%) aligned >1 times
95.68% overall alignment rate

# Example PE
1338 reads; of these:
  1338 (100.00%) were paired; of these:
    468 (34.98%) aligned concordantly 0 times
    839 (62.71%) aligned concordantly exactly 1 time
    31 (2.32%) aligned concordantly >1 times
    ----
    468 pairs aligned concordantly 0 times; of these:
      22 (4.70%) aligned discordantly 1 time
    ----
    446 pairs aligned 0 times concordantly or discordantly; of these:
      892 mates make up the pairs; of these:
        807 (90.47%) aligned 0 times
        74 (8.30%) aligned exactly 1 time
        11 (1.23%) aligned >1 times
69.84% overall alignment rate
"""
import re
from pathlib import Path

import pandas as pd

from ncbi_remap.parser import parse_hisat2


def main(file_name):
    srx, srr = parse_wildcards(file_name)

    try:
        df = parse_hisat2(file_name)
        df.index = pd.MultiIndex.from_tuples([(srx, srr)], names=["srx", "srr"])
        df.to_csv(file_name.replace(".fq.bam.log", ".hisat2.bam.tsv"), sep="\t")

        df.fillna(0, inplace=True)
        uniquely_aligned = (
            df.iloc[0, :]["num_concordant_reads_uniquely_aligned"]
            + df.iloc[0, :]["num_uniquely_aligned"]
        )

        if (df.iloc[0, :]["per_alignment"] < 1) | (uniquely_aligned < 1000):
            Path(Path(file_name).parent, "ALIGNMENT_BAD").touch()
    except:
        global CNT
        CNT += 1
        print(file_name)


def parse_wildcards(file_name):
    srx = re.findall(r"[SED]RX\d+", file_name)[0]
    srr = re.findall(r"[SED]RR\d+", file_name)[0]
    return srx, srr


if __name__ == "__main__":
    CNT = 0
    for file_name in Path("output/prealn-wf/samples").glob("**/**/*.hisat2.bam.log"):
        main(file_name.as_posix())
    print(f"{CNT:,} problems")
