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
import os
from pathlib import Path

import pandas as pd

from ncbi_remap.parser import parse_hisat2


def main():
    df = parse_hisat2(snakemake.input[0])
    df.to_csv(snakemake.output[0], sep="\t", index=False)

    df.fillna(0, inplace=True)
    uniquely_aligned = (
        df.iloc[0, :]["num_concordant_reads_uniquely_aligned"]
        + df.iloc[0, :]["num_uniquely_aligned"]
    )

    if (df.iloc[0, :]["per_alignment"] < 1) | (uniquely_aligned < 1000):
        Path(Path(snakemake.input[0]).parent, "ALIGNMENT_BAD").touch()


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.mock import MockSnake

        SRX, SRR = "SRX133010", "SRR450563"  # SE
        # SRX, SRR = "SRX2921832", "SRR5687046"  # PE
        # SRX, SRR = "SRX2260227", "SRR4441419"  # keep R1
        # SRX, SRR = "SRX010901", "SRR027010"  # keep R2

        snakemake = MockSnake(
            input=f"output/aln-wf/samples/{SRX}/{SRR}/{SRR}.fq.bam.log",
            wildcards=dict(srx=SRX, srr=SRR),
        )

    main()
