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
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_hisat2


class Hisat2Exception(Exception):
    """Basic exception when there are problems running Hisat2"""


def main():
    try:
        df = parse_hisat2(snakemake.input.log).fillna(0)
        df.index = pd.Index([snakemake.wildcards.srr], name="srr")
        df.to_parquet(snakemake.output[0])

        uniquely_aligned = (
            df.iloc[0, :]["num_concordant_reads_uniquely_aligned"]
            + df.iloc[0, :]["num_uniquely_aligned"]
        )

        if (df.iloc[0, :]["per_alignment"] < 1) | (uniquely_aligned < 1000):
            raise Hisat2Exception

    except Hisat2Exception:
        alignment_bad_path = Path(Path(snakemake.input.log).parents[1], "alignment_bad")
        alignment_bad_path.mkdir(exist_ok=True)
        Path(alignment_bad_path, snakemake.wildcards.srr).touch()

        # Remove hisat2 output if problems
        unlink(snakemake.input.bam)
        unlink(snakemake.input.bai)
        unlink(snakemake.input.log)


def unlink(file_name):
    if Path(file_name).exists():
        Path(file_name).unlink()


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
