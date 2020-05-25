"""Example Atropos output
# Example of SE
Reads                                 records  fraction
----------------------------------- --------- ---------
Total reads processed:                 85,319
Reads with adapter:                     4,725      5.5%
Reads that were too short:              1,879      2.2%
Reads written (passing filters):       83,440     97.8%

# Example of PE
Pairs                                records fraction
----------------------------------- -------- --------
Total read pairs processed:            1,354
  Read 1 with adapter:                   106     7.8%
  Read 2 with adapter:                   123     9.1%
Pairs that were too short:                16     1.2%
Pairs written (passing filters):       1,338    98.8%
"""
import os
import re
from pathlib import Path

import pandas as pd


class AtroposException(Exception):
    """Basic exception when there are problems running Atropos"""


def main():
    try:
        df = pd.DataFrame(
            [parse_atropos(snakemake.input.log)],
            index=pd.Index([snakemake.wildcards.srr], name="srr"),
            columns=["total_processed", "total_written", "too_short"],
        )
        df.to_parquet(snakemake.output[0])

        if df.total_written[0] < 1000:
            raise AtroposException

    except AtroposException:
        # Add to bad list
        atropos_bad_path = Path(Path(snakemake.output[0]).parents[3], "atropos_bad")
        atropos_bad_path.mkdir(exist_ok=True)
        Path(atropos_bad_path, snakemake.wildcards.srr).touch()

        # Delete atropos data if problems
        unlink(snakemake.input.R1)
        unlink(snakemake.input.R2)
        unlink(snakemake.input.log)


def parse_atropos(file_name):
    with open(file_name) as fh:
        string = fh.read().replace(",", "")

    if "ERROR" in string:
        raise AtroposException

    try:
        tot_processed = int(re.findall(r"Total read.*processed:\s+(\d+)", string)[0])
        tot_written = int(
            re.findall(r"[Read|Pair]s written \(passing filters\):\s+(\d+)", string)[0]
        )
        too_short = int(re.findall(r"[Read|Pair]s that were too short:\s+(\d+)", string)[0])
        return tot_processed, tot_written, too_short
    except IndexError:
        raise AtroposException


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
            input=dict(log=f"output/aln-wf/samples/{SRX}/{SRR}/{SRR}_1.trim.clean.fastq.gz.log"),
            wildcards=dict(srx=SRX, srr=SRR),
        )

    main()
