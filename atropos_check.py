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


def main(file_name):
    srx, srr = parse_wildcards(file_name)

    try:
        df = pd.DataFrame(
            [parse_atropos(file_name)], columns=["total_processed", "total_written", "too_short"]
        )
        df.index = pd.MultiIndex.from_tuples([(srx, srr)], names=["srx", "srr"])
        df.to_csv(file_name.replace("_1.trim.clean.fastq.gz.log", ".trim.clean.tsv"), sep="\t")

        if df.total_written[0] < 1000:
            Path(Path(file_name).parent, "ATROPOS_BAD").touch()
    except:
        print(file_name)
        global CNT
        CNT += 1


def parse_atropos(file_name):
    try:
        with open(file_name) as fh:
            string = fh.read().replace(",", "")
        tot_processed = int(re.findall(r"Total read.*processed:\s+(\d+)", string)[0])
        tot_written = int(re.findall(r"[Read|Pair]s written \(passing filters\):\s+(\d+)", string)[0])
        too_short = int(re.findall(r"[Read|Pair]s that were too short:\s+(\d+)", string)[0])
        return tot_processed, tot_written, too_short
    except IndexError:
        pass


def parse_wildcards(file_name):
    srx = re.findall(r"[SED]RX\d+", file_name)[0]
    srr = re.findall(r"[SED]RR\d+", file_name)[0]
    return srx, srr


if __name__ == "__main__":
    CNT = 0
    for file_name in Path("output/prealn-wf/samples").glob("**/**/*.trim.clean.fastq.gz.log"):
        main(file_name.as_posix())
    print(f"{CNT:,} Exceptions")
