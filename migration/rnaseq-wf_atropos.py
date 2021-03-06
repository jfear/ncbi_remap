import re
from collections import namedtuple
from pathlib import Path

import numpy as np
import pandas as pd

RNASEQ_PATH = Path("../output/rnaseq-wf/samples")
DTYPES = {
    "total_processed": np.int64,
    "total_written": np.int64,
    "too_short": np.int64,
}

Files = namedtuple("Files", "srr idx table log output")


class AtroposException(Exception):
    """Basic exception when there are problems running Atropos"""


def main():
    for srr_pth in RNASEQ_PATH.glob("**/SRR*"):
        if not srr_pth.is_dir():
            continue
        srr = srr_pth.name
        files = Files(
            srr,
            pd.Index([srr], name="srr"),
            (srr_pth / f"{srr}.trim.clean.tsv"),
            (srr_pth / f"{srr}_1.trim.clean.fastq.gz.log"),
            f"../output/rnaseq-wf/atropos/{srr}.parquet",
        )

        if files.table.exists():
            convert_tsv(files)
        elif files.log.exists():
            convert_log(files)


def convert_tsv(files):
    df = pd.read_table(files.table)[DTYPES.keys()].fillna(0).astype(DTYPES)
    df.index = files.idx
    df.to_parquet(files.output)


def convert_log(files):
    try:
        df = pd.DataFrame([parse_atropos(files.log)], index=files.idx, columns=DTYPES.keys())
        df.to_parquet(files.output)

        if df.total_written[0] < 1000:
            raise AtroposException(
                "Too few reads after Atropos: %d (%s)", df.total_written[0], files.log
            )

    except AtroposException:
        print(f"Problems with Atropos: {files.log}")
        atropos_bad_path = Path(f"../output/rnaseq-wf/atropos_bad")
        atropos_bad_path.mkdir(exist_ok=True)
        (atropos_bad_path / files.srr).touch()


def parse_atropos(file_name):
    with open(file_name) as fh:
        string = fh.read().replace(",", "")

    if "ERROR" in string:
        raise AtroposException("Problem with Atropos: %s", file_name)

    try:
        tot_processed = int(re.findall(r"Total read.*processed:\s+(\d+)", string)[0])
        tot_written = int(
            re.findall(r"[Read|Pair]s written \(passing filters\):\s+(\d+)", string)[0]
        )
        too_short = int(re.findall(r"[Read|Pair]s that were too short:\s+(\d+)", string)[0])
        return tot_processed, tot_written, too_short
    except IndexError:
        raise AtroposException("Problem with Atropos: %s", file_name)


if __name__ == "__main__":
    main()
