import sys
from pathlib import Path
from collections import namedtuple
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_featureCounts_counts

Files = namedtuple("Files", "srx file_name output")


def main():
    workers = Pool(24)
    work = []
    for srx_pth in Path("../output/rnaseq-wf/samples").iterdir():
        srx = srx_pth.name
        if not srx_pth.is_dir():
            continue

        files = Files(
            srx,
            srx_pth / f"{srx}.bam.intergenic.counts",
            f"../output/rnaseq-wf/intergenic_counts/{srx}.parquet"
        )

        if files.file_name.exists():
            work.append(files)

    workers.map(parse_srx, work)


def parse_table(data, srx):
    df = pd.read_table(data, index_col=0, squeeze=True).astype(np.uint32).to_frame().T
    df.index = pd.Index([srx], name="srx")
    return df


def parse_srx(files):
    # Check if file is from feature counts or is already parsed.
    with files.file_name.open() as fh:
        data = fh.read(10)

    if data.startswith("#"):
        df = parse_featureCounts_counts(files.file_name, files.srx)
    else:
        df = parse_table(files.file_name, files.srx)

    df.to_parquet(files.output)
    files.file_name.unlink()


if __name__ == "__main__":
    main()
