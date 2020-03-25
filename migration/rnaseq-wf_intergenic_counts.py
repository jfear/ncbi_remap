import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_featureCounts_counts


def main():
    for srx_pth in Path("../output/rnaseq-wf/samples").iterdir():
        if not srx_pth.is_dir():
            continue

        srx = srx_pth.name
        file_name = srx_pth / f"{srx}.bam.intergenic.counts"
        output = f"../output/rnaseq-wf/intergenic_counts/{srx}.parquet"

        if not file_name.exists():
            continue

        try:
            # Check if file is from feature counts or is already parsed.
            with file_name.open() as fh:
                data = fh.read(10)

            if data.startswith("#"):
                df = parse_featureCounts_counts(file_name, srx)
            else:
                df = parse_table(file_name, srx)

            df.to_parquet(output)
            file_name.unlink()
        except Exception as e:
            print(data[:10])
            print(file_name)
            raise e


def parse_table(data, srx):
    df = pd.read_table(data, index_col=0, squeeze=True).astype(np.uint32).to_frame().T
    df.index = pd.Index([srx], name="srx")
    return df


if __name__ == "__main__":
    main()
