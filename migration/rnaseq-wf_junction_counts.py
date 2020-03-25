import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_featureCounts_jcounts


def main():
    for srx_pth in Path("../output/rnaseq-wf/samples").iterdir():
        if not srx_pth.is_dir():
            continue

        srx = srx_pth.name
        file_name = srx_pth / f"{srx}.bam.counts.jcounts"
        output = f"../output/rnaseq-wf/junction_counts/{srx}.parquet"

        if not file_name.exists():
            continue

        try:
            df = parse_featureCounts_jcounts(file_name, srx)
            df.to_parquet(output)
            file_name.unlink()
        except Exception as e:
            print(file_name)
            raise e


if __name__ == "__main__":
    main()
