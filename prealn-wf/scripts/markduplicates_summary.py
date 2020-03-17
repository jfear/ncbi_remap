import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_picard_markduplicate_metrics

df = parse_picard_markduplicate_metrics(snakemake.input[0])[
    [
        "UNPAIRED_READS_EXAMINED",
        "READ_PAIRS_EXAMINED",
        "UNPAIRED_READ_DUPLICATES",
        "READ_PAIR_DUPLICATES",
        "READ_PAIR_OPTICAL_DUPLICATES",
        "PERCENT_DUPLICATION",
        "ESTIMATED_LIBRARY_SIZE",
    ]
].fillna(0)

# Fillna sometimes infers the wrong data type.
dtypes = {
    "UNPAIRED_READS_EXAMINED": np.int64,
    "READ_PAIRS_EXAMINED": np.int64,
    "UNPAIRED_READ_DUPLICATES": np.int64,
    "READ_PAIR_DUPLICATES": np.int64,
    "READ_PAIR_OPTICAL_DUPLICATES": np.int64,
    "PERCENT_DUPLICATION": np.float64,
    "ESTIMATED_LIBRARY_SIZE": np.int64,
}
df = df.astype(dtypes)

df.columns = [col.lower() for col in df.columns]
df.index = [snakemake.wildcards.srr]
df.to_parquet(snakemake.output[0])
