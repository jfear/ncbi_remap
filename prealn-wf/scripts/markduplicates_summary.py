import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_picard_markduplicate_metrics

dtypes = {
    "UNPAIRED_READS_EXAMINED": np.int64,
    "READ_PAIRS_EXAMINED": np.int64,
    "UNPAIRED_READ_DUPLICATES": np.int64,
    "READ_PAIR_DUPLICATES": np.int64,
    "PERCENT_DUPLICATION": np.float64,
    "ESTIMATED_LIBRARY_SIZE": np.int64,
}

df = parse_picard_markduplicate_metrics(snakemake.input[0])[dtypes.keys()].fillna(0).astype(dtypes)
df.PERCENT_DUPLICATION = df.PERCENT_DUPLICATION * 100

df.columns = [col.lower() for col in df.columns]
df.index = pd.Index([snakemake.wildcards.srr], name="srr")

df.to_parquet(snakemake.output[0])
