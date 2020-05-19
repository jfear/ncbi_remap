import sys
import os

sys.path.insert(0, os.getenv("SRC_PATH", "../src"))
from ncbi_remap.parser import parse_featureCounts_counts

df = parse_featureCounts_counts(snakemake.input[0], snakemake.wildcards.srx)
df.to_parquet(snakemake.output[0])
