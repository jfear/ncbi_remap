import sys

sys.path.insert(0, "../src")
from ncbi_remap.parser import parse_featureCounts_jcounts

df = parse_featureCounts_jcounts(snakemake.input[0], snakemake.wildcards.srx)
df.to_parquet(snakemake.output[0])
