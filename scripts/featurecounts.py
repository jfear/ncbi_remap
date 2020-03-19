import sys
from subprocess import SubprocessError

from snakemake.shell import shell
import pandas as pd

sys.path.insert(0, '../src')
from ncbi_remap.snakemake import get_flag

inputs = snakemake.input
outputs = snakemake.output
params = snakemake.params
log = snakemake.log_fmt_shell()

# Look up Layout
flag = pd.read_parquet(inputs.layout).layout[0]
if flag == 'PE':
    extra = params.extra_pe
else:
    extra = params.extra_se

# Look up strand
strand = get_flag(inputs.strand)
if strand == 'first_strand':
    extra += '-s 1'
elif strand == 'second_strand':
    extra += '-s 2'
else:
    extra += '-s 0'

shell(
    "featureCounts "
    "-T {snakemake.threads} "
    "{extra} "
    "-a {inputs.annotation} "
    "-o {outputs.counts} "
    "{inputs.bam} "
    "{log} "
)

with open(snakemake.log[0], "r") as fh:
    if not "Summary of counting results" in fh.read():
        raise SubprocessError(f"FeatureCounts log not complete: {snakemake.wildcards.srx}/{snakemake.wildcards.srr}")