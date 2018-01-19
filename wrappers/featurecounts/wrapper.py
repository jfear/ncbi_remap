__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

import sys
from snakemake.shell import shell

sys.path.insert(0, '../lib')
from ncbi_remap.snakemake import get_flag

inputs = snakemake.input
outputs = snakemake.output
params = snakemake.params
log = snakemake.log_fmt_shell()

# Look up Layout
flag = get_flag(inputs.layout)
if flag == 'PE':
    extra = params.extra_pe
elif flag == 'keep_R2':
    extra = params.extra_se
else:
    extra = params.extra_se

# Look up strand
try:
    strand = get_flag(inputs.strand)
except:
    strand = 'unstranded'

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
