__author__ = "Justin Fear"
__copyright__ = "Copyright 2016, Justin Fear"
__email__ = "justin.fear@nih.gov"
__license__ = "MIT"

import sys
from snakemake.shell import shell

sys.path.insert(0, '../lib')
from ncbi_remap.snakemake import get_flag

log = snakemake.log_fmt_shell()
inputs = snakemake.input
outputs = snakemake.output
params = snakemake.params

# Look up Layout
flag = get_flag(inputs.layout)
if flag == 'PE':
    in_R1 = inputs.R1
    in_R2 = inputs.R2
    extra = params.extra_pe
elif flag == 'keep_R2':
    in_R1 = None
    in_R2 = inputs.R2
    extra = params.extra_se
else:
    in_R1 = inputs.R1
    in_R2 = None
    extra = params.extra_se

# Get outputs
out_R1 = outputs.R1
out_R2 = outputs.R2


if (in_R1 is not None) & (in_R2 is not None):
    # Run paired end if both in_R1 and in_R2 are provided
    shell(
        "atropos trim "
        "--threads {snakemake.threads} "
        "{extra} "
        "-pe1 {in_R1} "
        "-pe2 {in_R2} "
        "-o {out_R1} "
        "-p {out_R2} "
        "{log}"
    )

elif (in_R1 is not None):
    # Run SE if only in_R1
    shell(
        "atropos trim "
        "{extra} "
        "--threads {snakemake.threads} "
        "-se {in_R1} "
        "-o {out_R1} "
        "{log} "
        "&& touch {out_R2}"
    )

elif (in_R2 is not None):
    # Run SE if only in_R2
    shell(
        "atropos trim "
        "{extra} "
        "--threads {snakemake.threads} "
        "-se {in_R2} "
        "-o {out_R2} "
        "{log} "
        "&& touch {out_R1}"
    )
