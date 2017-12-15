__author__ = "Justin Fear"
__copyright__ = "Copyright 2016, Justin Fear"
__email__ = "justin.fear@nih.gov"
__license__ = "MIT"

"""Example usage.

rule fastq_dump:
    output:
        fq1=patterns['fastq']['r1'],
        fq2=patterns['fastq']['r2'],
        flag=patterns['layout'],
        summary=patterns['fastq']['summary'],
    log: patterns['fastq']['r1'] + '.log'
    wrapper:
        wrapper_for('wrappers/fastq_dump')

"""

import os
import sys
import gzip
import shutil as sh
from pathlib import Path

import numpy as np
import pandas as pd

from snakemake.shell import shell

sys.path.insert(0, '../lib')
from ncbi_remap.fastq import check_fastq, md5sum, fastq_stats
from ncbi_remap.snakemake import put_flag

# Get TMPDIR
TMPDIR = os.environ['TMPDIR']

log = snakemake.log_fmt_shell()
output = snakemake.output
wildcards = snakemake.wildcards

# Get current sample id
sample = wildcards.srr

# Dump FASTQ to a TMPDIR
shell("fastq-dump -O {TMPDIR} -M 0 --split-files {sample} {log}")

# Summarize R1
t1 = os.path.join(TMPDIR, sample + '_1.fastq')
R1md5 = md5sum(t1)
R1Libsize, R1avgLen = fastq_stats(t1)

with open(t1, 'rb') as f_in:
    with gzip.open(output.fq1, 'wb') as f_out:
        sh.copyfileobj(f_in, f_out)

# Summarize R2 if it exists
t2 = os.path.join(TMPDIR, sample + '_2.fastq')
if check_fastq(t2):
    # Get md5sum
    R2md5 = md5sum(t2)
    R2Libsize, R2avgLen = fastq_stats(t2)

    with open(t2, 'rb') as f_in:
        with gzip.open(output.fq2, 'wb') as f_out:
            sh.copyfileobj(f_in, f_out)
else:
    R2md5 = R2Libsize = R2avgLen = None
    Path(output.fq2).touch()

# Save summaries
df = pd.DataFrame([[R1md5, R1Libsize, R1avgLen, R2md5, R2Libsize, R2avgLen]],
                  columns=['md5_R1', 'libsize_R1', 'avgLen_R1', 'md5_R2', 'libsize_R2', 'avgLen_R2'])

df.to_csv(output.summary, sep='\t', index=False)

# Figure out flags
if (R1Libsize > 1000) & (R2Libsize is None) & (R1avgLen > 10) & (R2avgLen is None):
    # R2 does not exists so SE
    put_flag(output.flag, 'SE')

elif (R2Libsize is None)  & (R2avgLen is None):
    # SE but R1 looks bad, flag as download bad
    fname = os.path.join(os.path.dirname(output.fq1), 'DOWNLOAD_BAD')
    Path(fname).touch()

elif (R1Libsize > 1000) & (R2Libsize > 1000) & (R1avgLen > 10) & (R2avgLen > 10):
    if R1Libsize == R2Libsize:
        # Both reads look good so PE
        put_flag(output.flag, 'PE')
    else:
        # There are an uneven number of reads between R1 and R2.
        # Instead of messing with this, just consider SE and use R1.
        put_flag(output.flag, 'keep_R1')

elif (R1Libsize > 1000) & (R1avgLen > 10):
    # Only R1 looks ok, so SE
    put_flag(output.flag, 'keep_R1')

elif (R2Libsize > 1000) & (R2avgLen > 10):
    # Only R2 looks ok, consider single-end
    put_flag(output.flag, 'keep_R2')

else:
    # R1 and R2 look bad, flag as download bad
    fname = os.path.join(os.path.dirname(output.fq1), 'DOWNLOAD_BAD')
    Path(fname).touch()
