__author__ = "Justin Fear"
__copyright__ = "Copyright 2016, Justin Fear"
__email__ = "justin.fear@nih.gov"
__license__ = "MIT"

import sys
from tempfile import NamedTemporaryFile
from snakemake.shell import shell
from lcdblib.snakemake import aligners

sys.path.insert(0, '../lib/python')
from ncbi_remap.snakemake import get_flag

inputs = snakemake.input
outputs = snakemake.output
hisat2_extra = snakemake.params.get('hisat2_extra', '')
samtools_view_extra = snakemake.params.get('samtools_view_extra', '')
samtools_sort_extra = snakemake.params.get('samtools_sort_extra', '')

log = snakemake.log_fmt_shell()

# Look up Layout
flag = get_flag(inputs.layout)

if flag == 'PE':
    in_R1 = inputs.R1
    in_R2 = inputs.R2
    fastqs = '-1 {0} -2 {1}'.format(inputs.R1, inputs.R2)
elif flag == 'keep_R2':
    fastqs = '-U {0}'.format(inputs.R2)
else:
    fastqs = '-U {0}'.format(inputs.R1)

# Look up strand
flag2 = get_flag(inputs.strand)

if (flag == 'PE') & (flag2 == 'first_strand'):
    strand = '--rna-strandness FR'
elif (flag == 'PE') & (flag2 == 'second_strand'):
    strand = '--rna-strandness RF'
elif flag2 == 'first_strand':
    strand = '--rna-strandness F'
elif flag2 == 'second_strand':
    strand = '--rna-strandness R'
else:
    strand = ''

# Grab index information
prefix = aligners.prefix_from_hisat2_index(snakemake.input.index)

# Create temporary files to store intermediates. Will use $TMDPIR if set.
sam = NamedTemporaryFile(suffix='.sam', delete=False).name
bam = NamedTemporaryFile(suffix='.bam', delete=False).name
sort_bam = NamedTemporaryFile(suffix='.sort.bam', delete=False).name

shell(
    "hisat2 "
    "-x {prefix} "
    "{fastqs} "
    "--threads {snakemake.threads} "
    "{hisat2_extra} "
    "{strand} "
    "-S {sam} "
    "{log}"
)

# hisat2 outputs SAM format so we convert to BAM here.
shell(
    "samtools view -Sb "
    "{samtools_view_extra} "
    "{sam} "
    "> {bam} && rm {sam}"
)

# sort the BAM and clean up
shell(
    "samtools sort "
    "-o {sort_bam} "
    "{samtools_sort_extra} "
    "-O BAM "
    "{bam} "
    "&& cp {sort_bam} {outputs.bam} "
    "&& rm {bam} "
    "&& rm {sort_bam} "
)
