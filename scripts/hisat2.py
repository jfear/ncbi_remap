from tempfile import NamedTemporaryFile

from snakemake.shell import shell
import pandas as pd

inputs = snakemake.input
outputs = snakemake.output
hisat2_extra = snakemake.params.get("hisat2_extra", "")
samtools_view_extra = snakemake.params.get("samtools_view_extra", "")
samtools_sort_extra = snakemake.params.get("samtools_sort_extra", "")

log = snakemake.log_fmt_shell()

# Look up Layout
layout = pd.read_parquet(inputs.layout).layout[0]
if layout == "PE":
    fastqs = "-1 {0} -2 {1}".format(inputs.R1, inputs.R2)
elif layout == "keep_R2":
    fastqs = "-U {0}".format(inputs.R2)
else:
    fastqs = "-U {0}".format(inputs.R1)

# Look up strand
strand_param = ""
if inputs.get("strand"):
    strand = pd.read_parquet(inputs.strand).strand[0]
    if (layout == "PE") & (strand == "first_strand"):
        strand_param = "--rna-strandness FR"
    elif (layout == "PE") & (strand == "second_strand"):
        strand_param = "--rna-strandness RF"
    elif strand == "first_strand":
        strand_param = "--rna-strandness F"
    elif strand == "second_strand":
        strand_param = "--rna-strandness R"


# Grab index information
prefix = '.'.join(inputs.index.split('.')[:-2])

# Create temporary files to store intermediates. Will use $TMDPIR if set.
sam = NamedTemporaryFile(suffix=".sam", delete=False).name
bam = NamedTemporaryFile(suffix=".bam", delete=False).name
sort_bam = NamedTemporaryFile(suffix=".sort.bam", delete=False).name

shell(
    "hisat2 "
    "-x {prefix} "
    "{fastqs} "
    "--threads {snakemake.threads} "
    "{hisat2_extra} "
    "{strand_param} "
    "-S {sam} "
    "{log}"
)

# hisat2 outputs SAM format so we convert to BAM here.
shell("samtools view -Sb {samtools_view_extra} {sam} > {bam} && rm {sam}")

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

# Make index
shell(
    "samtools index "
    "{outputs.bam} "
)
