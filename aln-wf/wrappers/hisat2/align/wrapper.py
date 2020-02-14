__author__ = "Justin Fear"
__copyright__ = "Copyright 2016, Justin Fear"
__email__ = "justin.fear@nih.gov"
__license__ = "MIT"

from tempfile import NamedTemporaryFile
from snakemake.shell import shell
from lcdblib.snakemake import aligners

from ncbi_remap.snakemake import get_flag

r1 = snakemake.input.r1
r2 = snakemake.input.r2
hisat2_extra = snakemake.params.get("hisat2_extra", "")
samtools_view_extra = snakemake.params.get("samtools_view_extra", "")
samtools_sort_extra = snakemake.params.get("samtools_sort_extra", "")

log = snakemake.log_fmt_shell()

# Look up Layout
layout = get_flag(snakemake.input.layout)

if layout == "PE":
    fastqs = "-1 {0} -2 {1}".format(r1, r2)
elif layout == "keep_R2":
    fastqs = "-U {0}".format(r2)
else:
    fastqs = "-U {0}".format(r1)

# Look up strand
strand = get_flag(snakemake.input.strand)

if (layout == "PE") & (strand == "first_strand"):
    strand_param = "--rna-strandness FR"
elif (layout == "PE") & (strand == "second_strand"):
    strand_param = "--rna-strandness RF"
elif strand == "first_strand":
    strand_param = "--rna-strandness F"
elif strand == "second_strand":
    strand_param = "--rna-strandness R"
else:
    strand_param = ""

# Grab index information
prefix = aligners.prefix_from_hisat2_index(snakemake.input.index)

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
    "&& cp {sort_bam} {snakemake.output.bam} "
    "&& rm {bam} "
    "&& rm {sort_bam} "
    "&& samtools index {snakemake.output.bam} "
)
