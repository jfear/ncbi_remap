import os
from subprocess import SubprocessError
from snakemake.shell import shell
import tempfile

# Pull in parameters
extra = snakemake.params.get("extra", "")
aligner = snakemake.params.get("aligner", "bowtie2")
subset = snakemake.params.get("subset", 100000)

log = snakemake.log_fmt_shell()

# snakemake.params.fastq_screen_config can be either a dict or a string. If
# string, interpret as a filename pointing to the fastq_screen config file.
# Otherwise, create a new tempfile out of the contents of the dict:

tmp = tempfile.NamedTemporaryFile(delete=False).name
with open(tmp, "w") as fout:
    for k, v in snakemake.input.items():
        if k == "layout":
            # Only included layout to make sure download checks were run.
            continue

        if k == "fastq":
            # skip because not reference
            continue

        label = k
        index = ".".join(v.replace(".rev", "").split(".")[:-2])
        fout.write("\t".join(["DATABASE", label, index, aligner.upper()]) + "\n")

    config_file = tmp

# fastq_screen hard-codes filenames according to this prefix. We will send
# hard-coded output to a temp dir, and then move them later.
prefix = os.path.basename(snakemake.input.fastq.split(".fastq")[0])
tempdir = tempfile.mkdtemp()

shell(
    "fastq_screen --outdir {tempdir} "
    "--force "
    "--aligner {aligner} "
    "--conf {config_file} "
    "--subset {subset} "
    "--threads {snakemake.threads} "
    "{extra} "
    "{snakemake.input.fastq} "
    "{log}"
)

# Make sure processing completed
with open(snakemake.log[0], "r") as fh:
    if not "Processing complete" in fh.read():
        raise SubprocessError(
            f"Fastq Screen log missing 'Processing complete': {snakemake.wildcards.srx}/{snakemake.wildcards.srr}"
        )

# Move output to the filenames specified by the rule
shell("cp {tempdir}/{prefix}_screen.txt {snakemake.output.txt}")

# Clean up temp
shell("rm -r {tempdir}")
shell("rm {tmp}")
