from subprocess import SubprocessError

from snakemake.shell import shell

log = snakemake.log_fmt_shell()
extra = snakemake.params.get('extra', '')
java_args = snakemake.params.get('java_args', '')
mem = int(snakemake.resources.get('mem_gb', ''))

if mem >= 10:
    mem = '-Xmx{}g'.format(mem - 2)
else:
    mem = '-Xmx{}g'.format(2)

shell(
    'picard '
    '{mem} '
    '{java_args} '
    'CollectRnaSeqMetrics '
    'REF_FLAT={snakemake.input.refflat} '
    'INPUT={snakemake.input.bam} '
    'OUTPUT={snakemake.output[0]} '
    '{extra} '
    '{log} '
)

with open(snakemake.log[0], "r") as fh:
    if not "CollectRnaSeqMetrics done" in fh.read():
        raise SubprocessError(f"CollectRnaSeqMetrics log not complete: {snakemake.wildcards.srx}/{snakemake.wildcards.srr}")