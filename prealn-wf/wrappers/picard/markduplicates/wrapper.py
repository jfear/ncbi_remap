from snakemake.shell import shell
log = snakemake.log_fmt_shell()
extra = snakemake.params.get('extra', '')
java_args = snakemake.params.get('java_args', '')
mem = snakemake.resources.get('mem_gb', '')

if mem:
    mem = '-Xmx{}g'.format(mem)

shell(
    'picard '
    '{mem} '
    '{java_args} '
    'MarkDuplicates '
    'INPUT={snakemake.input.bam} '
    'OUTPUT={snakemake.output.bam} '
    'METRICS_FILE={snakemake.output.metrics} '
    '{extra} '
    '{log} '
)
